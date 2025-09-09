import io
import os
import csv
import json
import time
import logging
import threading
from typing import Dict, List, Tuple, Optional, Any
from functools import lru_cache
from datetime import datetime, timedelta

import requests
from flask import Flask, request, jsonify, render_template
from werkzeug.exceptions import RequestEntityTooLarge
from werkzeug.utils import secure_filename


# -----------------------------------------------------------------------------
# App configuration
# -----------------------------------------------------------------------------
app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 20 * 1024 * 1024  # 20MB max file size
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'dev-key-change-in-production')

BASE_URL = 'https://www.subtiwiki.uni-goettingen.de/v5/api'
REQUEST_TIMEOUT_SEC = 20
MAX_RETRIES = 3
RETRY_BACKOFF_SEC = 1
ALLOWED_EXTENSIONS = {'csv', 'tsv', 'txt'}
MAX_ROWS = 100000


# -----------------------------------------------------------------------------
# Logging setup (file-based)
# -----------------------------------------------------------------------------
LOGS_DIR = os.path.join(os.path.dirname(__file__), 'logs')
os.makedirs(LOGS_DIR, exist_ok=True)
from logging.handlers import RotatingFileHandler

log_handler = RotatingFileHandler(
    os.path.join(LOGS_DIR, 'app.log'),
    maxBytes=10 * 1024 * 1024,
    backupCount=5
)
log_handler.setFormatter(logging.Formatter(
    '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
))
log_handler.setLevel(logging.DEBUG)
logging.getLogger().handlers = []
logging.getLogger().addHandler(log_handler)
logging.getLogger().setLevel(logging.DEBUG)
logging.info('SubtiWiki app startup')

# Simple TTL cache
class TTLCache:
    def __init__(self, ttl_seconds: int = 300):
        self._cache: Dict[str, Tuple[Any, datetime]] = {}
        self._lock = threading.Lock()
        self.ttl = ttl_seconds

    def get(self, key: str) -> Optional[Any]:
        with self._lock:
            if key in self._cache:
                value, expiry = self._cache[key]
                if datetime.now() < expiry:
                    return value
                del self._cache[key]
        return None

    def set(self, key: str, value: Any):
        with self._lock:
            self._cache[key] = (value, datetime.now() + timedelta(seconds=self.ttl))

graph_cache = TTLCache(ttl_seconds=900)
gene_cache = TTLCache(ttl_seconds=1800)


# -----------------------------------------------------------------------------
# Gene synonyms loading (Tools/static/api_cache/gene_synonyms.json)
# -----------------------------------------------------------------------------
_synonyms_lock = threading.Lock()
_synonyms_cache: Optional[Dict] = None


def _load_gene_synonyms() -> Dict:
    global _synonyms_cache
    with _synonyms_lock:
        if _synonyms_cache is not None:
            return _synonyms_cache

        # Resolve path to local file in this project folder
        # Current file is Tools/B.Sub_Graph/app.py → load B.Sub_Graph/gene_synonyms.json
        synonyms_path = os.path.join(os.path.dirname(__file__), 'gene_synonyms.json')
        if not os.path.exists(synonyms_path):
            logging.error(f'Missing gene_synonyms.json at {synonyms_path}')
            return {'synonyms': {}, 'reverse_synonyms': {}}

        try:
            with open(synonyms_path, 'r', encoding='utf-8') as f:
                _synonyms_cache = json.load(f)
            logging.info(f'Loaded {len(_synonyms_cache.get("synonyms", {}))} gene synonyms')
            return _synonyms_cache
        except Exception as e:
            logging.error(f'Error loading synonyms: {e}')
            return {'synonyms': {}, 'reverse_synonyms': {}}


def normalize_gene_name(raw_name: str) -> str:
    """Resolve a gene name to its primary SubtiWiki name using synonyms cache.

    - Tries exact key in 'synonyms' dict
    - Tries lowercase key in 'reverse_synonyms' dict
    - Falls back to original name if no mapping is found
    """
    logging.debug(f"Normalizing gene name for: '{raw_name}'")
    if not raw_name:
        logging.debug("Received empty raw_name, returning as is.")
        return raw_name

    data = _load_gene_synonyms()
    synonyms = data.get('synonyms', {})
    reverse = data.get('reverse_synonyms', {})

    # Direct hit: already primary
    if raw_name in synonyms:
        logging.debug(f"'{raw_name}' is already a primary name.")
        return raw_name

    # Reverse lookup using lowercase
    lower = raw_name.lower()
    primary = reverse.get(lower)
    if primary:
        logging.debug(f"Found primary name '{primary}' for synonym '{raw_name}'.")
        return primary
    else:
        logging.debug(f"No primary name found for '{raw_name}'. Returning original name.")
        return raw_name


# -----------------------------------------------------------------------------
# HTTP helpers
# -----------------------------------------------------------------------------
def success(data):
    resp = jsonify({'code': 200, 'is_success': True, 'data': data, 'message': None})
    resp.headers['Access-Control-Allow-Origin'] = '*'
    return resp


def error_response(status_code: int, type_: str, message: str, citation: str, recommended: str):
    payload = {
        'type': type_,
        'message': message,
        'citation': citation,
        'recommended_next_step': recommended
    }
    resp = jsonify({'code': status_code, 'is_success': False, 'data': None, 'message': payload})
    resp.status_code = status_code
    resp.headers['Access-Control-Allow-Origin'] = '*'
    return resp


# -----------------------------------------------------------------------------
# SubtiWiki API clients
# -----------------------------------------------------------------------------
def _retry_request(method, url, **kwargs):
    for attempt in range(MAX_RETRIES):
        logging.debug(f"Attempt {attempt + 1} of {MAX_RETRIES} to fetch {url}")
        try:
            r = requests.request(method, url, timeout=REQUEST_TIMEOUT_SEC, **kwargs)
            if r.status_code == 200:
                logging.debug(f"Successfully fetched {url}")
                return r
            if r.status_code == 429 or r.status_code >= 500:
                if attempt < MAX_RETRIES - 1:
                    sleep_time = RETRY_BACKOFF_SEC * (2 ** attempt)
                    logging.warning(f"Received status {r.status_code} from {url}. Retrying in {sleep_time}s...")
                    time.sleep(sleep_time)
                    continue
            r.raise_for_status()
        except requests.RequestException as e:
            logging.error(f"Request to {url} failed: {e}")
            if attempt == MAX_RETRIES - 1:
                raise
            sleep_time = RETRY_BACKOFF_SEC * (2 ** attempt)
            logging.warning(f"Retrying in {sleep_time}s...")
            time.sleep(sleep_time)

def fetch_full_interaction_graph() -> Dict:
    logging.debug("Fetching full interaction graph...")
    cached = graph_cache.get('full_graph')
    if cached:
        logging.debug("Returning cached full interaction graph.")
        return cached
    url = f"{BASE_URL}/interaction/graph"
    r = _retry_request('GET', url)
    body = r.json()
    if not isinstance(body, dict) or 'data' not in body:
        raise ValueError(json.dumps({
            'type': 'API_SCHEMA_ERROR',
            'message': 'Missing data field in response',
            'citation': 'Subtiwiki API docs',
            'recommended_next_step': 'Check API schema changes.'
        }))
    graph_cache.set('full_graph', body['data'])
    logging.debug("Successfully fetched and cached the full interaction graph.")
    return body['data']


def fetch_gene_id_by_name(primary_name: str) -> Optional[int]:
    logging.debug(f"Fetching gene ID for: {primary_name}")
    cache_key = f'gene:{primary_name}'
    cached = gene_cache.get(cache_key)
    if cached is not None:
        logging.debug(f"Returning cached gene ID for {primary_name}")
        return cached
    params = {'q': primary_name, 'category': 'Gene', 'mode': 'exact'}
    url = f"{BASE_URL}/search/"
    try:
        r = _retry_request('GET', url, params=params)
        body = r.json()
        exact = body.get('exact_hits') or []
        for hit in exact:
            if hit.get('category') == 'Gene' and str(hit.get('name')) == primary_name:
                gene_id = hit.get('id')
                logging.debug(f"Found gene ID {gene_id} for {primary_name}")
                gene_cache.set(cache_key, gene_id)
                return gene_id
        logging.debug(f"No exact match found for gene: {primary_name}")
        gene_cache.set(cache_key, None)
        return None
    except Exception as e:
        logging.warning(f'Gene ID lookup failed for {primary_name}: {e}')
        return None


# -----------------------------------------------------------------------------
# File parsing
# -----------------------------------------------------------------------------
def sniff_delimiter(sample: str) -> str:
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",	;|")
        return dialect.delimiter
    except Exception:
        # Default to comma
        return ','


def parse_tabular(file_storage) -> Tuple[List[str], List[List[str]]]:
    """Parse CSV/TSV-like files. Returns (headers, rows)."""
    logging.debug("Parsing tabular file...")
    content = file_storage.read()
    file_storage.seek(0)
    text = content.decode('utf-8', errors='replace')
    sample = '\n'.join(text.splitlines()[:5])
    delim = sniff_delimiter(sample)
    reader = csv.reader(io.StringIO(text), delimiter=delim)
    rows = list(reader)
    if not rows:
        logging.error("File is empty.")
        raise ValueError('Empty file')
    headers = [h.strip() for h in rows[0]]
    data_rows = [row for row in rows[1:] if any(col.strip() for col in row)]
    logging.debug(f"Parsed {len(headers)} headers and {len(data_rows)} data rows.")
    return headers, data_rows


# -----------------------------------------------------------------------------
# Endpoints
# -----------------------------------------------------------------------------
@app.after_request
def add_cors_headers(resp):
    resp.headers['Access-Control-Allow-Origin'] = '*'
    resp.headers['Access-Control-Allow-Headers'] = 'Content-Type'
    resp.headers['Access-Control-Allow-Methods'] = 'GET,POST,OPTIONS'
    return resp


@app.route('/api/health', methods=['GET'])
def health():
    return success({'status': 'ok'})

@app.route('/', methods=['GET'])
def ui_index():
    return render_template('index.html')

@app.route('/api/full_graph', methods=['GET'])
def full_graph():
    try:
        t0 = time.time()
        data = fetch_full_interaction_graph()
        molecules = data.get('molecules', [])
        interactions = data.get('interactions', [])
        logging.info('Fetched full graph: %d molecules, %d interactions in %.2fs',
                     len(molecules), len(interactions), time.time() - t0)
        return success({'molecules': molecules, 'interactions': interactions})
    except Exception as e:
        try:
            info = json.loads(str(e))
        except Exception:
            info = {
                'type': 'UNEXPECTED_ERROR',
                'message': str(e),
                'citation': 'Subtiwiki API documentation',
                'recommended_next_step': 'Check server logs and retry.'
            }
        logging.exception('full_graph error')
        return error_response(502, info['type'], info['message'], info['citation'], info['recommended_next_step'])


@app.route('/api/process_upload', methods=['POST'])
def process_upload():
    """Process uploaded file and resolve gene names → SubtiWiki gene IDs in parallel.

    Form fields:
      - file: CSV/TSV file
      - gene_col (optional): exact header for gene names
      - log2fc_col (optional): exact header for log2FC
      - pval_col (optional): exact header for p-value
    """
    try:
        if 'file' not in request.files:
            return error_response(400, 'INPUT_ERROR', 'Missing file in form-data', 'App input contract', 'Attach a CSV/TSV file as form field "file".')

        file_storage = request.files['file']
        filename = (file_storage.filename or '').lower()
        if not filename.endswith(('.csv', '.tsv', '.txt')):
            return error_response(415, 'UNSUPPORTED_FORMAT', 'Only CSV/TSV/TXT are supported at this time', 'File format policy', 'Export your data as CSV or TSV and upload again.')

        headers, rows = parse_tabular(file_storage)
        if not headers or not rows:
            return error_response(400, 'INPUT_ERROR', 'The file is empty or has no data rows', 'File parsing', 'Provide a non-empty CSV/TSV with a header row.')

        # Column selection
        gene_col = request.form.get('gene_col')
        log2fc_col = request.form.get('log2fc_col')
        pval_col = request.form.get('pval_col')

        # Auto-detect if not provided
        lower_headers = [h.lower() for h in headers]
        def guess(colnames: List[str]) -> Optional[str]:
            for cand in colnames:
                for h in headers:
                    if cand in h.lower():
                        return h
            return None

        if not gene_col:
            gene_col = guess(['gene', 'name', 'id'])
        if not log2fc_col:
            log2fc_col = guess(['log2fc', 'log2', 'fold'])
        # pval is optional in pipeline
        if not pval_col:
            pval_col = guess(['pval', 'p_value', 'p-value', 'padj', 'fdr'])

        if not gene_col or not log2fc_col:
            return error_response(400, 'COLUMN_DETECTION_ERROR', 'Could not detect gene or log2FC columns', 'File header requirements', 'Specify gene_col and log2fc_col explicitly in the request.')

        gene_idx = headers.index(gene_col)
        log2fc_idx = headers.index(log2fc_col)
        pval_idx = headers.index(pval_col) if pval_col in headers else None

        # Build raw records
        records: List[Tuple[str, Optional[float], Optional[float]]] = []
        for row in rows:
            try:
                gene_raw = (row[gene_idx] if gene_idx < len(row) else '').strip()
                if not gene_raw:
                    continue
                log2fc_val = row[log2fc_idx] if log2fc_idx < len(row) else ''
                log2fc = float(log2fc_val) if str(log2fc_val).strip() != '' else None
                pval = None
                if pval_idx is not None and pval_idx < len(row):
                    pval_str = str(row[pval_idx]).strip()
                    if pval_str not in ('', 'NA', 'NaN', 'nan', 'null', 'None'):
                        pval = float(pval_str)
                records.append((gene_raw, log2fc, pval))
            except Exception:
                # Skip malformed rows; continue processing
                continue
                
        if not records:
            return error_response(400, 'INPUT_ERROR', 'No valid data rows after parsing', 'File content validation', 'Verify the file has valid gene and log2FC values.')

        # Resolve primary names using synonyms
        resolved_primary: Dict[str, str] = {}
        for gene_raw, _, _ in records:
            resolved_primary[gene_raw] = normalize_gene_name(gene_raw)

        # Parallel gene ID resolution
        from concurrent.futures import ThreadPoolExecutor, as_completed
        unique_primary = sorted(set(resolved_primary.values()))
        id_map: Dict[str, Optional[int]] = {}
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=min(32, max(4, os.cpu_count() or 4))) as ex:
            futures = {ex.submit(fetch_gene_id_by_name, name): name for name in unique_primary}
            for fut in as_completed(futures):
                name = futures[fut]
                try:
                    gid = fut.result()
                except Exception:
                    gid = None
                id_map[name] = gid
        logging.info('Resolved %d gene IDs in %.2fs', len(id_map), time.time() - t0)

        # Build outputs
        data_map: Dict[str, float] = {}
        not_found: List[str] = []
        resolved_list: List[Dict] = []

        for gene_raw, log2fc, pval in records:
            primary = resolved_primary[gene_raw]
            gid = id_map.get(primary)
            resolved_list.append({'inputName': gene_raw, 'primaryName': primary, 'geneId': gid})
            if gid is None:
                not_found.append(gene_raw)
            if log2fc is not None:
                data_map[primary] = log2fc

        return success({
            'headers': headers,
            'selected_columns': {'gene': gene_col, 'log2fc': log2fc_col, 'pval': pval_col},
            'resolvedGenes': resolved_list,
            'notFoundGenes': sorted(set(not_found)),
            'dataMap': data_map
        })

    except Exception as e:
        logging.exception('process_upload error')
        try:
            info = json.loads(str(e))
            return error_response(400, info['type'], info['message'], info['citation'], info['recommended_next_step'])
        except Exception:
            return error_response(500, 'UNEXPECTED_ERROR', str(e), 'Backend processing', 'Check server logs and input file format.')


if __name__ == '__main__':
    # Get the logger instance
    logger = logging.getLogger()
    logger.info("Starting Flask development server on http://127.0.0.1:5000")
    # Setting debug=False is important to ensure our custom logging configuration is used.
    # Flask's debug mode has its own logger that can override our settings.
    app.run(debug=False, host='0.0.0.0', port=5000)