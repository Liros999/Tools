"""
Selection API - Handles gene and group selection for plot coloring
and provides local-cache derived pathway/complex data.
"""

import json
import os
import requests
from typing import Dict, List, Set, Tuple, Optional
from flask import current_app
from src.api.subtiwiki_api import (
    get_genes_in_category_by_name, 
    get_genes_in_category,
    get_all_categories,
    DETAILED_CATEGORIES,
    GENE_CATEGORY_MAP,
    GENE_NAME_TO_ID
)

# Cache for regulations and regulons
REGULATIONS_CACHE = {}
REGULONS_CACHE = {}

# In-process memoization for new features (acceptable for current scale)
_PATHWAY_CACHE = None  # {'pathwayInfo': {...}, 'geneToPathwayMap': {...}}
_COMPLEX_CACHE = None  # {'complexMap': {...}}


def _cache_dir() -> str:
    """Resolve the static cache directory."""
    return os.path.join(os.path.dirname(__file__), '..', 'data', 'cache')


def _read_json_safe(path: str):
    """Read JSON file defensively. Returns None if not found."""
    try:
        # Use utf-8-sig to auto-handle BOM if present (robust to Windows-saved files)
        with open(path, 'r', encoding='utf-8-sig') as f:
            return json.load(f)
    except FileNotFoundError:
        return None
    except json.JSONDecodeError as e:
        # Malformed or empty JSON should not crash the server – log and continue fallback
        try:
            from flask import current_app as _ca
            _ca.logger.error(f"[complexes] JSON decode error at {path}: {e}")
        except Exception:
            pass
        return None
    except Exception as e:
        try:
            from flask import current_app as _ca
            _ca.logger.error(f"[complexes] Unexpected error reading {path}: {e}")
        except Exception:
            pass
        return None


def _canonicalize_complex_map(raw) -> Dict[str, List[str]]:
    """Ensure complexMap is dict[str, list[str]]; coerce/clean defensively."""
    out: Dict[str, List[str]] = {}
    if not isinstance(raw, dict):
        return out
    for k, v in raw.items():
        try:
            name = str(k)
            genes: List[str] = []
            if v is None:
                genes = []
            elif isinstance(v, list):
                genes = [str(g) for g in v if g is not None]
            elif isinstance(v, dict):
                # Accept shapes like {"genes": [...]} or similar
                cand = v.get('genes') or v.get('members') or []
                if isinstance(cand, list):
                    genes = [str(g) for g in cand if g is not None]
            else:
                # Single string
                genes = [str(v)]
            if genes:
                # unique + sorted
                out[name] = sorted(list({g for g in genes}))
            else:
                out[name] = []
        except Exception:
            continue
    return out


def _build_pathway_cache_from_local() -> Dict[str, Dict]:
    """Build pathwayInfo and geneToPathwayMap from local caches.

    Primary source: metabolic_pathways.json
    Fallback: pathway_data.json + reaction_gene_map.json
    """
    cache_dir = _cache_dir()
    pathway_info: Dict[str, Dict] = {}
    gene_to_pathways: Dict[str, Set[str]] = {}

    # Preferred comprehensive source
    metabolic = _read_json_safe(os.path.join(cache_dir, 'metabolic_pathways.json'))
    if metabolic:
        # Normalize pathway container: dict (id->obj), list, or object with 'data'
        pathways = []
        if isinstance(metabolic, dict):
            if 'data' in metabolic and isinstance(metabolic['data'], list):
                pathways = metabolic['data']
            else:
                pathways = list(metabolic.values())
        elif isinstance(metabolic, list):
            pathways = metabolic
        for p in pathways:
            pid = str(p.get('id') or p.get('pathway_id') or p.get('uid') or '')
            pname = p.get('name') or p.get('pathway_name') or ''
            if not pid or not pname:
                continue
            pathway_info[pid] = {'name': pname}

            # Direct genes (accept both 'gene_name' and 'name')
            for g in p.get('genes', []):
                if isinstance(g, dict):
                    gname = (g.get('gene_name') or g.get('name'))
                else:
                    gname = g if isinstance(g, str) else None
                if gname:
                    gene_to_pathways.setdefault(gname, set()).add(pid)
            # Enzymes → gene names
            for enz in p.get('enzymes', []):
                gname = (enz.get('gene_name') or (enz.get('gene', {}) or {}).get('name'))
                if gname:
                    gene_to_pathways.setdefault(gname, set()).add(pid)
            # Reactions → genes/enzymes (accept 'gene_name' and 'name')
            for rxn in p.get('reactions', []):
                for g in rxn.get('genes', []):
                    if isinstance(g, dict):
                        gname = (g.get('gene_name') or g.get('name'))
                    else:
                        gname = g if isinstance(g, str) else None
                    if gname:
                        gene_to_pathways.setdefault(gname, set()).add(pid)
                for enz in rxn.get('enzymes', []):
                    gname = (enz.get('gene_name') or (enz.get('gene', {}) or {}).get('name'))
                    if gname:
                        gene_to_pathways.setdefault(gname, set()).add(pid)

    # Fallback derivation
    if not pathway_info or not gene_to_pathways:
        pathway_data = _read_json_safe(os.path.join(cache_dir, 'pathway_data.json')) or {}
        reaction_gene_map = _read_json_safe(os.path.join(cache_dir, 'reaction_gene_map.json')) or {}
        for pid, meta in pathway_data.items():
            pname = meta.get('name')
            if pname:
                pathway_info[str(pid)] = {'name': pname}
            for rid in meta.get('reactions', []):
                for g in reaction_gene_map.get(str(rid), []):
                    if isinstance(g, dict):
                        gname = (g.get('gene_name') or g.get('name'))
                    else:
                        gname = g if isinstance(g, str) else None
                    if gname:
                        gene_to_pathways.setdefault(gname, set()).add(str(pid))

    if not pathway_info:
        current_app.logger.warning("Pathway cache built with 0 pathways. Check src/data/cache/metabolic_pathways.json schema.")
    if not gene_to_pathways:
        current_app.logger.warning("Pathway cache built with 0 gene links. Verify reaction/pathway files.")

    return {
        'pathwayInfo': pathway_info,
        'geneToPathwayMap': {g: sorted(list(pids)) for g, pids in gene_to_pathways.items()}
    }


def _build_complex_cache_from_local() -> Dict[str, Dict]:
    """
    ==========================================================================================
    == FINAL VERSION: Builds the complete complex and adjacency map cache from local files. ==
    ==========================================================================================
    This function is the single source of truth for loading complex data. It ensures that
    the adjacency map (for neighborhood search) is ALWAYS built from the rich complex data.
    """
    from flask import current_app
    from src.api.synonym_service import synonym_service
    
    current_app.logger.info("[Cache Builder] Starting _build_complex_cache_from_local...")
    
    # === Step 1: Load the primary complex map from its JSON file ===
    complex_map_path = os.path.join(_cache_dir(), 'protein_complex_map.json')
    raw_map_data = _read_json_safe(complex_map_path)
    
    if not raw_map_data:
        current_app.logger.error("[Cache Builder] CRITICAL: protein_complex_map.json is missing or empty. Cannot build cache.")
        return {'complexMap': {}, 'adjacency': {}}

    # The canonical map of complex names to gene lists
    final_complex_map = _canonicalize_complex_map(raw_map_data)
    current_app.logger.info(f"[Cache Builder] Successfully loaded {len(final_complex_map)} complexes.")

    # === Step 2: Build the Adjacency Map from the Complex Data ===
    # NOTE FOR AI: This is the primary logic for creating the interaction network.
    # It treats every complex as a "clique" where all members interact with each other.
    # This logic MUST ALWAYS run if a complex map is present.
    
    adjacency_from_complexes: Dict[str, Set[str]] = {}
    for complex_name, genes in final_complex_map.items():
        if len(genes) >= 2:
            # Create pairwise, undirected edges between all members of the complex
            for i in range(len(genes)):
                for j in range(i + 1, len(genes)):
                    gene_a, gene_b = genes[i], genes[j]
                    if gene_a == gene_b: continue
                    adjacency_from_complexes.setdefault(gene_a, set()).add(gene_b)
                    adjacency_from_complexes.setdefault(gene_b, set()).add(gene_a)
    
    current_app.logger.info(f"[Cache Builder] Built pre-canonical adjacency map with {len(adjacency_from_complexes)} nodes.")

    # === Step 3: Canonicalize Gene Names in the Adjacency Map ===
    # NOTE FOR AI: This step is crucial for consistency. It converts all gene names
    # to their canonical form using the synonym service. Do not remove this.
    
    canonical_adjacency: Dict[str, Set[str]] = {}
    try:
        for gene, neighbors in adjacency_from_complexes.items():
            canonical_gene = synonym_service.resolve_gene_name(gene) or gene
            
            # CRITICAL FIX: Ensure consistent lowercase gene names for API compatibility
            # The API expects lowercase gene names (e.g., 'ccpA'), but the data may contain
            # mixed case (e.g., 'CcpA'). Convert to lowercase for consistency.
            canonical_gene = canonical_gene.lower()
            
            # Ensure the key exists for the canonical gene name
            canonical_adjacency.setdefault(canonical_gene, set())

            for neighbor in neighbors:
                canonical_neighbor = synonym_service.resolve_gene_name(neighbor) or neighbor
                # Apply same lowercase conversion to neighbors
                canonical_neighbor = canonical_neighbor.lower()
                if canonical_gene != canonical_neighbor:
                    canonical_adjacency[canonical_gene].add(canonical_neighbor)

    except Exception as e:
        current_app.logger.error(f"[Cache Builder] Synonym canonicalization failed: {e}. Using raw names.")
        # If synonym service fails, use the un-canonicalized map with lowercase conversion
        canonical_adjacency = {}
        for gene, neighbors in adjacency_from_complexes.items():
            canonical_gene = gene.lower()
            canonical_adjacency.setdefault(canonical_gene, set())
            for neighbor in neighbors:
                canonical_neighbor = neighbor.lower()
                if canonical_gene != canonical_neighbor:
                    canonical_adjacency[canonical_gene].add(canonical_neighbor)

    # Final conversion from sets to sorted lists for deterministic output
    final_adjacency = {gene: sorted(list(neighbors)) for gene, neighbors in canonical_adjacency.items()}
    current_app.logger.info(f"[Cache Builder] Finished. Final adjacency map has {len(final_adjacency)} nodes.")
    
    if 'ccpa' in final_adjacency:
        current_app.logger.info(f"[Cache Builder] Sanity Check Passed: 'ccpa' is in the final adjacency map with {len(final_adjacency['ccpa'])} neighbors.")
    else:
        current_app.logger.warning("[Cache Builder] Sanity Check FAILED: 'ccpa' was NOT found in the final adjacency map.")

    return {'complexMap': final_complex_map, 'adjacency': final_adjacency}


def get_pathways_cache() -> Dict[str, Dict]:
    """Public accessor used by blueprint to get cached pathway data."""
    global _PATHWAY_CACHE
    if _PATHWAY_CACHE is None:
        _PATHWAY_CACHE = _build_pathway_cache_from_local()
    return _PATHWAY_CACHE


def get_complexes_cache() -> Dict[str, Dict]:
    """Public accessor used by blueprint to get cached complex data."""
    global _COMPLEX_CACHE
    if _COMPLEX_CACHE is None:
        _COMPLEX_CACHE = _build_complex_cache_from_local()
    return _COMPLEX_CACHE


def reload_complexes_cache() -> Dict[str, Dict]:
    """Force rebuild complex cache from local files and replace memoized value."""
    global _COMPLEX_CACHE
    _COMPLEX_CACHE = _build_complex_cache_from_local()
    return _COMPLEX_CACHE

def get_regulations_data():
    """Get all regulations data from SubtiWiki API"""
    if REGULATIONS_CACHE:
        return REGULATIONS_CACHE
    
    try:
        base_url = "https://www.subtiwiki.uni-goettingen.de/v5/api"
        response = requests.get(f"{base_url}/regulation/graph", timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            REGULATIONS_CACHE.update(data)
            return data
        else:
            print(f"Error fetching regulations: HTTP {response.status_code}")
            print(f"Response text: {response.text[:200]}...")
            return {}
            
    except requests.exceptions.Timeout:
        print("Error fetching regulations data: Request timed out after 30 seconds")
        return {}
    except requests.exceptions.ConnectionError:
        print("Error fetching regulations data: Connection error - check internet connectivity")
        return {}
    except requests.exceptions.RequestException as e:
        print(f"Error fetching regulations data: Request failed - {e}")
        return {}
    except Exception as e:
        print(f"Error fetching regulations data: {e}")
        return {}

def get_regulons_data():
    """Get all regulons data from SubtiWiki API"""
    if REGULONS_CACHE:
        return REGULONS_CACHE
    
    try:
        base_url = "https://www.subtiwiki.uni-goettingen.de/v5/api"
        print(f"Fetching regulons from {base_url}/regulon/")
        
        # Add timeout and better error handling
        response = requests.get(f"{base_url}/regulon/", timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            print(f"Regulons API response keys: {list(data.keys())}")
            print(f"Number of regulons: {len(data.get('data', []))}")
            REGULONS_CACHE.update(data)
            return data
        else:
            print(f"Error fetching regulons: HTTP {response.status_code}")
            print(f"Response text: {response.text[:200]}...")
            return {}
            
    except requests.exceptions.Timeout:
        print("Error fetching regulons data: Request timed out after 30 seconds")
        return {}
    except requests.exceptions.ConnectionError:
        print("Error fetching regulons data: Connection error - check internet connectivity")
        return {}
    except requests.exceptions.RequestException as e:
        print(f"Error fetching regulons data: Request failed - {e}")
        return {}
    except Exception as e:
        print(f"Error fetching regulons data: {e}")
        return {}

def get_regulations_list():
    """Get list of all regulations for filter panel"""
    regulations_data = get_regulations_data()
    regulations = regulations_data.get('regulations', [])
    
    # Group regulations by regulator
    regulator_groups = {}
    
    for reg in regulations:
        regulator_id = reg.get('regulator_gene_id')
        regulator_name = reg.get('regulator_name', f"Regulator_{regulator_id}")
        
        if regulator_id not in regulator_groups:
            regulator_groups[regulator_id] = {
                'id': regulator_id,
                'name': regulator_name,
                'regulated_genes': set(),
                'regulation_count': 0
            }
        
        regulated_gene_id = reg.get('regulated_gene_id')
        if regulated_gene_id:
            regulator_groups[regulator_id]['regulated_genes'].add(regulated_gene_id)
            regulator_groups[regulator_id]['regulation_count'] += 1
    
    # Convert to list format
    regulations_list = []
    for regulator_id, data in regulator_groups.items():
        regulations_list.append({
            'id': regulator_id,
            'name': data['name'],
            'regulated_genes': list(data['regulated_genes']),
            'regulation_count': data['regulation_count']
        })
    
    # Sort by regulation count
    regulations_list.sort(key=lambda x: x['regulation_count'], reverse=True)
    return regulations_list

def get_regulons_list():
    """Get list of all regulons for filter panel from local mapping file"""
    try:
        # Load regulon mapping from local file
        mapping_file_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'regulon_gene_mapping.json')
        
        if os.path.exists(mapping_file_path):
            with open(mapping_file_path, 'r', encoding='utf-8') as f:
                mapping_data = json.load(f)
            
            regulon_to_genes = mapping_data.get('regulon_to_genes', {})
            print(f"Loaded {len(regulon_to_genes)} regulons from local mapping file")
            
            regulons_list = []
            for regulon_id, regulon_data in regulon_to_genes.items():
                regulons_list.append({
                    'id': regulon_data['id'],
                    'name': regulon_data['name'],
                    'regulator_gene_id': regulon_data.get('regulator_gene_id'),
                    'regulator_gene_name': regulon_data.get('regulator_gene_name'),
                    'regulated_genes': regulon_data['genes'],
                    'regulation_count': regulon_data['gene_count'],
                    'gene_count': regulon_data['gene_count'],
                    'function': regulon_data.get('function', '')
                })
            
            # Sort alphabetically by name
            regulons_list.sort(key=lambda x: x['name'].lower())
            print(f"Successfully processed {len(regulons_list)} regulons from local mapping")
            return regulons_list
        else:
            print("Regulon mapping file not found, falling back to API")
            return get_regulons_list_from_api()
            
    except Exception as e:
        print(f"Error loading regulons from local mapping: {e}")
        print("Falling back to API")
        return get_regulons_list_from_api()

def get_regulons_list_from_api():
    """Get list of all regulons for filter panel from API (fallback)"""
    regulons_data = get_regulons_data()
    regulons = regulons_data.get('data', []) if isinstance(regulons_data, dict) else regulons_data
    
    print(f"Processing {len(regulons)} regulons from API")
    
    # If no regulons data available, return empty list - no fallback data
    if not regulons:
        print("No regulons data available from API - returning empty list")
        return []
    
    regulons_list = []
    
    for regulon in regulons:
        try:
            regulon_id = regulon.get('id')
            regulator_gene = regulon.get('regulator_gene', {})
            
            # Handle different regulator name fields
            regulator_name = regulon.get('regulator_display_name')
            if not regulator_name and regulator_gene:
                regulator_name = regulator_gene.get('name')
            if not regulator_name:
                regulator_name = f"Regulon_{regulon_id}"
            
            # Get regulated genes from gene_regulations
            regulated_genes = set()
            gene_regulations = regulon.get('gene_regulations', [])
            
            for gene_reg in gene_regulations:
                if isinstance(gene_reg, dict):
                    gene = gene_reg.get('gene', {})
                    if isinstance(gene, dict):
                        gene_id = gene.get('id')
                        if gene_id:
                            regulated_genes.add(gene_id)
            
            # Get regulated genes from operon_regulations
            operon_regulations = regulon.get('operon_regulations', [])
            for operon_reg in operon_regulations:
                if isinstance(operon_reg, dict):
                    operon = operon_reg.get('operon', {})
                    if isinstance(operon, dict):
                        operon_genes = operon.get('genes', [])
                        for gene in operon_genes:
                            if isinstance(gene, dict):
                                gene_id = gene.get('id')
                                if gene_id:
                                    regulated_genes.add(gene_id)
            
            # Extract only the function field from the description
            function_description = ""
            description = regulon.get('description', '')
            if description:
                # Look for "Function:" in the description
                if 'Function:' in description:
                    # Extract the function part
                    function_start = description.find('Function:')
                    function_end = description.find(';', function_start)
                    if function_end == -1:
                        function_end = len(description)
                    function_description = description[function_start:function_end].replace('Function:', '').strip()
                else:
                    # If no "Function:" found, use the first sentence
                    sentences = description.split(';')
                    if sentences:
                        function_description = sentences[0].strip()
            
            regulons_list.append({
                'id': regulon_id,
                'name': regulator_name,
                'regulator_gene_id': regulator_gene.get('id') if regulator_gene else None,
                'regulator_gene_name': regulator_gene.get('name') if regulator_gene else None,
                'regulated_genes': list(regulated_genes),
                'regulation_count': len(regulated_genes),
                'gene_count': len(regulated_genes),  # ADDED: Frontend expects gene_count
                'function': function_description
            })
        except Exception as e:
            print(f"Error processing regulon {regulon.get('id', 'unknown')}: {e}")
            continue
    
    # Sort alphabetically by name
    regulons_list.sort(key=lambda x: x['name'].lower())
    print(f"Successfully processed {len(regulons_list)} regulons")
    return regulons_list

def get_genes_for_regulation(regulation_id: int) -> List[str]:
    """Get genes regulated by a specific regulation"""
    regulations_data = get_regulations_data()
    regulations = regulations_data.get('regulations', [])
    
    regulated_genes = []
    for reg in regulations:
        if reg.get('regulator_gene_id') == regulation_id:
            regulated_gene_id = reg.get('regulated_gene_id')
            if regulated_gene_id:
                # Convert gene ID to gene name
                gene_name = get_gene_name_by_id(regulated_gene_id)
                if gene_name:
                    regulated_genes.append(gene_name)
    
    return regulated_genes

def get_genes_for_regulon(regulon_id: int) -> List[str]:
    """Get genes regulated by a specific regulon from local mapping"""
    try:
        # Load regulon mapping from local file
        mapping_file_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'regulon_gene_mapping.json')
        
        if os.path.exists(mapping_file_path):
            with open(mapping_file_path, 'r', encoding='utf-8') as f:
                mapping_data = json.load(f)
            
            regulon_to_genes = mapping_data.get('regulon_to_genes', {})
            regulon_data = regulon_to_genes.get(str(regulon_id))
            
            if regulon_data:
                return regulon_data['genes']
            else:
                print(f"Regulon {regulon_id} not found in local mapping")
                return []
        else:
            print("Regulon mapping file not found, falling back to API")
            return get_genes_for_regulon_from_api(regulon_id)
            
    except Exception as e:
        print(f"Error loading regulon {regulon_id} from local mapping: {e}")
        print("Falling back to API")
        return get_genes_for_regulon_from_api(regulon_id)

def get_genes_for_regulon_from_api(regulon_id: int) -> List[str]:
    """Get genes regulated by a specific regulon from API (fallback)"""
    regulons_data = get_regulons_data()
    regulons = regulons_data.get('data', []) if isinstance(regulons_data, dict) else regulons_data
    
    for regulon in regulons:
        if regulon.get('id') == regulon_id:
            regulated_genes = set()  # Use set to prevent duplicates
            
            # Get genes from gene_regulations
            gene_regulations = regulon.get('gene_regulations', [])
            for gene_reg in gene_regulations:
                if isinstance(gene_reg, dict):
                    gene = gene_reg.get('gene', {})
                    if isinstance(gene, dict):
                        gene_name = gene.get('name')
                        if gene_name:
                            regulated_genes.add(gene_name)  # Use add() for set
            
            # Get genes from operon_regulations
            operon_regulations = regulon.get('operon_regulations', [])
            for operon_reg in operon_regulations:
                if isinstance(operon_reg, dict):
                    operon = operon_reg.get('operon', {})
                    if isinstance(operon, dict):
                        operon_genes = operon.get('genes', [])
                        for gene in operon_genes:
                            if isinstance(gene, dict):
                                gene_name = gene.get('name')
                                if gene_name:
                                    regulated_genes.add(gene_name)  # Use add() for set
            
            # Convert back to sorted list for consistent ordering
            return sorted(list(regulated_genes))
    
    return []

def get_gene_name_by_id(gene_id: int) -> Optional[str]:
    """Convert gene ID to gene name"""
    # This is a simplified lookup - in practice you might want to cache this mapping
    for gene_name, id_val in GENE_NAME_TO_ID.items():
        if id_val == gene_id:
            return gene_name
    return None

def get_genes_for_multi_selection(selection_data: Dict) -> Dict:
    """
    Get genes for multi-panel selection with intersection highlighting
    
    Args:
        selection_data: {
            'groups': [1, 3, 200],  // group IDs
            'genes': ['copA', 'dapF'],  // individual gene names
            'regulations': [1, 2],  // regulation IDs
            'regulons': [1, 3]  // regulon IDs
        }
    
    Returns:
        Dict with genes organized by category and intersections
    """
    all_genes = {}
    gene_sources = {}  # Track which categories each gene belongs to
    
    # Get plot data to check available genes
    from flask import session
    plots = []
    processed_file_id = session.get('processed_data_file', None)
    current_app.logger.info(f"DEBUG: processed_file_id from session: {processed_file_id}")
    
    if processed_file_id:
        processed_folder = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')
        file_path = os.path.join(processed_folder, processed_file_id)
        current_app.logger.info(f"DEBUG: Looking for plot data file: {file_path}")
        current_app.logger.info(f"DEBUG: File exists: {os.path.exists(file_path)}")
        
        try:
            with open(file_path, 'r') as f:
                plots = json.load(f)
            current_app.logger.info(f"DEBUG: Successfully loaded {len(plots)} plot traces")
        except Exception as e:
            current_app.logger.warning(f"Error reading processed data file: {e}")
    else:
        current_app.logger.warning("DEBUG: No processed_file_id found in session!")
    
    # Get genes from groups
    category_names = {}  # Store display names for categories
    if 'groups' in selection_data:
        for group_id in selection_data['groups']:
            group_genes = get_genes_in_category(group_id)
            if group_genes:
                # Get group display name
                group_cat = next((cat for cat in DETAILED_CATEGORIES if cat['id'] == group_id), None)
                group_name = group_cat['name'] if group_cat else f'Group {group_id}'
                category_key = f'group_{group_id}'
                category_names[category_key] = group_name
                
                all_genes[category_key] = group_genes
                for gene in group_genes:
                    if gene not in gene_sources:
                        gene_sources[gene] = set()
                    gene_sources[gene].add(category_key)
    
    # Get genes from individual selection
    if 'genes' in selection_data and selection_data['genes']:
        individual_genes = selection_data['genes']
        all_genes['individual'] = individual_genes
        for gene in individual_genes:
            if gene not in gene_sources:
                gene_sources[gene] = set()
            gene_sources[gene].add('individual')
    
    # Get genes from regulations
    if 'regulations' in selection_data:
        for regulation_id in selection_data['regulations']:
            regulation_genes = get_genes_for_regulation(regulation_id)
            if regulation_genes:
                all_genes[f'regulation_{regulation_id}'] = regulation_genes
                for gene in regulation_genes:
                    if gene not in gene_sources:
                        gene_sources[gene] = set()
                    gene_sources[gene].add(f'regulation_{regulation_id}')
    
    # Get genes from regulons
    if 'regulons' in selection_data and selection_data['regulons']:
        print(f"[DEBUG] Processing regulons: {selection_data['regulons']}")
        regulons_data = get_regulons_data()
        regulons = regulons_data.get('data', []) if isinstance(regulons_data, dict) else regulons_data
        for regulon_id in selection_data['regulons']:
            print(f"[DEBUG] Processing regulon ID: {regulon_id}")
            regulon_genes = get_genes_for_regulon(int(regulon_id))
            print(f"[DEBUG] Regulon {regulon_id} genes: {len(regulon_genes) if regulon_genes else 0}")
            if regulon_genes:
                # Get regulon display name
                regulon = next((r for r in regulons if r.get('id') == int(regulon_id)), None)
                regulon_name = regulon.get('regulator_display_name') or regulon.get('regulator_gene', {}).get('name') if regulon else f'Regulon {regulon_id}'
                category_key = f'regulon_{regulon_id}'
                category_names[category_key] = f'{regulon_name} regulon'
                
                all_genes[category_key] = regulon_genes
                for gene in regulon_genes:
                    if gene not in gene_sources:
                        gene_sources[gene] = set()
                    gene_sources[gene].add(category_key)
                print(f"[DEBUG] Added regulon category: {category_key} -> {regulon_name} regulon")
            else:
                print(f"[DEBUG] No genes found for regulon {regulon_id}")

    # Get genes from pathways (IDs)
    if 'pathways' in selection_data and selection_data['pathways']:
        try:
            path_map = get_genes_for_pathway_ids(selection_data['pathways'])
            pcache = get_pathways_cache()
            for pid, genes in path_map.items():
                category_key = f'pathway_{pid}'
                display = pcache['pathwayInfo'].get(str(pid), {}).get('name', f'Pathway {pid}')
                category_names[category_key] = display
                all_genes[category_key] = genes
                for gene in genes:
                    gene_sources.setdefault(gene, set()).add(category_key)
        except Exception as e:
            current_app.logger.warning(f"Pathways selection processing failed: {e}")

    # Get genes from complexes (names)
    if 'complexes' in selection_data and selection_data['complexes']:
        try:
            cx_map = get_genes_for_complex_names(selection_data['complexes'])
            for name, genes in cx_map.items():
                category_key = f'complex_{name}'
                category_names[category_key] = name
                all_genes[category_key] = genes
                for gene in genes:
                    gene_sources.setdefault(gene, set()).add(category_key)
        except Exception as e:
            current_app.logger.warning(f"Complex selection processing failed: {e}")
    
    # Get genes from expression clusters
    current_app.logger.info(f"DEBUG: Checking expression clusters in selection_data: {selection_data.get('expressionClusters', 'NOT_FOUND')}")
    if 'expressionClusters' in selection_data and selection_data['expressionClusters']:
        current_app.logger.info(f"DEBUG: Expression clusters found: {selection_data['expressionClusters']}")
        try:
            # Load the gene expression clusters JSON file
            clusters_file_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'gene_expression_clusters.json')
            current_app.logger.info(f"DEBUG: Looking for clusters file: {clusters_file_path}")
            current_app.logger.info(f"DEBUG: Clusters file exists: {os.path.exists(clusters_file_path)}")
            
            if os.path.exists(clusters_file_path):
                with open(clusters_file_path, 'r', encoding='utf-8') as f:
                    gene_clusters = json.load(f)
                
                # Load gene name mapping to resolve gene name mismatches
                mapping_file_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'gene_name_mapping.json')
                gene_mapping = {}
                if os.path.exists(mapping_file_path):
                    with open(mapping_file_path, 'r', encoding='utf-8') as f:
                        mapping_data = json.load(f)
                        gene_mapping = mapping_data.get('cluster_to_plot', {})
                        current_app.logger.info(f"Loaded gene name mapping with {len(gene_mapping)} mappings")
                
                # Get all genes present in the volcano plot data for validation
                plot_genes = set()
                for plot_data in plots:
                    if 'text' in plot_data:
                        plot_genes.update(plot_data['text'])
                
                # Also check gene_names field if available
                for plot_data in plots:
                    if 'gene_names' in plot_data:
                        plot_genes.update(plot_data['gene_names'])
                
                current_app.logger.info(f"DEBUG: Total unique genes in plot: {len(plot_genes)}")
                current_app.logger.info(f"DEBUG: Sample plot genes: {list(plot_genes)[:5]}")
                
                # Process each selected cluster
                for cluster_id in selection_data['expressionClusters']:
                    current_app.logger.info(f"DEBUG: Processing cluster {cluster_id}")
                    cluster_genes = []
                    for gene_name, clusters in gene_clusters.items():
                        if cluster_id in clusters:
                            cluster_genes.append(gene_name)
                    
                    current_app.logger.info(f"DEBUG: Cluster {cluster_id} contains {len(cluster_genes)} genes")
                    current_app.logger.info(f"DEBUG: Sample cluster genes: {cluster_genes[:5]}")
                    
                    if cluster_genes:
                        # Use gene mapping to convert cluster gene names to plot gene names
                        mapped_genes = []
                        for cluster_gene in cluster_genes:
                            # Check if we have a mapping for this gene
                            if cluster_gene in gene_mapping:
                                mapped_gene = gene_mapping[cluster_gene]
                                if mapped_gene in plot_genes:
                                    mapped_genes.append(mapped_gene)
                                else:
                                    current_app.logger.debug(f"DEBUG: Mapped gene '{mapped_gene}' not found in plot data")
                            else:
                                # Try direct match if no mapping available
                                if cluster_gene in plot_genes:
                                    mapped_genes.append(cluster_gene)
                                else:
                                    current_app.logger.debug(f"DEBUG: No mapping found for cluster gene '{cluster_gene}'")
                        
                        current_app.logger.info(f"DEBUG: Cluster {cluster_id} has {len(mapped_genes)} genes that exist in plot after mapping")
                        
                        if mapped_genes:
                            category_key = f'expression_cluster_{cluster_id}'
                            category_names[category_key] = f'Expression Cluster {cluster_id}'
                            all_genes[category_key] = mapped_genes
                            for gene in mapped_genes:
                                gene_sources.setdefault(gene, set()).add(category_key)
                            
                            # Log mapping statistics
                            mapped_count = len([g for g in cluster_genes if g in gene_mapping])
                            direct_count = len([g for g in cluster_genes if g in plot_genes])
                            current_app.logger.info(f"Cluster {cluster_id}: {len(mapped_genes)} genes found in plot ({mapped_count} via mapping, {direct_count} direct)")
                        else:
                            current_app.logger.warning(f"Cluster {cluster_id}: No genes from this cluster are present in the volcano plot data after mapping. Cluster contains: {cluster_genes[:5]}{'...' if len(cluster_genes) > 5 else ''}")
            else:
                current_app.logger.warning("Expression clusters data file not found")
        except Exception as e:
            current_app.logger.warning(f"Expression clusters selection processing failed: {e}")
    
    # Find intersections (genes that appear in ALL selected categories)
    intersection_genes = []
    total_categories = len(all_genes)
    
    for gene, sources in gene_sources.items():
        # Only include genes that appear in ALL selected categories
        if len(sources) == total_categories and total_categories > 1:
            intersection_genes.append({
                'gene': gene,
                'sources': list(sources),
                'source_count': len(sources)
            })
    
    # Sort intersections by number of sources (most intersections first)
    intersection_genes.sort(key=lambda x: x['source_count'], reverse=True)
    
    # Convert sets to lists for JSON serialization
    gene_sources_serializable = {}
    for gene, sources in gene_sources.items():
        gene_sources_serializable[gene] = list(sources)
    
    return {
        'categories': all_genes,
        'category_names': category_names,  # ADDED: Display names for legend
        'intersections': intersection_genes,
        'gene_sources': gene_sources_serializable,
        'total_unique_genes': len(gene_sources)
    }


# === ADD: Listing helpers for pathways/complexes ===
def _invert_gene_to_pathway(g2p: Dict[str, List[str]], selected_ids) -> Dict[str, List[str]]:
    pid_to_genes: Dict[str, Set[str]] = {str(pid): set() for pid in selected_ids}
    for gene, pids in g2p.items():
        for pid in pids:
            s = str(pid)
            if s in pid_to_genes:
                pid_to_genes[s].add(gene)
    return {pid: sorted(list(genes)) for pid, genes in pid_to_genes.items()}


def get_genes_for_pathway_ids(pathway_ids: List[str]) -> Dict[str, List[str]]:
    pc = get_pathways_cache()
    # Debug: log available keys and requested IDs
    req = [str(pid) for pid in pathway_ids]
    have = list(pc['pathwayInfo'].keys())
    current_app.logger.info(f"[pathways] requested={req[:10]}... total_req={len(req)} available_keys={len(have)}")
    result = _invert_gene_to_pathway(pc['geneToPathwayMap'], req)
    # Per-ID debug
    for pid in req[:10]:
        genes = result.get(pid, [])
        current_app.logger.info(f"[pathways] pid={pid} -> genes={len(genes)}")
    return result


def get_genes_for_complex_names(complex_names: List[str]) -> Dict[str, List[str]]:
    cache = get_complexes_cache()
    cm = cache.get('complexMap', {})
    # Debug: log counts
    current_app.logger.info(f"[complexes] available_complexes={len(cm)} requested={len(complex_names)} sample_req={complex_names[:5]}")
    out = {name: sorted(list(set(cm.get(name, [])))) for name in complex_names}
    for n in complex_names[:5]:
        current_app.logger.info(f"[complexes] name={n} -> genes={len(out.get(n, []))}")
    return out


def list_pathways() -> List[Dict]:
    pc = get_pathways_cache()
    info = pc['pathwayInfo']
    pid_to_genes = _invert_gene_to_pathway(pc['geneToPathwayMap'], info.keys())
    items = [
        {'id': pid, 'name': meta.get('name', f'Pathway {pid}'), 'gene_count': len(pid_to_genes.get(pid, []))}
        for pid, meta in info.items()
    ]
    # Filter out zero-gene pathways to avoid confusing UX
    return [it for it in items if it['gene_count'] > 0]


def list_complexes() -> List[Dict]:
    cm = get_complexes_cache()['complexMap']
    # Include full genes list so frontend can show accurate counts without extra lookups
    return [
        {'id': name, 'name': name, 'genes': genes, 'gene_count': len(genes)}
        for name, genes in cm.items()
    ]


# =====================
# Metabolites (cache-backed)
# =====================
try:
    from src.api.metabolite_cache import (
        get_gene_to_metabolites,
        get_metabolite_to_genes,
        get_metabolite_categories,
        get_metabolite_to_categories,
    )
except Exception:
    def get_gene_to_metabolites():
        return {}
    def get_metabolite_to_genes():
        return {}
    def get_metabolite_categories():
        return []
    def get_metabolite_to_categories():
        return {}


def list_metabolite_categories() -> List[Dict]:
    """Return cached metabolite categories (flat list with children ids)."""
    cats = get_metabolite_categories() or []
    # Normalize keys and ensure determinism
    out = []
    for c in cats:
        out.append({
            'id': c.get('id'),
            'name': c.get('name'),
            'dot_notation': c.get('dot_notation'),
            'parent_id': c.get('parent_id'),
            'children': c.get('children') or [],
        })
    return out


def _metabolite_id_to_name() -> Dict[str, str]:
    """Build metabolite id → name map from gene_to_metabolites cache."""
    id_to_name: Dict[str, str] = {}
    g2m = get_gene_to_metabolites() or {}
    for entries in g2m.values():
        if not isinstance(entries, list):
            continue
        for e in entries:
            mid = e.get('metabolite_id')
            nm = e.get('metabolite_name')
            if mid is not None and nm:
                id_to_name[str(mid)] = nm
    return id_to_name


def list_metabolites_by_category(category_id: int) -> List[Dict]:
    """Return metabolites (id, name) belonging to a category id."""
    m2c = get_metabolite_to_categories() or {}
    id_to_name = _metabolite_id_to_name()
    items: List[Dict] = []
    for mid_key, cats in (m2c.items() if isinstance(m2c, dict) else []):
        try:
            mid_str = str(mid_key)
            for c in cats or []:
                if int(c.get('id')) == int(category_id):
                    items.append({'id': int(mid_str), 'name': id_to_name.get(mid_str, f'Metabolite {mid_str}')})
                    break
        except Exception:
            continue
    # Deduplicate & sort
    seen = set()
    uniq = []
    for it in items:
        if it['id'] not in seen:
            uniq.append(it)
            seen.add(it['id'])
    return sorted(uniq, key=lambda x: x['name'].lower())

def get_main_groups():
    """
    Return main parent groups for initial flat list selection.
    Focuses on Level 1 and Level 2 categories with significant gene counts.
    """
    main_groups = []
    
    # Get level 1 and 2 categories that have genes (via children)
    for cat in DETAILED_CATEGORIES:
        dot = cat.get('dot_notation', '')
        if not dot:
            continue
            
        level = len(dot.split('.'))
        
        # Level 1 (top categories) and Level 2 (major subcategories)
        if level <= 2:
            # Get total gene count for this category (including children)
            total_genes = get_genes_in_category(cat['id'])
            gene_count = len(total_genes) if total_genes else 0
            
            # Only include categories with substantial gene counts
            if gene_count >= 10:
                main_groups.append({
                    'id': cat['id'],
                    'name': cat['name'],
                    'gene_count': gene_count,
                    'dot_notation': dot,
                    'level': level
                })
    
    # Sort by gene count (largest first) and return top groups
    main_groups.sort(key=lambda x: x['gene_count'], reverse=True)
    return main_groups[:20]  # Top 20 most gene-rich categories

def get_hierarchical_groups():
    """
    Return hierarchical tree structure of all groups for tree view.
    Organizes categories by their dot notation hierarchy.
    """
    # Get all categories with dot notation
    categories_with_dot = [cat for cat in DETAILED_CATEGORIES if cat.get('dot_notation')]
    
    # Sort by dot notation for proper hierarchy
    categories_with_dot.sort(key=lambda x: [int(n) for n in x['dot_notation'].split('.')])
    
    # Build hierarchical structure
    tree = {}
    
    for cat in categories_with_dot:
        dot_parts = cat['dot_notation'].split('.')
        level = len(dot_parts)
        
        # Get gene count for this category
        total_genes = get_genes_in_category(cat['id'])
        gene_count = len(total_genes) if total_genes else 0
        
        # Create category node
        cat_node = {
            'id': cat['id'],
            'name': cat['name'],
            'dot_notation': cat['dot_notation'],
            'level': level,
            'gene_count': gene_count,
            'children': {}
        }
        
        # Place in tree structure
        current_level = tree
        for i, part in enumerate(dot_parts[:-1]):
            # Find the parent node
            parent_key = '.'.join(dot_parts[:i+1])
            if parent_key not in current_level:
                # Create parent node if it doesn't exist
                parent_cat = next((c for c in categories_with_dot if c['dot_notation'] == parent_key), None)
                if parent_cat:
                    parent_genes = get_genes_in_category(parent_cat['id'])
                    parent_gene_count = len(parent_genes) if parent_genes else 0
                    current_level[parent_key] = {
                        'id': parent_cat['id'],
                        'name': parent_cat['name'],
                        'dot_notation': parent_cat['dot_notation'],
                        'level': len(parent_key.split('.')),
                        'gene_count': parent_gene_count,
                        'children': {}
                    }
                else:
                    # This shouldn't happen if data is well-formed
                    break
            current_level = current_level[parent_key]['children']
        
        current_level[cat['dot_notation']] = cat_node
    
    return tree

def get_category_children(category_id):
    """
    Get direct children of a category for dynamic tree loading.
    """
    try:
        category_id_int = int(category_id)
    except (TypeError, ValueError):
        return []
    
    # Find the category
    parent_cat = next((cat for cat in DETAILED_CATEGORIES if cat['id'] == category_id_int), None)
    if not parent_cat:
        return []
    
    children = []
    children_ids = parent_cat.get('children_ids', [])
    
    for child_id in children_ids:
        child_cat = next((cat for cat in DETAILED_CATEGORIES if cat['id'] == child_id), None)
        if child_cat:
            # Get gene count
            total_genes = get_genes_in_category(child_id)
            gene_count = len(total_genes) if total_genes else 0
            
            children.append({
                'id': child_cat['id'],
                'name': child_cat['name'],
                'dot_notation': child_cat.get('dot_notation', ''),
                'level': len(child_cat.get('dot_notation', '').split('.')) if child_cat.get('dot_notation') else 0,
                'gene_count': gene_count,
                'has_children': len(child_cat.get('children_ids', [])) > 0
            })
    
    # Sort by dot notation
    children.sort(key=lambda x: x['dot_notation'])
    return children

def get_genes_for_selection(selection_data):
    """
    Get genes based on selection data (groups and/or individual genes).
    
    Args:
        selection_data (dict): {
            'groups': [group_id1, group_id2, ...],
            'genes': [gene_name1, gene_name2, ...]
        }
    
    Returns:
        dict: {
            'genes': [list of gene names],
            'sources': {gene_name: [source1, source2, ...]},
            'total_count': int
        }
    """
    all_selected_genes = set()
    gene_sources = {}  # Track where each gene came from
    
    # Process group selections
    for group_id in selection_data.get('groups', []):
        group_genes = get_genes_in_category(group_id)
        if group_genes:
            group_name = next((cat['name'] for cat in DETAILED_CATEGORIES if cat['id'] == int(group_id)), f"Group {group_id}")
            
            for gene in group_genes:
                all_selected_genes.add(gene)
                if gene not in gene_sources:
                    gene_sources[gene] = []
                gene_sources[gene].append(group_name)
    
    # Process individual gene selections
    for gene_name in selection_data.get('genes', []):
        if gene_name in GENE_NAME_TO_ID:
            all_selected_genes.add(gene_name)
            if gene_name not in gene_sources:
                gene_sources[gene_name] = []
            gene_sources[gene_name].append('Manual selection')
    
    return {
        'genes': list(all_selected_genes),
        'sources': gene_sources,
        'total_count': len(all_selected_genes)
    }


def get_genes_for_metabolite_selection(payload: Dict) -> Dict[str, List[str]]:
    """Compute genes for selected metabolite ids and/or category ids.

    payload keys:
      - metabolite_ids: List[int]
      - metabolite_category_ids: List[int]
      - interaction_types: Optional[List[str]] (Binding, Cofactor, Effector, Transport, Transport:Import, Transport:Export)
    """
    metabolite_ids = set(payload.get('metabolite_ids') or [])
    category_ids = set(payload.get('metabolite_category_ids') or [])
    interaction_types = set(payload.get('interaction_types') or [])

    result: Dict[str, List[str]] = {}

    # Expand category ids to metabolite ids
    if category_ids:
        expanded: List[int] = []
        for cid in category_ids:
            for m in list_metabolites_by_category(int(cid)):
                expanded.append(int(m['id']))
        for mid in expanded:
            metabolite_ids.add(int(mid))

    m2g = get_metabolite_to_genes() or {}
    g2m = get_gene_to_metabolites() or {}

    def gene_has_interaction_with(mid: int, gene: str) -> bool:
        if not interaction_types:
            return True
        for e in g2m.get(gene, []) or []:
            if e.get('metabolite_id') == mid:
                etype = e.get('type')
                if etype == 'Transport' and (f"Transport:{e.get('transport_type')}" in interaction_types):
                    return True
                if etype in interaction_types:
                    return True
        return False

    # Buckets per metabolite id
    for mid in metabolite_ids:
        key = f"metabolite_{int(mid)}"
        genes = []
        raw = m2g.get(str(int(mid))) or m2g.get(int(mid)) or []
        for g in raw:
            if gene_has_interaction_with(int(mid), g):
                genes.append(g)
        result[key] = sorted(list(set(genes)))

    # Buckets per category id (aggregated)
    for cid in category_ids:
        key = f"metabolite_category_{int(cid)}"
        agg: Set[str] = set()
        for m in list_metabolites_by_category(int(cid)):
            raw = m2g.get(str(int(m['id']))) or m2g.get(int(m['id'])) or []
            for g in raw:
                if gene_has_interaction_with(int(m['id']), g):
                    agg.add(g)
        result[key] = sorted(list(agg))

    return result

def search_genes(query, limit=50):
    """
    Search for genes matching the query for autocomplete.
    
    Args:
        query (str): Search term
        limit (int): Maximum results to return
    
    Returns:
        list: Matching gene names
    """
    if not query or len(query) < 2:
        return []
    
    query_lower = query.lower()
    matches = []
    
    # Search through all genes
    for gene_name in GENE_NAME_TO_ID.keys():
        if query_lower in gene_name.lower():
            matches.append(gene_name)
            if len(matches) >= limit:
                break
    
    # Sort by relevance (exact start matches first)
    matches.sort(key=lambda x: (
        not x.lower().startswith(query_lower),  # Starts with query first
        len(x),  # Shorter names first
        x.lower()  # Alphabetical
    ))
    
    return matches

def get_gene_info_detailed(gene_name):
    """
    Get detailed information about a gene including all groups it belongs to.
    
    Args:
        gene_name (str): Gene name
    
    Returns:
        dict: Gene information with groups
    """
    if gene_name not in GENE_NAME_TO_ID:
        return None
    
    # Find all groups this gene belongs to
    gene_groups = []
    for cat_name, genes in GENE_CATEGORY_MAP.items():
        if gene_name in genes:
            # Get category details
            cat_details = next((cat for cat in DETAILED_CATEGORIES if cat['name'] == cat_name), None)
            if cat_details:
                gene_groups.append({
                    'id': cat_details['id'],
                    'name': cat_name,
                    'dot_notation': cat_details.get('dot_notation', ''),
                    'level': len(cat_details.get('dot_notation', '').split('.')) if cat_details.get('dot_notation') else 0
                })
    
    # Sort groups by hierarchy level (top-level first)
    gene_groups.sort(key=lambda x: (x['level'], x['name']))
    
    return {
        'gene_name': gene_name,
        'gene_id': GENE_NAME_TO_ID[gene_name],
        'groups': gene_groups,
        'group_count': len(gene_groups)
    }
