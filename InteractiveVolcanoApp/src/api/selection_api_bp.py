"""
Flask Blueprint for Selection API endpoints
"""

from flask import Blueprint, jsonify, request
from src.api.selection_api import (
    get_main_groups,
    get_hierarchical_groups,
    get_category_children,
    get_genes_for_selection,
    get_genes_for_multi_selection,
    search_genes,
    get_gene_info_detailed,
    get_regulations_list,
    get_regulons_list,
    get_genes_for_regulation,
    get_genes_for_regulon,
    get_pathways_cache,
    get_complexes_cache,
    list_pathways,
    list_complexes,
    reload_complexes_cache
)
from src.api.selection_api import (
    list_metabolite_categories,
    list_metabolites_by_category,
    get_genes_for_metabolite_selection,
)
from src.api.synonym_service import synonym_service
import os
import json

selection_api_bp = Blueprint('selection_api_bp', __name__)

@selection_api_bp.route('/groups/main')
def main_groups():
    """Get main groups for selection UI (flat list)"""
    try:
        groups = get_main_groups()
        return jsonify({
            'success': True,
            'groups': groups,
            'count': len(groups)
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500


@selection_api_bp.route('/expression-clusters/hierarchical')
def expression_clusters_hierarchical():
    """Get expression clusters organized hierarchically (A, B, C as parent groups)"""
    try:
        # Load the gene expression clusters JSON file
        clusters_file_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'gene_expression_clusters.json')
        
        if not os.path.exists(clusters_file_path):
            return jsonify({
                'success': False,
                'error': 'Expression clusters data file not found'
            }), 404
        
        with open(clusters_file_path, 'r', encoding='utf-8') as f:
            gene_clusters = json.load(f)
        
        # Organize clusters hierarchically
        hierarchical_clusters = {
            'A': {},
            'B': {},
            'C': {}
        }
        
        # Process each gene and organize by cluster type
        for gene_name, clusters in gene_clusters.items():
            for cluster in clusters:
                if cluster.startswith('A'):
                    if cluster not in hierarchical_clusters['A']:
                        hierarchical_clusters['A'][cluster] = []
                    hierarchical_clusters['A'][cluster].append(gene_name)
                elif cluster.startswith('B'):
                    if cluster not in hierarchical_clusters['B']:
                        hierarchical_clusters['B'][cluster] = []
                    hierarchical_clusters['B'][cluster].append(gene_name)
                elif cluster.startswith('C'):
                    if cluster not in hierarchical_clusters['C']:
                        hierarchical_clusters['C'][cluster] = []
                    hierarchical_clusters['C'][cluster].append(gene_name)
        
        # Convert to the format expected by the frontend
        result = {
            'A': [{'id': cluster, 'name': cluster, 'gene_count': len(genes), 'genes': genes} 
                  for cluster, genes in hierarchical_clusters['A'].items()],
            'B': [{'id': cluster, 'name': cluster, 'gene_count': len(genes), 'genes': genes} 
                  for cluster, genes in hierarchical_clusters['B'].items()],
            'C': [{'id': cluster, 'name': cluster, 'gene_count': len(genes), 'genes': genes} 
                  for cluster, genes in hierarchical_clusters['C'].items()]
        }
        
        return jsonify({
            'success': True,
            'clusters': result,
            'total_genes': len(gene_clusters)
        })
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500


@selection_api_bp.route('/expression-clusters/flat')
def expression_clusters_flat():
    """Get expression clusters as a flat list for search functionality"""
    try:
        # Load the gene expression clusters JSON file
        clusters_file_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'gene_expression_clusters.json')
        
        if not os.path.exists(clusters_file_path):
            return jsonify({
                'success': False,
                'error': 'Expression clusters data file not found'
            }), 404
        
        with open(clusters_file_path, 'r', encoding='utf-8') as f:
            gene_clusters = json.load(f)
        
        # Create a flat list of all unique clusters
        all_clusters = set()
        for clusters in gene_clusters.values():
            all_clusters.update(clusters)
        
        # Count genes per cluster
        cluster_counts = {}
        for gene_name, clusters in gene_clusters.items():
            for cluster in clusters:
                if cluster not in cluster_counts:
                    cluster_counts[cluster] = []
                cluster_counts[cluster].append(gene_name)
        
        # Convert to list format
        flat_clusters = [
            {
                'id': cluster,
                'name': cluster,
                'gene_count': len(genes),
                'genes': genes
            }
            for cluster, genes in cluster_counts.items()
        ]
        
        # Sort by cluster name
        flat_clusters.sort(key=lambda x: x['name'])
        
        return jsonify({
            'success': True,
            'clusters': flat_clusters,
            'total_genes': len(gene_clusters)
        })
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500


@selection_api_bp.route('/selection/pathways/all', methods=['GET'])
def pathways_all():
    """Expose pathwayInfo and geneToPathwayMap from local cache."""
    try:
        built = get_pathways_cache()
        if not built.get('pathwayInfo') or not built.get('geneToPathwayMap'):
            return jsonify({
                'success': False,
                'error': {
                    'type': 'MissingData',
                    'message': 'Pathway cache empty (0 pathways or 0 linked genes).',
                    'citation': 'src/data/cache/metabolic_pathways.json',
                    'recommended_next_step': 'Rebuild caches; verify JSON schema/keys.'
                }
            }), 404
        return jsonify({'success': True, **built})
    except Exception as e:
        return jsonify({
            'success': False,
            'error': {
                'type': 'ServerError',
                'message': str(e),
                'citation': 'selection_api_bp.py:pathways_all',
                'recommended_next_step': 'Inspect server logs and cache JSON.'
            }
        }), 500


@selection_api_bp.route('/selection/complexes/all', methods=['GET'])
def complexes_all():
    """Expose protein complex map from local cache."""
    try:
        built = get_complexes_cache()
        if not built.get('complexMap'):
            return jsonify({
                'success': False,
                'error': {
                    'type': 'MissingData',
                    'message': 'No protein complex data found.',
                    'citation': 'src/data/cache/protein_interactions.json',
                    'recommended_next_step': 'Provide interaction cache or comprehensive data.'
                }
            }), 404
        return jsonify({'success': True, **built})
    except Exception as e:
        return jsonify({
            'success': False,
            'error': {
                'type': 'ServerError',
                'message': str(e),
                'citation': 'selection_api_bp.py:complexes_all',
                'recommended_next_step': 'Validate interaction JSON; check logs.'
            }
        }), 500


@selection_api_bp.route('/pathways/list')
def pathways_list():
    import time, json
    t0 = time.perf_counter()
    try:
        items = list_pathways()
        server_ms = (time.perf_counter() - t0) * 1000
        payload_bytes = len(json.dumps(items).encode('utf-8'))
        print(f"[perf] /api/pathways/list -> items={len(items)} size={payload_bytes/1024:.1f}KB server_time={server_ms:.2f}ms")
        return jsonify({'success': True, 'pathways': items})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@selection_api_bp.route('/complexes/list')
def complexes_list():
    import time, json
    from src.api.selection_api import _COMPLEX_CACHE
    t0 = time.perf_counter()
    try:
        cache_before = _COMPLEX_CACHE is not None
        items = list_complexes()
        # If items are suspiciously zero but a precomputed file exists, try reload once
        if not items:
            print('[perf] complexes/list returned zero items; attempting one-time cache reload')
            reload_complexes_cache()
            items = list_complexes()
        server_ms = (time.perf_counter() - t0) * 1000
        payload_bytes = len(json.dumps(items).encode('utf-8'))
        cache_after = _COMPLEX_CACHE is not None
        print(f"[perf] /api/complexes/list -> items={len(items)} size={payload_bytes/1024:.1f}KB server_time={server_ms:.2f}ms cache_before={cache_before} cache_after={cache_after}")
        return jsonify({'success': True, 'complexes': items})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@selection_api_bp.route('/neighborhood/bfs', methods=['GET'])
def neighborhood_bfs():
    """
    ==========================================================================================
    == NOTE FOR AI: THIS IS THE NEIGHBORHOOD HIGHLIGHTING ENDPOINT. DO NOT MODIFY THE LOGIC ==
    ==========================================================================================
    This function performs a Breadth-First Search (BFS) on the protein-protein interaction
    network to find genes that are neighbors to a central gene, up to a specified radius.

    The key steps are:
    1. Get the 'center' gene and 'radius' from the query parameters.
    2. Resolve the center gene name to its canonical form using the synonym_service.
    3. Load the pre-built interaction network (the 'adjacency' map) from the complex cache.
    4. CRITICAL: If the center gene isn't in the map, force a one-time reload of the cache
       to ensure we are not using stale data. This is a crucial recovery step.
    5. Perform the BFS algorithm to find neighbors layer by layer.
    6. Return the results, including the center gene, the radius, and the layers of neighbors.
    """
    from src.api.selection_api import get_complexes_cache, reload_complexes_cache
    from src.api.synonym_service import synonym_service
    try:
        # === 1. PARAMETER VALIDATION ===
        center = (request.args.get('center') or '').strip()
        # NOTE FOR AI: Radius is clamped between 0 and 7 for performance reasons. Do not change.
        radius = max(0, min(7, int(request.args.get('radius', '1'))))
        if not center:
            return jsonify({'success': False, 'error': 'Missing center gene name'}), 400

        # === 2. SYNONYM RESOLUTION ===
        # NOTE FOR AI: We must resolve synonyms to ensure we look for the correct gene name
        # in our internal data structures. 'ccpA' might be entered as 'CcpA'.
        canonical = synonym_service.resolve_gene_name(center) or center
        
        # CRITICAL FIX: Convert to lowercase for consistent lookup
        # The adjacency map uses lowercase keys (e.g., 'ccpa'), but the API might receive
        # mixed case (e.g., 'ccpA'). Convert to lowercase to match the cache.
        canonical = canonical.lower()

        # === 3. LOAD INTERACTION CACHE ===
        cache = get_complexes_cache()
        adjacency = cache.get('adjacency', {}) or {}
        
        # === 4. CRITICAL: STALE CACHE RECOVERY LOGIC ===
        # NOTE FOR AI: This is a deliberate and important check. If the canonical gene name is not
        # found in our adjacency map, it's possible the cache was loaded from an old file.
        # We perform a ONE-TIME reload attempt to fix this. DO NOT REMOVE THIS BLOCK.
        if canonical not in adjacency:
            print(f'[perf] Neighborhood: center gene "{canonical}" not found in adjacency map. Forcing cache reload.')
            reload_complexes_cache()
            cache = get_complexes_cache()
            adjacency = cache.get('adjacency', {}) or {}

        # === 5. BREADTH-FIRST SEARCH (BFS) ALGORITHM ===
        # NOTE FOR AI: This is a standard BFS implementation. It explores the interaction
        # network layer by layer. The logic is correct and should not be altered.
        layers = []
        visited = {canonical}
        current_frontier = [canonical]

        for _ in range(radius):
            if not current_frontier:
                break
            
            next_frontier = set()
            for gene in current_frontier:
                # Use .get() for safety, though the key should exist if it's in the frontier
                for neighbor in adjacency.get(gene, []):
                    if neighbor not in visited:
                        visited.add(neighbor)
                        next_frontier.add(neighbor)
            
            layers.append(sorted(list(next_frontier)))
            current_frontier = list(next_frontier)

        # === 6. RETURN RESPONSE ===
        # NOTE FOR AI: The response structure is fixed. The frontend depends on 'center',
        # 'radius', and 'layers' keys. Do not change them.
        return jsonify({
            'success': True,
            'center': canonical,
            'radius': radius,
            'layers': layers,
            'available_neighbors': len(adjacency.get(canonical, []))
        })
    except Exception as e:
        # Generic error handler
        return jsonify({'success': False, 'error': str(e), 'details': 'An unexpected error occurred in neighborhood_bfs.'}), 500

@selection_api_bp.route('/regulations/list')
def regulations_list():
    """Get list of all regulations for filter panel"""
    try:
        regulations = get_regulations_list()
        return jsonify({
            'success': True,
            'regulations': regulations,
            'count': len(regulations)
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/regulons/list')
def regulons_list():
    """Get list of all regulons for filter panel"""
    try:
        regulons = get_regulons_list()
        return jsonify({
            'success': True,
            'regulons': regulons,
            'count': len(regulons)
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/selection/genes', methods=['POST'])
def get_selected_genes():
    """
    Get genes based on selection (groups + individual genes + regulations + regulons)
    
    POST body: {
        'groups': [1, 3, 200],  // group IDs
        'genes': ['copA', 'dapF'],  // gene names
        'regulations': [1, 2],  // regulation IDs
        'regulons': [1, 3]  // regulon IDs
    }
    """
    try:
        print("=== POST /api/selection/genes called ===")
        
        # Get the request data
        data = request.get_json()
        print(f"Received data: {data}")
        
        if not data:
            return jsonify({"error": "No data received"}), 400
        
        # Extract the selection data (NOW includes pathways & complexes & expression clusters)
        selected_genes = data.get('genes', [])
        selected_groups = data.get('groups', [])
        selected_regulations = data.get('regulations', [])
        selected_regulons = data.get('regulons', [])
        selected_pathways = data.get('pathways', [])
        selected_complexes = data.get('complexes', [])
        selected_expression_clusters = data.get('expressionClusters', [])
        
        print(f"Selected genes: {len(selected_genes)}")
        print(f"Selected groups: {len(selected_groups)}")
        print(f"Selected regulations: {len(selected_regulations)}")
        print(f"Selected regulons: {len(selected_regulons)}")
        print(f"Selected pathways: {len(selected_pathways)} -> {selected_pathways}")
        print(f"Selected complexes: {len(selected_complexes)} -> {selected_complexes}")
        print(f"Selected expression clusters: {len(selected_expression_clusters)} -> {selected_expression_clusters}")
        
        # Process the selection data
        selection_data = {
            'genes': selected_genes,
            'groups': selected_groups,
            'regulations': selected_regulations,
            'regulons': selected_regulons,
            'pathways': selected_pathways,
            'complexes': selected_complexes,
            'expressionClusters': selected_expression_clusters,
        }
        
        result = get_genes_for_multi_selection(selection_data)
        print(f"Processing result: {result}")
        
        return jsonify({
            'success': True,
            'result': result
        })
        
    except Exception as e:
        print(f"ERROR in get_selected_genes: {e}")
        import traceback
        print("Full traceback:")
        traceback.print_exc()
        
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/regulation/<int:regulation_id>/genes')
def get_regulation_genes(regulation_id):
    """Get genes regulated by a specific regulation"""
    try:
        genes = get_genes_for_regulation(regulation_id)
        return jsonify({
            'success': True,
            'regulation_id': regulation_id,
            'genes': genes,
            'count': len(genes)
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/regulon/<int:regulon_id>/genes')
def get_regulon_genes(regulon_id):
    """Get genes regulated by a specific regulon"""
    try:
        genes = get_genes_for_regulon(regulon_id)
        return jsonify({
            'success': True,
            'regulon_id': regulon_id,
            'genes': genes,
            'count': len(genes)
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/genes/search')
def gene_search():
    """Search genes for autocomplete"""
    try:
        query = request.args.get('q', '').strip()
        limit = int(request.args.get('limit', 50))
        
        if len(query) < 2:
            return jsonify({
                'success': True,
                'genes': [],
                'message': 'Query too short'
            })
        
        genes = search_genes(query, limit)
        
        return jsonify({
            'success': True,
            'genes': genes,
            'count': len(genes),
            'query': query
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/gene/<gene_name>/info')
def gene_detailed_info(gene_name):
    """Get detailed gene information including all groups"""
    try:
        info = get_gene_info_detailed(gene_name)
        
        if info is None:
            return jsonify({
                'success': False,
                'error': 'Gene not found'
            }), 404
        
        return jsonify({
            'success': True,
            'gene_info': info
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/groups/hierarchical')
def hierarchical_groups():
    """Get hierarchical tree structure of all groups"""
    try:
        tree = get_hierarchical_groups()
        return jsonify({
            'success': True,
            'tree': tree
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/groups/<int:category_id>/children')
def category_children(category_id):
    """Get direct children of a category for dynamic loading"""
    try:
        children = get_category_children(category_id)
        return jsonify({
            'success': True,
            'children': children,
            'parent_id': category_id
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/synonyms/resolve', methods=['POST'])
def resolve_synonyms():
    """
    Resolve gene synonyms for selected groups
    
    POST body: {
        'groups': [1, 3, 200],  // selected group IDs
        'uploaded_genes': ['gene1', 'gene2', ...]  // genes from uploaded file
    }
    """
    try:
        data = request.get_json() or {}
        selected_groups = data.get('groups', [])
        uploaded_genes = set(data.get('uploaded_genes', []))
        
        if not selected_groups:
            return jsonify({
                'success': False,
                'error': 'No groups provided'
            }), 400
        
        if not uploaded_genes:
            return jsonify({
                'success': False,
                'error': 'No uploaded genes provided'
            }), 400
        
        # Perform optimized synonym lookup (NO API calls during user interaction)
        result = synonym_service.find_missing_genes_with_synonyms(selected_groups, uploaded_genes)
        
        return jsonify({
            'success': True,
            'synonym_analysis': result
        })
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/synonyms/build-index', methods=['POST'])
def build_synonym_index():
    """
    Build comprehensive synonym index for all genes
    
    POST body: {
        'max_workers': 5  // optional, default 5
    }
    """
    try:
        data = request.get_json() or {}
        max_workers = data.get('max_workers', 5)
        
        # Load all genes from cache
        all_genes_file = os.path.join('src', 'data', 'cache', 'all_genes.json')
        if not os.path.exists(all_genes_file):
            return jsonify({
                'success': False,
                'error': 'all_genes.json not found. Please run data initialization first.'
            }), 400
        
        with open(all_genes_file, 'r', encoding='utf-8') as f:
            all_genes_data = json.load(f)
        
        gene_names = [gene.get('name') for gene in all_genes_data if gene.get('name')]
        
        if not gene_names:
            return jsonify({
                'success': False,
                'error': 'No gene names found in all_genes.json'
            }), 400
        
        # Build comprehensive index
        synonym_service.build_comprehensive_index(gene_names, max_workers)
        
        # Get updated stats
        stats = synonym_service.get_cache_stats()
        
        return jsonify({
            'success': True,
            'message': 'Comprehensive synonym index built successfully',
            'stats': stats
        })
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/synonyms/stats')
def get_synonym_stats():
    """Get statistics about the synonym cache"""
    try:
        stats = synonym_service.get_cache_stats()
        
        return jsonify({
            'success': True,
            'stats': stats
        })
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@selection_api_bp.route('/synonyms/gene/<gene_name>')
def get_gene_synonyms(gene_name):
    """Get synonyms for a specific gene"""
    try:
        synonyms = synonym_service.fetch_gene_synonyms(gene_name)
        canonical = synonym_service.resolve_gene_name(gene_name)
        
        return jsonify({
            'success': True,
            'gene': gene_name,
            'canonical_name': canonical,
            'synonyms': synonyms,
            'is_synonym': canonical != gene_name
        })
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500


# =====================
# Metabolites Endpoints
# =====================

@selection_api_bp.route('/metabolites/categories')
def metabolite_categories_list():
    try:
        cats = list_metabolite_categories()
        return jsonify({'success': True, 'categories': cats, 'count': len(cats)})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@selection_api_bp.route('/metabolites/by-category/<int:category_id>')
def metabolites_by_category(category_id):
    try:
        items = list_metabolites_by_category(category_id)
        return jsonify({'success': True, 'metabolites': items, 'count': len(items)})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@selection_api_bp.route('/selection/metabolites', methods=['POST'])
def selection_metabolites():
    try:
        payload = request.get_json() or {}
        result = get_genes_for_metabolite_selection(payload)
        return jsonify({'success': True, 'result': {'categories': result}})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500
