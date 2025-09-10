"""
Cache-backed SubtiWiki accessors using local files in src/data/cache.
"""

import os
import json

# Load cached data
CACHE_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'cache')

with open(os.path.join(CACHE_DIR, 'all_genes.json'), 'r', encoding='utf-8') as f:
    ALL_GENES = json.load(f)

with open(os.path.join(CACHE_DIR, 'all_categories.json'), 'r', encoding='utf-8') as f:
    ALL_CATEGORIES = json.load(f)

# Load detailed categories with hierarchical information
with open(os.path.join(CACHE_DIR, 'detailed_categories.json'), 'r', encoding='utf-8') as f:
    DETAILED_CATEGORIES = json.load(f)

with open(os.path.join(CACHE_DIR, 'gene_category_map.json'), 'r', encoding='utf-8') as f:
    GENE_CATEGORY_MAP = json.load(f)

# Create a name-to-id map for genes for quick lookup
GENE_NAME_TO_ID = {gene['name']: gene['id'] for gene in ALL_GENES}

def get_all_genes():
    """Return all genes from cache."""
    return ALL_GENES

def get_all_categories():
    """Return all gene categories from cache."""
    return ALL_CATEGORIES

def get_genes_in_category(category_id):
    """
    Return genes (names) in a category by ID from cache.
    If the category has no direct genes but has children, fetches genes from all subcategories recursively.
    Returns None if category not found.
    """
    # Ensure ID comparison uses the same type as cached IDs
    try:
        category_id_int = int(category_id)
    except (TypeError, ValueError):
        return None

    # Find the detailed category by ID
    category_details = None
    for cat in DETAILED_CATEGORIES:
        if cat['id'] == category_id_int:
            category_details = cat
            break
    
    if not category_details:
        return None
    
    category_name = category_details['name']
    
    # Try to get direct genes for this category
    direct_genes = GENE_CATEGORY_MAP.get(category_name, [])
    
    # If no direct genes but has children, fetch genes from all children recursively
    if not direct_genes and category_details.get('children_ids'):
        child_genes = []
        for child_id in category_details['children_ids']:
            child_genes_result = get_genes_in_category(child_id)  # Recursive call
            if child_genes_result:
                child_genes.extend(child_genes_result)
        
        # Remove duplicates while preserving order
        if child_genes:
            seen = set()
            unique_genes = []
            for gene in child_genes:
                if gene not in seen:
                    seen.add(gene)
                    unique_genes.append(gene)
            return unique_genes
    
    return direct_genes if direct_genes else None

def get_gene_info(gene_name):
    """Return cached gene object by name; None if not found."""
    gene_id = GENE_NAME_TO_ID.get(gene_name)
    if gene_id:
        for gene in ALL_GENES:
            if gene['id'] == gene_id:
                return gene
    return None

def get_category_id_by_name(category_name: str):
    """Resolve a category name to its numeric ID (case-insensitive)."""
    if not category_name:
        return None
    name_norm = category_name.strip().lower()
    for cat in DETAILED_CATEGORIES:
        if isinstance(cat, dict) and cat.get('name', '').strip().lower() == name_norm:
            return cat.get('id')
    return None

def get_genes_in_category_by_name(category_name: str):
    """Return genes (names) by category name; None if not found."""
    cat_id = get_category_id_by_name(category_name)
    if cat_id is None:
        return None
    return get_genes_in_category(cat_id)
