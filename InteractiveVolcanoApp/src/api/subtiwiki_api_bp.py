"""
This module defines the blueprint for the SubtiWiki API.
"""

from flask import Blueprint, jsonify
from src.api.subtiwiki_api import (
    get_gene_info,
    get_all_categories,
    get_genes_in_category,
    get_all_genes,
    get_genes_in_category_by_name,
)

subtiwiki_api_bp = Blueprint('subtiwiki_api_bp', __name__)

@subtiwiki_api_bp.route('/gene/<gene_name>')
def gene_info(gene_name):
    """Fetch gene information from the SubtiWiki API."""
    info = get_gene_info(gene_name)
    if info:
        return jsonify(info)
    else:
        return jsonify({'error': 'Gene not found'}), 404

@subtiwiki_api_bp.route('/gene')
def all_genes():
    """Fetch all genes from the SubtiWiki API."""
    genes = get_all_genes()
    if genes:
        return jsonify(genes)
    else:
        return jsonify({'error': 'Could not fetch genes'}), 500

@subtiwiki_api_bp.route('/gene-category')
def all_categories():
    """Fetch all gene categories from the SubtiWiki API."""
    categories = get_all_categories()
    if categories:
        return jsonify(categories)
    else:
        return jsonify({'error': 'Could not fetch categories'}), 500

@subtiwiki_api_bp.route('/gene-category/<category_id>')
def genes_in_category(category_id):
    """Fetch genes in a specific category from the SubtiWiki API."""
    genes = get_genes_in_category(category_id)
    if genes:
        return jsonify(genes)
    else:
        return jsonify({'error': 'Category not found'}), 404

@subtiwiki_api_bp.route('/gene-category/by-name/<path:category_name>')
def genes_in_category_by_name(category_name):
    """Resolve category by name and return gene names (uses cache)."""
    genes = get_genes_in_category_by_name(category_name)
    if genes:
        return jsonify(genes)
    else:
        return jsonify({'error': 'Category not found'}), 404
