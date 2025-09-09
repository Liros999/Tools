/**
 * Frontend API Module
 *
 * This module handles all communication with the application's Flask backend,
 * which in turn communicates with the Subtiwiki API.
 */

/**
 * Fetches the entire cached interaction graph from the backend.
 * @returns {Promise<Object>} A promise that resolves to the full graph data.
 */
export async function fetchFullGraph() {
    try {
        const response = await fetch('/api/full_graph');
        if (!response.ok) {
            const errorData = await response.json();
            throw new Error(errorData.error || 'Failed to fetch full graph');
        }
        const jsonData = await response.json();
        if (!jsonData.data) {
            throw new Error('Invalid graph data received from server.');
        }
        return jsonData.data;
    } catch (error) {
        console.error('Error fetching full graph:', error);
        throw error;
    }
}

/**
 * Sends the user's data file and settings to the backend for processing.
 * @param {File} file - The user-uploaded CSV/TSV file.
 * @param {Object} settings - An object containing column mappings and filter options.
 * @returns {Promise<Object>} A promise that resolves to the processed data map for coloring.
 */
export async function processFile(file, settings) {
    const formData = new FormData();
    formData.append('file', file);
    formData.append('settings', JSON.stringify(settings));

    try {
        const response = await fetch('/api/process_file', {
            method: 'POST',
            body: formData,
        });

        if (!response.ok) {
            const errorData = await response.json();
            throw new Error(errorData.error || `Server responded with status ${response.status}`);
        }

        return await response.json();

    } catch (error) {
        console.error('Failed to process file:', error);
        throw error;
    }
}

/**
 * Fetches rich context for a specific gene from the backend for tooltips.
 * @param {number} geneId - The Subtiwiki ID of the gene.
 * @returns {Promise<Object>} A promise that resolves with the gene's context.
 */
export async function getGeneContext(geneId) {
    // Simple client-side cache to prevent refetching during the same session hover events
    if (!window.geneContextCache) {
        window.geneContextCache = new Map();
    }
    if (window.geneContextCache.has(geneId)) {
        return window.geneContextCache.get(geneId);
    }

    try {
        // This endpoint needs to be created in the Flask backend
        const response = await fetch(`/api/gene_context/${geneId}`);
        if (!response.ok) {
            throw new Error('Gene context not found.');
        }
        const context = await response.json();
        window.geneContextCache.set(geneId, context);
        return context;
    } catch (error) {
        console.error(`Failed to get context for gene ${geneId}:`, error);
        throw error;
    }
}

/**
 * Searches for genes on the backend.
 * @param {string} query - The search term.
 * @returns {Promise<string[]>} A promise that resolves to an array of matching gene names.
 */
export async function searchGenes(query) {
    try {
        const response = await fetch(`/api/search?query=${encodeURIComponent(query)}`);
        if (!response.ok) {
            throw new Error('Search failed');
        }
        return await response.json();
    } catch (error) {
        console.error('Gene search error:', error);
        return []; // Return empty array on error
    }
}
