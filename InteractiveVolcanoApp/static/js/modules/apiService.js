// apiService.js
import * as state from './state.js';
const API_BASE = '/api';

async function fetchAPI(endpoint, options = {}) {
    const response = await fetch(`${API_BASE}${endpoint}`, options);
    if (!response.ok) {
        const errorText = await response.text();
        throw new Error(`HTTP error! status: ${response.status}, message: ${errorText}`);
    }
    return response.json();
}

export async function loadGroups() {
    const data = await fetchAPI('/groups/main');
    if (data.success) return data.groups;
    throw new Error(data.error || 'Failed to load groups');
}

export async function loadHierarchicalGroups() {
    console.log('📡 API: Starting loadHierarchicalGroups...');
    try {
        const data = await fetchAPI('/groups/hierarchical');
        console.log('📡 API: Received response:', data);
        console.log('📡 API: Response success:', data.success);
        console.log('📡 API: Response tree keys:', Object.keys(data.tree || {}));
        
        if (data.success) {
            console.log('✅ API: Returning tree data successfully');
            return data.tree;
        }
        
        console.error('❌ API: Response not successful:', data.error);
        throw new Error(data.error || 'Failed to load hierarchical groups');
    } catch (error) {
        console.error('❌ API: Exception in loadHierarchicalGroups:', error);
        throw error;
    }
}

export async function loadRegulons() {
    console.log('📡 API REGULONS: Starting loadRegulons...');
    try {
        const data = await fetchAPI('/regulons/list');
        console.log('📡 API REGULONS: Received response:', data);
        console.log('📡 API REGULONS: Response success:', data.success);
        console.log('📡 API REGULONS: Regulons count:', data.regulons?.length);
        
        if (data.success) {
            console.log('✅ API REGULONS: Returning regulons data successfully');
            return data.regulons;
        }
        
        console.error('❌ API REGULONS: Response not successful:', data.error);
        throw new Error(data.error || 'Failed to load regulons');
    } catch (error) {
        console.error('❌ API REGULONS: Exception in loadRegulons:', error);
        throw error;
    }
}

export async function loadExpressionClusters() {
    console.log('📡 API EXPRESSION CLUSTERS: Starting loadExpressionClusters...');
    try {
        const data = await fetchAPI('/expression-clusters/hierarchical');
        console.log('📡 API EXPRESSION CLUSTERS: Received response:', data);
        console.log('📡 API EXPRESSION CLUSTERS: Response success:', data.success);
        console.log('📡 API EXPRESSION CLUSTERS: Clusters keys:', Object.keys(data.clusters || {}));
        
        if (data.success) {
            console.log('✅ API EXPRESSION CLUSTERS: Returning clusters data successfully');
            return data.clusters;
        }
        
        console.error('❌ API EXPRESSION CLUSTERS: Response not successful:', data.error);
        throw new Error(data.error || 'Failed to load expression clusters');
    } catch (error) {
        console.error('❌ API EXPRESSION CLUSTERS: Exception in loadExpressionClusters:', error);
        throw error;
    }
}

export async function searchGenes(query, limit = 50) {
    const data = await fetchAPI(`/genes/search?q=${encodeURIComponent(query)}&limit=${limit}`);
    if (data.success) return data.genes;
    throw new Error(data.error || 'Failed to search genes');
}

export async function applyMultiPanelSelection() {
    const selectionData = {
        groups: Array.from(state.selectedGroups.get()),
        regulons: Array.from(state.selectedRegulons.get()),
        genes: Array.from(state.selectedGenes.get()),
        regulations: Array.from(state.selectedRegulations.get()),
        pathways: Array.from(state.selectedPathways.get ? state.selectedPathways.get() : []),
        complexes: Array.from(state.selectedComplexes.get ? state.selectedComplexes.get() : []),
        expressionClusters: Array.from(state.selectedExpressionClusters.get()),
    };
    console.log('📡 API: Dispatching selection payload:', selectionData);

    const data = await fetchAPI('/selection/genes', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(selectionData),
    });

    console.log('📡 API: Received response:', data);
    if (data.success) {
        // If metabolites were selected, merge their categories into the result
        const selectedMetaboliteIds = Array.from(state.selectedMetabolites.get ? state.selectedMetabolites.get() : []);
        const selectedMetaboliteCategoryIds = Array.from(state.selectedMetaboliteCategories.get ? state.selectedMetaboliteCategories.get() : []);
        if ((selectedMetaboliteIds.length > 0) || (selectedMetaboliteCategoryIds.length > 0)) {
            try {
                const metab = await applyMetaboliteSelection({
                    metabolite_ids: selectedMetaboliteIds.map(Number),
                    metabolite_category_ids: selectedMetaboliteCategoryIds.map(Number),
                    interaction_types: [],
                });
                if (metab && metab.categories) {
                    data.result.categories = {
                        ...data.result.categories,
                        ...metab.categories,
                    };
                }
            } catch (e) {
                console.warn('Metabolite selection merge failed:', e);
            }
        }
        console.log('📡 API: Selection successful, categories:', Object.keys(data.result.categories));
        Object.entries(data.result.categories).forEach(([key, genes]) => {
            console.log(`📡 API: Category ${key}: ${genes ? genes.length : 0} genes`);
        });
        return data;
    }
    throw new Error(data.error || 'Failed to apply multi-panel selection');
}

export async function fetchNeighborhood(centerGene, radius) {
    const q = `/neighborhood/bfs?center=${encodeURIComponent(centerGene)}&radius=${encodeURIComponent(radius)}`;
    const data = await fetchAPI(q);
    if (data.success) return data;
    throw new Error(data.error || 'Failed to fetch neighborhood');
}

// Metabolites API
export async function loadMetaboliteCategories() {
    const data = await fetchAPI('/metabolites/categories');
    if (data.success) return data.categories || [];
    throw new Error(data.error || 'Failed to load metabolite categories');
}

export async function loadMetabolitesByCategory(categoryId) {
    const data = await fetchAPI(`/metabolites/by-category/${encodeURIComponent(categoryId)}`);
    if (data.success) return data.metabolites || [];
    throw new Error(data.error || 'Failed to load metabolites for category');
}

export async function applyMetaboliteSelection(payload) {
    const data = await fetchAPI('/selection/metabolites', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
    });
    if (data.success) return data.result;
    throw new Error(data.error || 'Failed to apply metabolite selection');
}
