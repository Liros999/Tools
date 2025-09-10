// uiManager.js
import * as state from './state.js';
import * as api from './apiService.js';
import * as plotManager from './plotManager.js';
import { pathwayInfo, selectedPathways } from './state.js';

export function setupEventHandlers() {
    const applyButton = document.getElementById('apply-selection');
    if (applyButton) {
        applyButton.addEventListener('click', handleApplyMultiPanelSelection);
    }
    
    // Setup tab handlers for regulons and other tabs
    setupTabHandlers();
    
    // Set up regulon search
    setupRegulonSearch();
    
    // Setup plot control buttons
    setupPlotControlButtons();
    
    // Setup clear all button
    setupClearAllButton();
}

// Coloring mode removed; integrated filters
export function initializeColoringUI() { /* no-op for integrated filters */ }

function attachSearch(inputId, listRootSelector) {
    const input = document.getElementById(inputId);
    if (!input) return;
    input.addEventListener('input', () => {
        const q = input.value.toLowerCase();
        document.querySelectorAll(`${listRootSelector} .list-group-item`).forEach(li => {
            li.style.display = li.textContent.toLowerCase().includes(q) ? '' : 'none';
        });
    });
}

export function buildPathwayFilter() {
    const list = document.getElementById('pathway-selection-list') || document.getElementById('pathway-list');
    const container = document.getElementById('pathways-panel') || document.getElementById('pathway-filter-container');
    if (!list || !container) return;
    const info = pathwayInfo.get() || {};
    const entries = Object.entries(info);
    if (entries.length === 0) { return; }
    list.innerHTML = '';
    entries.sort((a,b)=> (a[1]?.name||'').localeCompare(b[1]?.name||''))
        .forEach(([pid, meta]) => {
            const id = `pathway-${pid}`;
            const div = document.createElement('div');
            div.className = 'list-group-item';
            div.innerHTML = `
              <input class="form-check-input me-1" type="checkbox" value="${pid}" id="${id}">
              <label class="form-check-label stretched-link" for="${id}">${meta?.name || pid}</label>
            `;
            list.appendChild(div);
            div.querySelector('input').addEventListener('change', ev => {
                const v = ev.target.value;
                if (ev.target.checked) selectedPathways.add(v); else selectedPathways.delete(v);
                // Unified UX: defer plot update until user clicks Apply
            });
        });
    attachSearch('pathway-search', '#pathway-selection-list');
}

// Complexes UI removed per updated design

async function loadPathwaysPanel() {
    try {
        const list = document.getElementById('pathway-selection-list');
        if (!list) return;
        if (list.children.length > 0 && !list.querySelector('.text-muted')) return;
        const info = pathwayInfo.get() || {};
        if (Object.keys(info).length > 0) { buildPathwayFilter(); return; }
        const r = await fetch('/api/pathways/list');
        if (r.ok) {
            const d = await r.json();
            if (d.success && Array.isArray(d.pathways)) {
                const infoMap = {};
                d.pathways.forEach(p => { if (p && p.id != null) infoMap[String(p.id)] = { name: p.name || `Pathway ${p.id}` }; });
                pathwayInfo.set(infoMap);
                buildPathwayFilter();
            }
        }
    } catch (e) { console.warn('Pathways panel load failed:', e); }
}

async function loadExpressionClustersPanel() {
    console.log('🧬 EXPRESSION CLUSTERS: Starting loadExpressionClustersPanel...');
    const container = document.getElementById('expression-cluster-selection-tree');
    if (!container) {
        console.error('❌ EXPRESSION CLUSTERS: Expression cluster selection tree container not found');
        return;
    }
    
    console.log('✅ EXPRESSION CLUSTERS: Container found, loading expression clusters...');
    
    try {
        container.innerHTML = '<div class="text-center p-3"><i class="fas fa-spinner fa-spin"></i> Loading expression clusters...</div>';
        
        console.log('📡 EXPRESSION CLUSTERS: Calling api.loadExpressionClusters()...');
        const clusters = await api.loadExpressionClusters();
        console.log(`📋 EXPRESSION CLUSTERS: Loaded clusters from API`);
        console.log('🔍 EXPRESSION CLUSTERS: Cluster types:', Object.keys(clusters));
        
        // Render hierarchical tree structure with collapsible parent groups
        let html = '';
        
        // Render A clusters
        if (clusters.A && clusters.A.length > 0) {
            html += '<div class="cluster-category mb-2">';
            html += '<div class="cluster-parent fw-bold text-primary" style="cursor: pointer;" onclick="toggleClusterCategory(\'A\')">';
            html += '<i class="fas fa-chevron-down me-2" id="toggle-A"></i>A Clusters</div>';
            html += '<div class="cluster-children" id="children-A">';
            clusters.A.forEach(cluster => {
                html += `
                    <div class="form-check ms-3 mb-1">
                        <input class="form-check-input expression-cluster-checkbox" type="checkbox" 
                               data-id="${cluster.id}" data-name="${cluster.name}" 
                               data-type="A" id="cluster-${cluster.id}">
                        <label class="form-check-label small" for="cluster-${cluster.id}">
                            <strong>${cluster.name}</strong>
                            <span class="text-muted"> (${cluster.gene_count} genes)</span>
                        </label>
                    </div>`;
            });
            html += '</div></div>';
        }
        
        // Render B clusters
        if (clusters.B && clusters.B.length > 0) {
            html += '<div class="cluster-category mb-2">';
            html += '<div class="cluster-parent fw-bold text-success" style="cursor: pointer;" onclick="toggleClusterCategory(\'B\')">';
            html += '<i class="fas fa-chevron-down me-2" id="toggle-B"></i>B Clusters</div>';
            html += '<div class="cluster-children" id="children-B">';
            clusters.B.forEach(cluster => {
                html += `
                    <div class="form-check ms-3 mb-1">
                        <input class="form-check-input expression-cluster-checkbox" type="checkbox" 
                               data-id="${cluster.id}" data-name="${cluster.name}" 
                               data-type="B" id="cluster-${cluster.id}">
                        <label class="form-check-label small" for="cluster-${cluster.id}">
                            <strong>${cluster.name}</strong>
                            <span class="text-muted"> (${cluster.gene_count} genes)</span>
                        </label>
                    </div>`;
            });
            html += '</div></div>';
        }
        
        // Render C clusters
        if (clusters.C && clusters.C.length > 0) {
            html += '<div class="cluster-category mb-2">';
            html += '<div class="cluster-parent fw-bold text-warning" style="cursor: pointer;" onclick="toggleClusterCategory(\'C\')">';
            html += '<i class="fas fa-chevron-down me-2" id="toggle-C"></i>C Clusters</div>';
            html += '<div class="cluster-children" id="children-C">';
            clusters.C.forEach(cluster => {
                html += `
                    <div class="form-check ms-3 mb-1">
                        <input class="form-check-input expression-cluster-checkbox" type="checkbox" 
                               data-id="${cluster.id}" data-name="${cluster.name}" 
                               data-type="C" id="cluster-${cluster.id}">
                        <label class="form-check-label small" for="cluster-${cluster.id}">
                            <strong>${cluster.name}</strong>
                            <span class="text-muted"> (${cluster.gene_count} genes)</span>
                        </label>
                    </div>`;
            });
            html += '</div></div>';
        }
        
        console.log(`🏗️ EXPRESSION CLUSTERS: Generated HTML length: ${html.length}`);
        container.innerHTML = html;
        setupExpressionClusterHandlers(container);
        setupExpressionClusterSearch();
        console.log('✅ EXPRESSION CLUSTERS: Panel loaded successfully');
        
    } catch (error) {
        console.error('❌ EXPRESSION CLUSTERS: Error loading expression clusters:', error);
        console.error('❌ EXPRESSION CLUSTERS: Error stack:', error.stack);
        container.innerHTML = '<div class="text-danger p-3">Error loading expression clusters: ' + error.message + '</div>';
    }
}

// Function to toggle cluster category visibility (expand/collapse)
function toggleClusterCategory(category) {
    const childrenContainer = document.getElementById(`children-${category}`);
    const toggleIcon = document.getElementById(`toggle-${category}`);
    
    if (!childrenContainer || !toggleIcon) {
        console.error(`❌ Toggle elements not found for category ${category}`);
        return;
    }
    
    const isExpanded = childrenContainer.style.display !== 'none';
    
    if (isExpanded) {
        // Collapse
        childrenContainer.style.display = 'none';
        toggleIcon.classList.remove('fa-chevron-down');
        toggleIcon.classList.add('fa-chevron-right');
        console.log(`🔼 Collapsed ${category} clusters`);
    } else {
        // Expand
        childrenContainer.style.display = 'block';
        toggleIcon.classList.remove('fa-chevron-right');
        toggleIcon.classList.add('fa-chevron-down');
        console.log(`🔽 Expanded ${category} clusters`);
    }
}

// Make the function globally accessible for onclick handlers
window.toggleClusterCategory = toggleClusterCategory;

// Complexes loading removed per updated design

function setupPlotControlButtons() {
    // Labels toggle button
    const labelsButton = document.getElementById('toggle-labels');
    if (labelsButton) {
        labelsButton.addEventListener('click', function() {
            console.log('🏷️ Labels button clicked');
            plotManager.toggleGeneLabels();
            // Button state is now managed within toggleGeneLabels function
        });
    }
    
    // Download PNG button
    const downloadButton = document.getElementById('download-png');
    if (downloadButton) {
        downloadButton.addEventListener('click', function() {
            console.log('💾 Download PNG button clicked');
            plotManager.downloadPlotAsPNG();
        });
    }
    
    // Reset axes button
    const resetButton = document.getElementById('reset-axes');
    if (resetButton) {
        resetButton.addEventListener('click', function() {
            console.log('🔄 Reset axes button clicked');
            plotManager.resetPlotScale();
        });
    }
}

// Clear All button functionality
function setupClearAllButton() {
    const clearButton = document.getElementById('clear-selection');
    if (!clearButton) {
        console.error('❌ Clear All button not found');
        return;
    }
    
    clearButton.addEventListener('click', function() {
        console.log('🧹 Clear All button clicked - clearing everything');
        
        // Clear all state selections
        state.selectedGroups.clear();
        state.selectedRegulons.clear();
        state.selectedGenes.clear();
        state.selectedRegulations.clear();
        state.selectedExpressionClusters.clear();
        // NEW: clear pathways and complexes as part of global reset
        if (state.selectedPathways && state.selectedPathways.clear) state.selectedPathways.clear();
        if (state.selectedComplexes && state.selectedComplexes.clear) state.selectedComplexes.clear();
        
        // Clear all checkboxes in Groups panel
        const groupCheckboxes = document.querySelectorAll('#group-selection-tree input[type="checkbox"]');
        groupCheckboxes.forEach(cb => cb.checked = false);
        
        // Clear all checkboxes in Regulons panel
        const regulonCheckboxes = document.querySelectorAll('#regulon-selection-list input[type="checkbox"]');
        regulonCheckboxes.forEach(cb => cb.checked = false);

        // NEW: clear any Pathway/Complex checkboxes if panels exist
        const pathwayCheckboxes = document.querySelectorAll('#pathway-selection-list input[type="checkbox"], #pathway-list input[type="checkbox"]');
        pathwayCheckboxes.forEach(cb => cb.checked = false);
        const complexCheckboxes = document.querySelectorAll('#complex-selection-list input[type="checkbox"], #complex-list input[type="checkbox"]');
        complexCheckboxes.forEach(cb => cb.checked = false);
        
        // Clear expression cluster checkboxes
        const expressionClusterCheckboxes = document.querySelectorAll('#expression-cluster-selection-tree input[type="checkbox"]');
        expressionClusterCheckboxes.forEach(cb => cb.checked = false);
        
        // Clear gene search input
        const geneSearchInput = document.getElementById('gene-search');
        if (geneSearchInput) geneSearchInput.value = '';
        
        // Clear regulon search input
        const regulonSearchInput = document.getElementById('regulon-search');
        if (regulonSearchInput) regulonSearchInput.value = '';
        // NEW: clear pathway/complex search inputs if present
        const pathwaySearchInput = document.getElementById('pathway-search');
        if (pathwaySearchInput) pathwaySearchInput.value = '';
        const complexSearchInput = document.getElementById('complex-search');
        if (complexSearchInput) complexSearchInput.value = '';
        
        // Clear expression cluster search input
        const expressionClusterSearchInput = document.getElementById('expression-cluster-search');
        if (expressionClusterSearchInput) expressionClusterSearchInput.value = '';
        
        // Reset plot to default coloring (all grey)
        plotManager.resetPlotColors();
        
        // Remove any labels
        plotManager.removeGeneLabels();
        
        // Reset labels button to inactive state
        const labelsButton = document.getElementById('toggle-labels');
        if (labelsButton) {
            labelsButton.classList.remove('btn-primary');
            labelsButton.classList.add('btn-outline-primary');
        }
        
        console.log('✅ All selections and visualizations cleared');
    });
}

// Regulon search functionality
function setupRegulonSearch() {
    const searchInput = document.getElementById('regulon-search');
    if (!searchInput) {
        console.error('❌ Regulon search input not found');
        return;
    }
    
    searchInput.addEventListener('input', function(e) {
        const searchTerm = e.target.value.toLowerCase().trim();
        console.log('🔍 Searching regulons for:', searchTerm);
        
        const regulonCheckboxes = document.querySelectorAll('#regulon-selection-list .form-check');
        regulonCheckboxes.forEach(checkbox => {
            const label = checkbox.querySelector('label');
            if (label) {
                const text = label.textContent.toLowerCase();
                if (searchTerm === '' || text.includes(searchTerm)) {
                    checkbox.style.display = '';
                } else {
                    checkbox.style.display = 'none';
                }
            }
        });
    });
}

function setupTabHandlers() {
    const regulonsTab = document.getElementById('regulons-tab');
    if (regulonsTab) {
        regulonsTab.addEventListener('click', function() {
            console.log('🧬 Regulons tab clicked - loading regulons...');
            loadRegulonsPanel();
        });
    }
    
    const groupsTab = document.getElementById('groups-tab');
    if (groupsTab) {
        groupsTab.addEventListener('click', function() {
            console.log('📁 Groups tab clicked - ensuring groups are loaded...');
            // Groups are already loaded in renderHierarchicalTree
        });
    }

    const pathwaysTab = document.getElementById('pathways-tab');
    if (pathwaysTab) {
        pathwaysTab.addEventListener('click', function() {
            console.log('🧪 Pathways tab clicked - loading pathways list...');
            loadPathwaysPanel();
        });
    }

    const expressionClustersTab = document.getElementById('expression-clusters-tab');
    if (expressionClustersTab) {
        expressionClustersTab.addEventListener('click', function() {
            console.log('🧬 Expression Clusters tab clicked - loading expression clusters...');
            loadExpressionClustersPanel();
        });
    }

    // Complexes tab removed

    const metabolitesTab = document.getElementById('metabolites-tab');
    if (metabolitesTab) {
        metabolitesTab.addEventListener('click', function() {
            console.log('🧪 Metabolites tab clicked - loading metabolite categories...');
            loadMetabolitesPanel();
        });
    }
}

async function loadRegulonsPanel() {
    console.log('🧬 REGULONS: Starting loadRegulonsPanel...');
    const container = document.getElementById('regulon-selection-list');
    if (!container) {
        console.error('❌ REGULONS: Regulon selection list container not found');
        return;
    }
    
    console.log('✅ REGULONS: Container found, loading regulons...');
    
    try {
        container.innerHTML = '<div class="text-center p-3"><i class="fas fa-spinner fa-spin"></i> Loading regulons...</div>';
        
        console.log('📡 REGULONS: Calling api.loadRegulons()...');
        const regulons = await api.loadRegulons();
        console.log(`📋 REGULONS: Loaded ${regulons.length} regulons from API`);
        console.log('🔍 REGULONS: Sample regulon:', regulons[0]);
        
        let html = '';
        regulons.forEach(regulon => {
            html += `
                <div class="form-check mb-2">
                    <input class="form-check-input regulon-checkbox" type="checkbox" 
                           data-id="${regulon.id}" data-name="${regulon.name}" 
                           id="regulon-${regulon.id}">
                    <label class="form-check-label small" for="regulon-${regulon.id}">
                        <strong>${regulon.name}</strong>
                        <br><span class="text-muted">${regulon.gene_count || 0} genes</span>
                        ${regulon.description ? `<br><small class="text-secondary">${regulon.description}</small>` : ''}
                    </label>
                </div>`;
        });
        
        console.log(`🏗️ REGULONS: Generated HTML length: ${html.length}`);
        container.innerHTML = html;
        setupRegulonHandlers(container);
        console.log('✅ REGULONS: Panel loaded successfully');
        
    } catch (error) {
        console.error('❌ REGULONS: Error loading regulons:', error);
        console.error('❌ REGULONS: Error stack:', error.stack);
        
        let errorMessage = 'Error loading regulons: ' + error.message;
        
        // Provide more specific error messages based on error type
        if (error.message.includes('Failed to fetch')) {
            errorMessage = 'Network error: Unable to connect to the server. Please check your internet connection and try again.';
        } else if (error.message.includes('timeout')) {
            errorMessage = 'Request timeout: The server took too long to respond. Please try again.';
        } else if (error.message.includes('404')) {
            errorMessage = 'Service not found: The regulon service is currently unavailable.';
        }
        
        container.innerHTML = '<div class="text-danger p-3">' + errorMessage + '</div>';
    }
}

function setupRegulonHandlers(container) {
    const checkboxes = container.querySelectorAll('.regulon-checkbox');
    checkboxes.forEach(checkbox => {
        checkbox.addEventListener('change', function() {
            const regulonId = parseInt(this.dataset.id);
            const regulonName = this.dataset.name;
            
            if (this.checked) {
                state.selectedRegulons.add(regulonId);
                console.log(`✓ Selected regulon: ${regulonName} (ID: ${regulonId})`);
            } else {
                state.selectedRegulons.delete(regulonId);
                console.log(`✗ Deselected regulon: ${regulonName} (ID: ${regulonId})`);
            }
            
            console.log('Current selected regulons:', Array.from(state.selectedRegulons.get()));
        });
    });
}

function setupExpressionClusterHandlers(container) {
    const checkboxes = container.querySelectorAll('.expression-cluster-checkbox');
    checkboxes.forEach(checkbox => {
        checkbox.addEventListener('change', function() {
            const clusterId = this.dataset.id;
            const clusterName = this.dataset.name;
            const clusterType = this.dataset.type;
            
            if (this.checked) {
                state.selectedExpressionClusters.add(clusterId);
                console.log(`✓ Selected expression cluster: ${clusterName} (Type: ${clusterType})`);
            } else {
                state.selectedExpressionClusters.delete(clusterId);
                console.log(`✗ Deselected expression cluster: ${clusterName} (Type: ${clusterType})`);
            }
            
            console.log('Current selected expression clusters:', Array.from(state.selectedExpressionClusters.get()));
        });
    });
}

function setupExpressionClusterSearch() {
    const searchInput = document.getElementById('expression-cluster-search');
    if (!searchInput) {
        console.error('❌ Expression cluster search input not found');
        return;
    }
    
    searchInput.addEventListener('input', function(e) {
        const searchTerm = e.target.value.toLowerCase().trim();
        console.log('🔍 Searching expression clusters for:', searchTerm);
        
        const clusterCheckboxes = document.querySelectorAll('#expression-cluster-selection-tree .form-check');
        clusterCheckboxes.forEach(checkbox => {
            const label = checkbox.querySelector('label');
            if (label) {
                const text = label.textContent.toLowerCase();
                if (searchTerm === '' || text.includes(searchTerm)) {
                    checkbox.style.display = '';
                } else {
                    checkbox.style.display = 'none';
                }
            }
        });
    });
}

export function setupGeneSearch() {
    console.log('Setting up gene search...');
    const searchInput = document.getElementById('gene-search');
    if (!searchInput) {
        console.error('Gene search input not found');
        return;
    }

    // Create results container if it doesn't exist
    let resultsContainer = document.getElementById('gene-search-results');
    if (!resultsContainer) {
        resultsContainer = document.createElement('div');
        resultsContainer.id = 'gene-search-results';
        resultsContainer.className = 'dropdown-menu';
        resultsContainer.style.cssText = 'position: absolute; top: 100%; left: 0; right: 0; z-index: 1000; display: none; max-height: 200px; overflow-y: auto;';
        
        // Insert after the search input
        searchInput.parentNode.style.position = 'relative';
        searchInput.parentNode.appendChild(resultsContainer);
    }

    let searchTimeout;
    searchInput.addEventListener('input', function(e) {
        const query = e.target.value.trim();
        
        // Clear previous timeout
        if (searchTimeout) {
            clearTimeout(searchTimeout);
        }
        
        if (query.length < 2) {
            hideSearchResults();
            return;
        }
        
        // Debounce search
        searchTimeout = setTimeout(async () => {
            try {
                const genes = await api.searchGenes(query);
                displayGeneSearchResults(genes, resultsContainer, searchInput);
            } catch (error) {
                console.error('Gene search error:', error);
                hideSearchResults();
            }
        }, 300);
    });

    // Hide results when clicking outside
    document.addEventListener('click', function(e) {
        if (!searchInput.contains(e.target) && !resultsContainer.contains(e.target)) {
            hideSearchResults();
        }
    });

    function hideSearchResults() {
        resultsContainer.style.display = 'none';
    }

    function displayGeneSearchResults(genes, container, input) {
        container.innerHTML = '';
        
        if (!genes || genes.length === 0) {
            container.style.display = 'none';
            return;
        }
        
        genes.slice(0, 10).forEach(gene => {
            const item = document.createElement('div');
            item.className = 'dropdown-item';
            item.style.cssText = 'cursor: pointer; padding: 8px 12px; border-bottom: 1px solid #eee;';
            item.textContent = gene;
            
            item.addEventListener('click', function() {
                addGeneToSelection(gene);
                input.value = gene; // keep gene for neighborhood center
                container.style.display = 'none';
            });
            
            item.addEventListener('mouseenter', function() {
                this.style.backgroundColor = '#f8f9fa';
            });
            
            item.addEventListener('mouseleave', function() {
                this.style.backgroundColor = 'white';
            });
            
            container.appendChild(item);
        });
        
        container.style.display = 'block';
    }
}

// Wire persistent neighborhood controls
document.addEventListener('DOMContentLoaded', () => {
    const depthSlider = document.getElementById('neighborhood-depth');
    const depthValue = document.getElementById('neighborhood-depth-value');
    const geneInput = document.getElementById('gene-search');
    if (depthSlider && depthValue && geneInput) {
        const handler = async () => {
            depthValue.textContent = String(depthSlider.value);
            const center = (geneInput.value || '').trim();
            if (!center) return; // Wait until a gene is entered
            const r = Math.max(0, Math.min(7, parseInt(depthSlider.value || '1', 10)));
            try {
                const data = await api.fetchNeighborhood(center, r);
                plotManager.applyNeighborhoodColoring(data.center, data.layers);
            } catch (e) {
                console.error('Neighborhood fetch failed:', e);
            }
        };
        depthSlider.addEventListener('input', handler);
        depthSlider.addEventListener('change', handler);
    }
});

function addGeneToSelection(geneName) {
    state.selectedGenes.add(geneName);
    updateSelectedGenesDisplay();
}

function updateSelectedGenesDisplay() {
    const container = document.getElementById('selected-genes');
    if (!container) return;
    
    container.innerHTML = '';
    const selectedGenes = Array.from(state.selectedGenes.get());
    
    selectedGenes.forEach(gene => {
        const tag = document.createElement('span');
        tag.className = 'badge bg-primary me-1 mb-1';
        tag.textContent = gene;
        
        const removeBtn = document.createElement('button');
        removeBtn.type = 'button';
        removeBtn.className = 'btn-close btn-close-white ms-1';
        removeBtn.style.fontSize = '0.5em';
        removeBtn.addEventListener('click', () => {
            state.selectedGenes.delete(gene);
            updateSelectedGenesDisplay();
        });
        
        tag.appendChild(removeBtn);
        container.appendChild(tag);
    });
}

export function renderHierarchicalTree() {
    console.log('🌳 STARTING renderHierarchicalTree...');
    const container = document.getElementById('group-selection-tree');
    if (!container) {
        console.error('❌ group-selection-tree element not found');
        return;
    }
    
    console.log('📡 Calling api.loadHierarchicalGroups()...');
    container.innerHTML = '<div class="text-center p-3"><i class="fas fa-spinner fa-spin"></i> Loading groups...</div>';
    
    api.loadHierarchicalGroups()
        .then(treeData => {
            console.log('✅ Received tree data:', treeData);
            console.log('🔍 Tree data type:', typeof treeData);
            console.log('🔍 Tree data keys:', Object.keys(treeData || {}));
            
            const htmlResult = buildTreeHTML(treeData);
            console.log('🏗️ Built HTML length:', htmlResult.length);
            
            container.innerHTML = htmlResult;
            setupTreeToggleHandlers(container);
            console.log('🎯 Tree rendering completed successfully');
        })
        .catch(error => {
            console.error('❌ CRITICAL ERROR in renderHierarchicalTree:', error);
            console.error('❌ Error stack:', error.stack);
            container.innerHTML = '<div class="text-danger p-3">Error loading groups: ' + error.message + '</div>';
        });
}

function setupTreeToggleHandlers(container) {
    // Handle checkbox changes for group selection
    const checkboxes = container.querySelectorAll('.group-checkbox');
    checkboxes.forEach(checkbox => {
        checkbox.addEventListener('change', function() {
            const groupId = parseInt(this.dataset.id);
            const groupName = this.dataset.name;
            
            if (this.checked) {
                state.selectedGroups.add(groupId);
                console.log(`✓ Selected group: ${groupName} (ID: ${groupId})`);
            } else {
                state.selectedGroups.delete(groupId);
                console.log(`✗ Deselected group: ${groupName} (ID: ${groupId})`);
            }
            
            console.log('Current selected groups:', Array.from(state.selectedGroups.get()));
        });
    });
    
    // Handle tree expand/collapse functionality
    const toggles = container.querySelectorAll('.tree-toggle');
    toggles.forEach(toggle => {
        toggle.addEventListener('click', function(e) {
            e.stopPropagation();
            const nodeId = this.dataset.id;
            const childrenContainer = container.querySelector(`.tree-children[data-parent-id="${nodeId}"]`);
            
            if (childrenContainer) {
                const isCollapsed = childrenContainer.style.display === 'none';
                
                if (isCollapsed) {
                    // Expand
                    childrenContainer.style.display = 'block';
                    this.style.transform = 'rotate(90deg)';
                    this.classList.remove('fa-chevron-right');
                    this.classList.add('fa-chevron-down');
                    console.log(`🔽 Expanded group node: ${nodeId}`);
    } else {
                    // Collapse
                    childrenContainer.style.display = 'none';
                    this.style.transform = 'rotate(0deg)';
                    this.classList.remove('fa-chevron-down');
                    this.classList.add('fa-chevron-right');
                    console.log(`🔼 Collapsed group node: ${nodeId}`);
                }
            }
        });
    });
    
    // Handle label clicks for expand/collapse
    const labels = container.querySelectorAll('.tree-label');
    labels.forEach(label => {
        label.addEventListener('click', function(e) {
            const nodeElement = this.closest('.tree-node');
            const toggle = nodeElement.querySelector('.tree-toggle');
            if (toggle) {
                toggle.click();
            }
        });
    });
}

// HIERARCHICAL TREE BUILDER - FIXED FOR ACTUAL DATA STRUCTURE
function buildTreeHTML(data) {
    console.log('🏗️ BUILDING TREE: Starting buildTreeHTML with data:', data);
    console.log('🏗️ BUILDING TREE: Data type:', typeof data);
    console.log('🏗️ BUILDING TREE: Data keys count:', Object.keys(data || {}).length);
    
    if (!data || typeof data !== 'object') {
        console.error('❌ BUILDING TREE: Invalid data received');
        return '<div class="text-danger p-3">Invalid tree data received</div>';
    }
    
    let html = '';
    
    function buildNodeHTML(node, level = 0) {
        const indentPx = level * 20;
        const hasChildren = node.children && Object.keys(node.children).length > 0;
        
        let nodeHTML = `
            <div class="tree-node" data-id="${node.id}" data-level="${level}" data-dot-notation="${node.dot_notation}">
                <div class="d-flex align-items-center" style="margin-left: ${indentPx}px;">
                    ${hasChildren ? 
                        `<i class="fas fa-chevron-right tree-toggle me-2" data-id="${node.id}" style="cursor: pointer; font-size: 12px; transition: transform 0.2s;"></i>` : 
                        `<span class="me-2" style="width: 16px;"></span>`
                    }
                    <input type="checkbox" class="form-check-input me-2 group-checkbox" data-id="${node.id}" data-name="${node.name}">
                    <span class="tree-label" style="cursor: pointer; user-select: none;">${node.name} (${node.gene_count || 0})</span>
                </div>
                ${hasChildren ? `<div class="tree-children" data-parent-id="${node.id}" style="display: none;"></div>` : ''}
            </div>`;
            
        // Add children if they exist
        if (hasChildren) {
            let childrenHTML = '';
            Object.values(node.children).forEach(child => {
                childrenHTML += buildNodeHTML(child, level + 1);
            });
            nodeHTML = nodeHTML.replace(`<div class="tree-children" data-parent-id="${node.id}" style="display: none;"></div>`, 
                `<div class="tree-children" data-parent-id="${node.id}" style="display: none;">${childrenHTML}</div>`);
        }
        
        return nodeHTML;
    }
    
    // Build HTML for all top-level nodes
    Object.values(data).forEach(node => {
        html += buildNodeHTML(node);
    });
    
    console.log('🏗️ BUILDING TREE: Generated HTML length:', html.length);
    return html;
}

async function handleApplyMultiPanelSelection() {
    console.log('Applying selection...');
    console.log('📋 Currently selected groups:', Array.from(state.selectedGroups.get()));
    console.log('📋 Currently selected regulons:', Array.from(state.selectedRegulons.get()));
    console.log('📋 Currently selected genes:', Array.from(state.selectedGenes.get()));
    console.log('📋 Currently selected expression clusters:', Array.from(state.selectedExpressionClusters.get()));
    
    try {
        const result = await api.applyMultiPanelSelection();
        if (result && result.result && result.result.categories) {
            console.log('API call succeeded. Updating plot with full result object.');
            console.log('📊 Categories:', Object.keys(result.result.categories));
            console.log('🏷️ Category names:', result.result.category_names);
            
            // Check for expression clusters with no matching genes
            const expressionClusters = result.result.categories || {};
            const emptyClusters = [];
            for (const [key, genes] of Object.entries(expressionClusters)) {
                if (key.startsWith('expression_cluster_') && genes.length === 0) {
                    const clusterId = key.replace('expression_cluster_', '');
                    emptyClusters.push(clusterId);
                }
            }
            
            if (emptyClusters.length > 0) {
                console.warn(`⚠️ Warning: ${emptyClusters.length} expression cluster(s) have no genes in the current plot: ${emptyClusters.join(', ')}`);
                // You can add a user notification here if you have a notification system
            }
            
            plotManager.updatePlotWithSelection(result.result);  // Pass full result object
        } else {
            console.error('API call succeeded but returned no categories.');
        }
    } catch (error) {
        console.error('Error during apply selection:', error);
    }
}

// =====================
// Metabolites UI
// =====================
async function loadMetabolitesPanel() {
    const container = document.getElementById('metabolite-category-tree');
    if (!container) return;
    try {
        container.innerHTML = '<div class="text-center p-3"><i class="fas fa-spinner fa-spin"></i> Loading metabolite categories...</div>';
        const cats = await api.loadMetaboliteCategories();
        if (!Array.isArray(cats) || cats.length === 0) {
            container.innerHTML = '<div class="text-muted p-3">No metabolite categories available.</div>';
            return;
        }
        const tree = buildMetaboliteTree(cats);
        container.innerHTML = buildMetaboliteTreeHTML(tree);
        setupMetaboliteTreeHandlers(container);
        // Refresh metabolite list if some categories were previously selected
        await refreshSelectedMetabolitesList();
    } catch (e) {
        console.error('Metabolite categories load failed:', e);
        container.innerHTML = '<div class="text-danger p-3">Failed to load metabolite categories</div>';
    }
}

function buildMetaboliteTree(cats) {
    // Build lookup by id
    const nodesById = new Map();
    cats.forEach(c => nodesById.set(String(c.id), { id: c.id, name: c.name, dot: c.dot_notation || '', parent_id: c.parent_id, children: [] }));
    // Attach children
    const roots = {};
    nodesById.forEach(node => {
        const pid = node.parent_id;
        if (pid != null && nodesById.has(String(pid))) {
            nodesById.get(String(pid)).children.push(node);
        } else {
            roots[node.dot || String(node.id)] = node;
        }
    });
    return roots;
}

function buildMetaboliteTreeHTML(tree) {
    function render(node, level) {
        const hasChildren = node.children && node.children.length > 0;
        const indentPx = level * 20;
        const childrenHTML = hasChildren ? node.children.map(ch => render(ch, level + 1)).join('') : '';
        return `
            <div class="tree-node" data-id="${node.id}" data-level="${level}">
                <div class="d-flex align-items-center" style="margin-left: ${indentPx}px;">
                    ${hasChildren ? `<i class="fas fa-chevron-right tree-toggle me-2" data-id="${node.id}" style="cursor: pointer; font-size: 12px; transition: transform 0.2s;"></i>` : `<span class="me-2" style="width: 16px;"></span>`}
                    <input type="checkbox" class="form-check-input me-2 metabolite-group-checkbox" data-id="${node.id}" data-name="${node.name}">
                    <span class="tree-label" style="cursor: pointer; user-select: none;">${node.name}</span>
                </div>
                ${hasChildren ? `<div class="tree-children" data-parent-id="${node.id}" style="display: none;">${childrenHTML}</div>` : ''}
            </div>`;
    }
    return Object.values(tree).map(root => render(root, 0)).join('');
}

function setupMetaboliteTreeHandlers(container) {
    // Multi-select categories
    const boxes = container.querySelectorAll('.metabolite-group-checkbox');
    boxes.forEach(cb => {
        cb.addEventListener('change', async function() {
            const id = parseInt(this.dataset.id);
            if (this.checked) state.selectedMetaboliteCategories.add(id); else state.selectedMetaboliteCategories.delete(id);
            await refreshSelectedMetabolitesList();
        });
    });
    // Expand/collapse
    const toggles = container.querySelectorAll('.tree-toggle');
    toggles.forEach(t => {
        t.addEventListener('click', function(e) {
            e.stopPropagation();
            const pid = this.dataset.id;
            const children = container.querySelector(`.tree-children[data-parent-id="${pid}"]`);
            if (!children) return;
            const isCollapsed = children.style.display === 'none';
            children.style.display = isCollapsed ? 'block' : 'none';
            this.classList.toggle('fa-chevron-right', !isCollapsed);
            this.classList.toggle('fa-chevron-down', isCollapsed);
        });
    });
    // Label toggles expand/collapse
    const labels = container.querySelectorAll('.tree-label');
    labels.forEach(l => {
        l.addEventListener('click', function() {
            const node = this.closest('.tree-node');
            const toggle = node.querySelector('.tree-toggle');
            if (toggle) toggle.click();
        });
    });
}

async function refreshSelectedMetabolitesList() {
    const target = document.getElementById('metabolite-selection-list');
    if (!target) return;
    const selectedIds = Array.from(state.selectedMetaboliteCategories.get ? state.selectedMetaboliteCategories.get() : []);
    if (selectedIds.length === 0) {
        target.innerHTML = '<div class="text-muted p-3">Select categories to load metabolites…</div>';
        return;
    }
    target.innerHTML = '<div class="text-center p-3"><i class="fas fa-spinner fa-spin"></i> Loading metabolites…</div>';
    try {
        const all = [];
        for (const cid of selectedIds) {
            const items = await api.loadMetabolitesByCategory(cid);
            if (Array.isArray(items)) all.push(...items);
        }
        // De-duplicate by id
        const byId = new Map();
        all.forEach(m => { if (m && m.id != null) byId.set(String(m.id), m); });
        const items = Array.from(byId.values());
        if (items.length === 0) {
            target.innerHTML = '<div class="text-muted p-3">No metabolites for selected categories.</div>';
            return;
        }
        const html = items.map(m => {
            const id = `metab-${m.id}`;
            return `<div class=\"form-check mb-1\">\n                <input class=\"form-check-input metabolite-checkbox\" type=\"checkbox\" data-id=\"${m.id}\" id=\"${id}\">\n                <label class=\"form-check-label small\" for=\"${id}\">${m.name || ('Metabolite ' + m.id)}</label>\n            </div>`;
        }).join('');
        target.innerHTML = html;
        const boxes = target.querySelectorAll('.metabolite-checkbox');
        boxes.forEach(cb => {
            cb.addEventListener('change', function() {
                const mid = parseInt(this.dataset.id);
                if (this.checked) state.selectedMetabolites.add(mid); else state.selectedMetabolites.delete(mid);
            });
        });
    } catch (e) {
        console.error('refreshSelectedMetabolitesList failed:', e);
        target.innerHTML = '<div class="text-danger p-3">Failed to load metabolites.</div>';
    }
}

