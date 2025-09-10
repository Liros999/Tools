/**
 * Interactive Volcano Plot with Gene & Group Selection
 * Phase 1 MVP Implementation
 */

// Import plotManager
import * as plot from './modules/plotManager.js';

// Global state
let plotData = [];
let selectedGenes = new Set();
let selectedGroups = new Set();
let availableGroups = [];
let currentColorMap = new Map(); // gene -> color mapping
let labelsVisible = false;
let thresholdLabelsVisible = false;
let legendVisible = true;
let currentAnnotations = []; // Store gene label annotations
let groupHierarchyPaths = new Map(); // group -> path markings

// Expanded color palette for many groups
const SELECTION_COLORS = [
  '#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6',
  '#1abc9c', '#e67e22', '#34495e', '#f1c40f', '#e91e63',
  '#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4', '#ffeaa7',
  '#dda0dd', '#98d8c8', '#f7dc6f', '#bb8fce', '#85c1e9',
  '#f8c471', '#82e0aa', '#f1948a', '#85929e', '#d5a6bd',
  '#a9cce3', '#f9e79f', '#d7bde2', '#a3e4d7', '#fadbd8',
  '#aed6f1', '#fff2cc', '#e8daef', '#d1f2eb', '#fdebd0',
  '#ccd1d1', '#fcf3cf', '#ebdef0', '#d0ece7', '#fdf2e9'
];

class VolcanoPlotApp {
  constructor() {
    this.plotInstances = [];
    this.init();
  }

  async init() {
    console.log('🚀 Initializing Volcano Plot App');
    
    // Wait for DOM to be ready
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', () => this.start());
    } else {
      this.start();
    }
  }

  async start() {
    try {
      // Load available groups
      await this.loadGroups();
      
      // Setup UI event handlers
      this.setupEventHandlers();
      
      // Setup gene search autocomplete
      this.setupGeneSearch();
      
      // Render plots if data is available
      if (window.plotData && Array.isArray(window.plotData)) {
        plotData = window.plotData;
        // Use the plotManager instead of local renderPlots
        if (typeof plot !== 'undefined' && plot.initializePlots) {
          plot.initializePlots();
        } else {
          console.warn('PlotManager not available, using fallback renderPlots');
          this.renderPlots();
        }
      }
      
      console.log('✅ App initialized successfully');
    } catch (error) {
      console.error('❌ Failed to initialize app:', error);
    }
  }

  async loadGroups() {
    try {
      // Load flat groups
      const response = await fetch('/api/groups/main');
      const data = await response.json();
      
      if (data.success) {
        availableGroups = data.groups;
        // Note: Flat view removed, only tree view used
      } else {
        throw new Error(data.error || 'Failed to load groups');
      }
      
      // Load hierarchical tree data
      await this.loadHierarchicalGroups();
      
    } catch (error) {
      console.error('Error loading groups:', error);
      document.getElementById('group-selection-tree').innerHTML = 
        '<div class="text-danger text-center p-3">Error loading groups</div>';
    }
  }

  async loadHierarchicalGroups() {
    try {
      const response = await fetch('/api/groups/hierarchical');
      const data = await response.json();
      
      if (data.success) {
        this.cachedTreeData = data.tree; // Store for later use
        this.renderHierarchicalTree(data.tree);
      } else {
        throw new Error(data.error || 'Failed to load hierarchical groups');
      }
    } catch (error) {
      console.error('Error loading hierarchical groups:', error);
      document.getElementById('group-selection-tree').innerHTML = 
        '<div class="text-danger text-center p-3">Error loading tree</div>';
    }
  }

  // renderGroupSelection method removed - only tree view now

  renderHierarchicalTree(tree) {
    const container = document.getElementById('group-selection-tree');
    
    if (!tree || Object.keys(tree).length === 0) {
      container.innerHTML = '<div class="text-muted text-center p-3">No hierarchical data available</div>';
      return;
    }

    let html = this.buildTreeHTML(tree, 0);
    container.innerHTML = html;
  }

  buildTreeHTML(node, level) {
    let html = '';
    
    if (typeof node === 'object' && node.id !== undefined) {
      // This is a category node
      const hasChildren = node.children && Object.keys(node.children).length > 0;
      const indent = '  '.repeat(level);
      const geneDisplay = node.gene_count > 0 ? ` (${node.gene_count} genes)` : '';
      
      html += `
        <div class="tree-node" style="margin-left: ${level * 15}px;">
          <div class="d-flex align-items-center mb-1">
            ${hasChildren ? `
              <button class="btn btn-sm btn-link p-0 me-1 tree-toggle" 
                      data-target="tree-children-${node.id}" 
                      style="font-size: 12px; text-decoration: none;">
                <i class="fas fa-chevron-right"></i>
              </button>
            ` : '<span style="width: 20px; display: inline-block;"></span>'}
            
            <div class="form-check">
              <input class="form-check-input group-checkbox" type="checkbox" 
                     value="${node.id}" id="group-tree-${node.id}">
              <label class="form-check-label" for="group-tree-${node.id}">
                <small><strong>${node.name}</strong>${geneDisplay}</small><br>
                <small class="text-muted">${node.dot_notation} • Level ${node.level}</small>
              </label>
            </div>
          </div>
          
          ${hasChildren ? `
            <div class="tree-children" id="tree-children-${node.id}" style="display: none;">
              ${this.buildTreeHTML(node.children, level + 1)}
            </div>
          ` : ''}
        </div>
      `;
    } else {
      // This is a children object, iterate through its properties
      for (const [key, child] of Object.entries(node)) {
        html += this.buildTreeHTML(child, level);
      }
    }
    
    return html;
  }

  setupEventHandlers() {
    // Apply selection button
    document.getElementById('apply-selection').addEventListener('click', () => {
      this.applySelection();
    });

    // Clear selection button
    document.getElementById('clear-selection').addEventListener('click', () => {
      this.clearSelection();
    });

    // Removed flat view toggle buttons

    // Removed flat view event listeners

    // Group checkboxes for tree view (event delegation)
    document.getElementById('group-selection-tree').addEventListener('change', (e) => {
      if (e.target.classList.contains('group-checkbox')) {
        this.handleGroupSelection(e);
      }
    });

            // Tree expand/collapse buttons (event delegation)
        document.getElementById('group-selection-tree').addEventListener('click', (e) => {
          if (e.target.classList.contains('tree-toggle') || e.target.closest('.tree-toggle')) {
            const button = e.target.closest('.tree-toggle') || e.target;
            const targetId = button.getAttribute('data-target');
            const targetElement = document.getElementById(targetId);
            const icon = button.querySelector('i');
            const nodeContainer = button.closest('.tree-node');
            
            if (targetElement.style.display === 'none') {
              targetElement.style.display = 'block';
              icon.className = 'fas fa-chevron-down';
              // Add visual indicator for opened group
              nodeContainer.style.backgroundColor = '#f8f9fa';
              nodeContainer.style.border = '2px solid #007bff';
              nodeContainer.style.borderRadius = '8px';
              nodeContainer.style.padding = '4px';
            } else {
              targetElement.style.display = 'none';
              icon.className = 'fas fa-chevron-right';
              // Remove visual indicator
              nodeContainer.style.backgroundColor = '';
              nodeContainer.style.border = '';
              nodeContainer.style.borderRadius = '';
              nodeContainer.style.padding = '';
            }
          }
        });

    // Label toggle button
    document.getElementById('toggle-labels').addEventListener('click', () => {
      this.toggleGeneLabels();
    });

    // Threshold label toggle button
    document.getElementById('toggle-threshold-labels').addEventListener('click', () => {
      this.toggleThresholdLabels();
    });

    // Legend toggle button
    document.getElementById('toggle-legend').addEventListener('click', () => {
      this.toggleLegend();
    });

    // Export genes button
    document.getElementById('export-genes').addEventListener('click', () => {
      this.exportSelectedGenes();
    });

    // Synonym checking is now automatic in applySelection()
  }

  // switchView method removed - only tree view now

  handleGroupSelection(e) {
    const groupId = parseInt(e.target.value);
    if (e.target.checked) {
      selectedGroups.add(groupId);
      // Track hierarchy path for this group
      this.trackGroupHierarchy(groupId);
    } else {
      selectedGroups.delete(groupId);
      // Remove hierarchy path for this group
      this.clearGroupHierarchy(groupId);
    }
    this.updateSelectionSummary();
    this.updateTreePathMarkings();
  }

  setupGeneSearch() {
    const searchInput = document.getElementById('gene-search');
    let searchTimeout;

    searchInput.addEventListener('input', (e) => {
      clearTimeout(searchTimeout);
      const query = e.target.value.trim();

      if (query.length < 2) {
        return;
      }

      // Debounce search
      searchTimeout = setTimeout(() => {
        this.searchGenes(query);
      }, 300);
    });

    searchInput.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') {
        e.preventDefault();
        const query = e.target.value.trim();
        if (query.length >= 2) {
          this.addGeneToSelection(query);
          e.target.value = '';
        }
      }
    });
  }

  async searchGenes(query) {
    try {
      const response = await fetch(`/api/genes/search?q=${encodeURIComponent(query)}&limit=10`);
      const data = await response.json();
      
      if (data.success && data.genes.length > 0) {
        this.showGeneSearchResults(data.genes);
      }
    } catch (error) {
      console.error('Gene search error:', error);
    }
  }

  showGeneSearchResults(genes) {
    // Create dropdown for search results
    let dropdown = document.getElementById('gene-search-dropdown');
    if (!dropdown) {
      dropdown = document.createElement('div');
      dropdown.id = 'gene-search-dropdown';
      dropdown.className = 'list-group position-absolute w-100';
      dropdown.style.zIndex = '1000';
      dropdown.style.maxHeight = '200px';
      dropdown.style.overflowY = 'auto';
      
      const searchInput = document.getElementById('gene-search');
      searchInput.parentNode.style.position = 'relative';
      searchInput.parentNode.appendChild(dropdown);
    }

    let html = '';
    genes.forEach(gene => {
      html += `
        <button type="button" class="list-group-item list-group-item-action" 
                onclick="app.addGeneToSelection('${gene}')">
          ${gene}
        </button>
      `;
    });

    dropdown.innerHTML = html;

    // Hide dropdown when clicking outside
    setTimeout(() => {
      document.addEventListener('click', function hideDropdown(e) {
        if (!dropdown.contains(e.target) && e.target.id !== 'gene-search') {
          dropdown.innerHTML = '';
          document.removeEventListener('click', hideDropdown);
        }
      });
    }, 100);
  }

  addGeneToSelection(geneName) {
    if (!selectedGenes.has(geneName)) {
      selectedGenes.add(geneName);
      this.renderSelectedGenes();
      this.updateSelectionSummary();
      
      // Clear search
      document.getElementById('gene-search').value = '';
      const dropdown = document.getElementById('gene-search-dropdown');
      if (dropdown) dropdown.innerHTML = '';
    }
  }

  removeGeneFromSelection(geneName) {
    selectedGenes.delete(geneName);
    this.renderSelectedGenes();
    this.updateSelectionSummary();
  }

  renderSelectedGenes() {
    const container = document.getElementById('selected-genes');
    let html = '';

    selectedGenes.forEach(gene => {
      html += `
        <span class="badge bg-secondary me-1 mb-1">
          ${gene}
          <button type="button" class="btn-close btn-close-white ms-1" 
                  style="font-size: 0.7em;" 
                  onclick="app.removeGeneFromSelection('${gene}')"
                  aria-label="Remove ${gene}"></button>
        </span>
      `;
    });

    container.innerHTML = html;
  }

  async applySelection() {
    const selectionData = {
      groups: Array.from(selectedGroups),
      genes: Array.from(selectedGenes)
    };

    if (selectionData.groups.length === 0 && selectionData.genes.length === 0) {
      // Reset to grey if no selection
      this.resetPlotColors();
      return;
    }

    try {
      // Auto-run synonym check first if groups are selected
      let synonymAnalysis = null;
      if (selectedGroups.size > 0) {
        synonymAnalysis = await this.runSynonymCheck(selectionData);
      }

      const response = await fetch('/api/selection/genes', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(selectionData)
      });

      const data = await response.json();
      
      if (data.success) {
        this.colorSelectedGenesWithGroups(data.result, Array.from(selectedGroups));
        this.updateSelectionSummary(data.result.genes.length);
        
        // Show synonym results if any
        if (synonymAnalysis && synonymAnalysis.found_via_synonyms > 0) {
          this.displaySynonymResults(synonymAnalysis);
        }
      } else {
        throw new Error(data.error);
      }
    } catch (error) {
      console.error('Error applying selection:', error);
      alert('Error applying selection. Please try again.');
    }
  }

  colorSelectedGenesWithGroups(result, selectedGroupIds) {
    currentColorMap.clear();
    
    // Create group-to-color and group-to-name mapping
    const groupColorMap = new Map();
    const groupNameMap = new Map();
    const groupGeneMap = new Map(); // Track which genes belong to which groups
    
    // Initialize group mappings
    selectedGroupIds.forEach((groupId, index) => {
      groupColorMap.set(groupId, SELECTION_COLORS[index % SELECTION_COLORS.length]);
      groupGeneMap.set(groupId, []);
    });
    
    // Get group names from the tree data immediately  
    selectedGroupIds.forEach(groupId => {
      const groupName = this.getGroupNameFromTree(groupId);
      if (groupName) {
        groupNameMap.set(groupId, groupName);
      } else {
        groupNameMap.set(groupId, `Group ${groupId}`); // Fallback
      }
    });
    
    // Color genes based on their source groups and track group membership
    result.genes.forEach(gene => {
      const sourcesInfo = result.sources[gene] || [];
      let assignedColor = SELECTION_COLORS[9]; // Default for individual genes
      let assignedGroupId = null;
      
      // Find the first group this gene belongs to
      if (sourcesInfo.length > 0 && typeof sourcesInfo[0] === 'string' && sourcesInfo[0] !== 'Individual Selection') {
        // Handle string-based source (group name) - need to map back to group ID
        const groupName = sourcesInfo[0];
        for (const groupId of selectedGroupIds) {
          const expectedGroupName = groupNameMap.get(groupId);
          if (expectedGroupName === groupName) {
            assignedColor = groupColorMap.get(groupId);
            assignedGroupId = groupId;
            break;
          }
        }
      }
      
      // If no group found, assign to first available group (fallback)
      if (!assignedGroupId && selectedGroupIds.length > 0) {
        assignedGroupId = selectedGroupIds[0];
        assignedColor = groupColorMap.get(assignedGroupId);
      }
      
      currentColorMap.set(gene, assignedColor);
      
      // Track which genes belong to which group
      if (assignedGroupId && groupGeneMap.has(assignedGroupId)) {
        groupGeneMap.get(assignedGroupId).push(gene);
      }
    });

    // Store group information for legend creation
    this.currentGroupData = {
      groupColorMap,
      groupNameMap,
      groupGeneMap,
      individualGenes: Array.from(selectedGenes)
    };

    // Update all plots with new colors and proper legend
    this.updatePlotColorsWithDynamicLegend();
  }

  updatePlotColors() {
    this.plotInstances.forEach((plotDiv, index) => {
      if (plotData[index]) {
        this.recolorPlot(plotDiv, plotData[index]);
      }
    });
  }

  updatePlotColorsWithDynamicLegend() {
    this.plotInstances.forEach((plotDiv, index) => {
      if (plotData[index]) {
        const plot = plotData[index];
        const traces = [];
        
        // Background trace for uncolored genes
        const uncoloredGenes = { x: [], y: [], text: [] };
        
        for (let i = 0; i < plot.x.length; i++) {
          const geneName = plot.text[i];
          if (!currentColorMap.has(geneName)) {
            uncoloredGenes.x.push(plot.x[i]);
            uncoloredGenes.y.push(plot.y[i]);
            uncoloredGenes.text.push(plot.text[i]);
          }
        }
        
        // Add background trace
        if (uncoloredGenes.x.length > 0) {
          traces.push({
            x: uncoloredGenes.x,
            y: uncoloredGenes.y,
            text: uncoloredGenes.text,
            mode: 'markers',
            type: 'scattergl',
            marker: {
              size: 6,
              color: '#c0c0c0',
              line: { width: 0.5, color: '#333' }
            },
            hoverinfo: 'text',
            name: 'Unselected',
            showlegend: false
          });
        }
        
        // Create separate traces for each group (for proper legend)
        if (this.currentGroupData) {
          this.currentGroupData.groupColorMap.forEach((color, groupId) => {
            const groupGenes = this.currentGroupData.groupGeneMap.get(groupId) || [];
            if (groupGenes.length > 0) {
              const groupTrace = { x: [], y: [], text: [] };
              
              // Find genes for this group in the plot
              for (let i = 0; i < plot.x.length; i++) {
                const geneName = plot.text[i];
                if (groupGenes.includes(geneName)) {
                  groupTrace.x.push(plot.x[i]);
                  groupTrace.y.push(plot.y[i]);
                  groupTrace.text.push(plot.text[i]);
                }
              }
              
              if (groupTrace.x.length > 0) {
                const groupName = this.currentGroupData.groupNameMap.get(groupId) || `Group ${groupId}`;
                traces.push({
                  x: groupTrace.x,
                  y: groupTrace.y,
                  text: groupTrace.text,
                  mode: 'markers',
                  type: 'scattergl',
                  marker: {
                    size: 8,
                    color: color,
                    line: { width: 1, color: '#000' }
                  },
                  hoverinfo: 'text',
                  name: `${groupName} (${groupTrace.x.length})`,
                  showlegend: legendVisible
                });
              }
            }
          });
        }
        
        // Add individual genes trace if any
        if (this.currentGroupData && this.currentGroupData.individualGenes.length > 0) {
          const individualTrace = { x: [], y: [], text: [] };
          
          for (let i = 0; i < plot.x.length; i++) {
            const geneName = plot.text[i];
            if (this.currentGroupData.individualGenes.includes(geneName)) {
              individualTrace.x.push(plot.x[i]);
              individualTrace.y.push(plot.y[i]);
              individualTrace.text.push(plot.text[i]);
            }
          }
          
          if (individualTrace.x.length > 0) {
            traces.push({
              x: individualTrace.x,
              y: individualTrace.y,
              text: individualTrace.text,
              mode: 'markers',
              type: 'scattergl',
              marker: {
                size: 8,
                color: SELECTION_COLORS[9],
                line: { width: 1, color: '#000' }
              },
              hoverinfo: 'text',
              name: `Individual (${individualTrace.x.length})`,
              showlegend: legendVisible
            });
          }
        }
        
            // Preserve existing layout (including title) while updating traces
            const updatedLayout = {
                ...plotDiv.layout,  // Preserve existing layout properties
                // Only override specific properties we want to change
            };
            
            Plotly.react(plotDiv, traces, updatedLayout, {
                displayModeBar: true,
                modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d'],
                displaylogo: false
            }).then(() => {
                // CRITICAL: Re-setup tooltip handlers after plot update to maintain validation logic
                console.log(`🔧 Re-setting up tooltip handlers for ${plotDiv.id} after updatePlotColorsWithDynamicLegend`);
                // Import and call setupTooltipHandlers from plotManager
                if (window.plotManager && window.plotManager.setupTooltipHandlers) {
                    window.plotManager.setupTooltipHandlers();
                }
            });
      }
    });
  }

  updatePlotColorsWithLayering() {
    // Fallback method for backward compatibility
    this.updatePlotColorsWithDynamicLegend();
  }

  getGroupNameFromTree(targetId) {
    if (!this.cachedTreeData) return null;
    
    const found = this.findGroupInTree(this.cachedTreeData, targetId);
    return found ? found.name : null;
  }

  findGroupInTree(node, targetId) {
    if (typeof node === 'object' && node.id !== undefined) {
      // This is a category node
      if (node.id === targetId) {
        return node;
      }
      
      if (node.children) {
        const result = this.findGroupInTree(node.children, targetId);
        if (result) return result;
      }
    } else if (typeof node === 'object') {
      // This is a children object, iterate through its properties
      for (const [key, child] of Object.entries(node)) {
        const result = this.findGroupInTree(child, targetId);
        if (result) return result;
      }
    }
    
    return null;
  }

  trackGroupHierarchy(groupId) {
    // Find the path from root to this group
    const path = this.findPathToGroup(this.cachedTreeData, groupId, []);
    if (path.length > 0) {
      // Get color for this group (based on selection order)
      const groupIndex = Array.from(selectedGroups).indexOf(groupId);
      const groupColor = SELECTION_COLORS[groupIndex % SELECTION_COLORS.length];
      groupHierarchyPaths.set(groupId, { path, color: groupColor });
    }
  }

  clearGroupHierarchy(groupId) {
    groupHierarchyPaths.delete(groupId);
  }

  findPathToGroup(node, targetId, currentPath) {
    if (typeof node === 'object' && node.id !== undefined) {
      // This is a category node
      const newPath = [...currentPath, node.id];
      
      if (node.id === targetId) {
        return newPath;
      }
      
      if (node.children) {
        const result = this.findPathToGroup(node.children, targetId, newPath);
        if (result.length > 0) return result;
      }
    } else if (typeof node === 'object') {
      // This is a children object, iterate through its properties
      for (const [key, child] of Object.entries(node)) {
        const result = this.findPathToGroup(child, targetId, currentPath);
        if (result.length > 0) return result;
      }
    }
    
    return [];
  }

  updateTreePathMarkings() {
    // Clear all existing path markings
    document.querySelectorAll('.tree-node').forEach(node => {
      const checkbox = node.querySelector('.group-checkbox');
      if (checkbox) {
        const nodeId = parseInt(checkbox.value);
        
        // Reset styling
        node.style.borderLeft = '';
        node.style.paddingLeft = '';
        
        // Apply path markings for selected groups
        groupHierarchyPaths.forEach((pathData, groupId) => {
          if (pathData.path.includes(nodeId)) {
            // Add V-mark indicator for path
            node.style.borderLeft = `4px solid ${pathData.color}`;
            node.style.paddingLeft = '8px';
            
            // Add visual indicator icon if this is part of the path
            const existingIcon = node.querySelector('.path-indicator');
            if (!existingIcon) {
              const pathIcon = document.createElement('span');
              pathIcon.className = 'path-indicator';
              pathIcon.innerHTML = ' <i class="fas fa-chevron-right" style="font-size: 8px;"></i>';
              pathIcon.style.color = pathData.color;
              
              const label = node.querySelector('.form-check-label');
              if (label) {
                label.appendChild(pathIcon);
              }
            }
          }
        });
      }
    });
  }

  recolorPlot(plotDiv, plot) {
    const colors = plot.x.map((_, i) => {
      const geneName = plot.text[i]; // Assuming text contains gene names
      return currentColorMap.get(geneName) || 'grey';
    });

    const update = {
      'marker.color': [colors]
    };

    Plotly.restyle(plotDiv, update, 0);
  }

  resetPlotColors() {
    currentColorMap.clear();
    this.plotInstances.forEach((plotDiv, index) => {
      if (plotData[index]) {
        const colors = new Array(plotData[index].x.length).fill('#c0c0c0'); // Brighter grey
        const update = { 'marker.color': [colors] };
        Plotly.restyle(plotDiv, update, 0);
      }
    });
    this.updateSelectionSummary(0);
  }

  clearSelection() {
    selectedGenes.clear();
    selectedGroups.clear();
    groupHierarchyPaths.clear(); // Clear hierarchy paths
    
    // Uncheck all group checkboxes in tree view
    document.querySelectorAll('.group-checkbox').forEach(cb => {
      cb.checked = false;
    });
    
    // Close all opened tree nodes and remove visual indicators
    document.querySelectorAll('.tree-node').forEach(node => {
      node.style.backgroundColor = '';
      node.style.border = '';
      node.style.borderRadius = '';
      node.style.padding = '';
      node.style.borderLeft = ''; // Clear path markings
      node.style.paddingLeft = '';
      
      // Remove path indicator icons
      const pathIndicator = node.querySelector('.path-indicator');
      if (pathIndicator) {
        pathIndicator.remove();
      }
    });
    
    document.querySelectorAll('.tree-children').forEach(children => {
      children.style.display = 'none';
    });
    
    document.querySelectorAll('.tree-toggle i').forEach(icon => {
      icon.className = 'fas fa-chevron-right';
    });
    
    this.renderSelectedGenes();
    this.resetPlotColors();
    
    // Hide synonym results
    document.getElementById('synonym-results').style.display = 'none';
  }

  updateSelectionSummary(totalGenes = null) {
    const container = document.getElementById('selection-summary');
    const groupCount = selectedGroups.size;
    const geneCount = selectedGenes.size;
    
    if (totalGenes !== null) {
      container.textContent = `${totalGenes} genes selected from ${groupCount} groups + ${geneCount} individual genes`;
    } else if (groupCount > 0 || geneCount > 0) {
      container.textContent = `${groupCount} groups + ${geneCount} genes selected (click Apply to see results)`;
    } else {
      container.textContent = 'No genes selected';
    }
  }

    renderPlots() {
        const plotContainer = document.getElementById('plot-container');
    plotContainer.innerHTML = '';

    plotData.forEach((plot, index) => {
            const plotDiv = document.createElement('div');
      plotDiv.className = 'mb-3';
      plotDiv.id = 'plot-' + index;
      plotDiv.style.height = plotData.length === 1 ? '600px' : '500px';
            plotContainer.appendChild(plotDiv);

            const trace = {
                x: plot.x,
                y: plot.y,
                text: plot.text,
                mode: 'markers',
        type: 'scattergl',
                marker: { 
          size: 6,
          color: '#c0c0c0', // Brighter grey
          line: { width: 0.5, color: '#333' }
                },
        hoverinfo: 'skip',
                    name: 'Non-selected Genes'
            };

            const layout = {
        title: {
          text: plot.title || ('Plot ' + (index + 1)),
          font: { size: 16, family: 'Arial, sans-serif' }
        },
        xaxis: { 
          title: { text: 'Log2 Fold Change', font: { size: 14 } },
          zeroline: true,
          zerolinecolor: '#ddd'
        },
        yaxis: { 
          title: { text: '-log10(p-value)', font: { size: 14 } },
          zeroline: true,
          zerolinecolor: '#ddd'
        },
        hovermode: 'closest', // Changed from 'x unified' to 'closest' for more precise tooltip behavior
        shapes: this.createThresholdLines(plot),
        showlegend: legendVisible,
        legend: {
          x: 1.02,
          y: 1,
          bgcolor: 'rgba(255,255,255,0.8)',
          bordercolor: '#ddd',
          borderwidth: 1
        },
        annotations: [], // Will be populated with gene labels
        margin: { l: 60, r: 120, t: 60, b: 60 }
      };

      Plotly.newPlot(plotDiv, [trace], layout, { 
        displayModeBar: true,
        modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian'],
        modeBarButtonsToAdd: [
          {
            name: 'Reset to original view',
            icon: Plotly.Icons.home,
            click: function(gd) {
              Plotly.relayout(gd, {
                'xaxis.autorange': true,
                'yaxis.autorange': true
              });
            }
          },
          {
            name: 'Download plot as PNG',
            icon: Plotly.Icons.camera,
            click: function(gd) {
              Plotly.downloadImage(gd, {
                format: 'png',
                width: 1200,
                height: 800,
                filename: `volcano_plot_${index + 1}`
              });
            }
          }
        ],
        displaylogo: false,
        responsive: true
      });
      
      // Setup hover events for gene info
      this.setupPlotHover(plotDiv);
      
      this.plotInstances.push(plotDiv);
    });
  }

  createThresholdLines(plot) {
    const shapes = [];
    
    if (plot.p_thresh_log != null) {
      shapes.push({
        type: 'line', xref: 'paper', x0: 0, x1: 1,
                        y0: plot.p_thresh_log, y1: plot.p_thresh_log,
                        line: { color: 'grey', width: 1, dash: 'dot' },
                        hoverinfo: 'skip',
                        hoveron: 'points', // Only hover on actual data points, not on the line
                        hovermode: 'closest' // Use closest mode for this specific shape
      });
    }
    
    if (plot.fc_thresh != null) {
      shapes.push({
        type: 'line', yref: 'paper', y0: 0, y1: 1,
                        x0: plot.fc_thresh, x1: plot.fc_thresh,
                        line: { color: 'grey', width: 1, dash: 'dot' },
                        hoverinfo: 'skip',
                        hoveron: 'points', // Only hover on actual data points, not on the line
                        hovermode: 'closest' // Use closest mode for this specific shape
      });
      shapes.push({
        type: 'line', yref: 'paper', y0: 0, y1: 1,
                        x0: -plot.fc_thresh, x1: -plot.fc_thresh,
                        line: { color: 'grey', width: 1, dash: 'dot' },
                        hoverinfo: 'skip',
                        hoveron: 'points', // Only hover on actual data points, not on the line
                        hovermode: 'closest' // Use closest mode for this specific shape
      });
    }
    
    return shapes;
  }

  setupPlotHover(plotDiv) {
    // REMOVED: Duplicate tooltip implementation - using plotManager.js tooltips instead
    
    // Only handle click events for gene info modal
    plotDiv.on('plotly_click', async (data) => {
      // VALIDATION: Only process clicks on valid data points
      if (!data || !data.points || !data.points[0] || !data.points[0].text) {
        return; // Don't process clicks on invalid data
      }
      
      const point = data.points[0];
      const geneName = point.text;
      
      // Additional validation for coordinates
      if (typeof point.x !== 'number' || typeof point.y !== 'number') {
        return; // Don't process clicks with invalid coordinates
      }
      
      if (geneName && geneName.trim() !== '') {
        await this.showGeneInfo(geneName);
      }
    });
  }

  async showGeneInfo(geneName) {
    try {
      const response = await fetch(`/api/gene/${geneName}/info`);
      const data = await response.json();
      
      if (data.success) {
        this.displayGeneInfoModal(data.gene_info);
      } else {
        throw new Error(data.error);
      }
    } catch (error) {
      console.error('Error fetching gene info:', error);
    }
  }

  displayGeneInfoModal(geneInfo) {
    const modalBody = document.getElementById('gene-info-content');
    const modalTitle = document.getElementById('geneInfoModalLabel');
    
    modalTitle.textContent = `Gene: ${geneInfo.gene_name}`;
    
    let html = `
      <p><strong>Gene ID:</strong> ${geneInfo.gene_id}</p>
      <p><strong>Number of Groups:</strong> ${geneInfo.group_count}</p>
      <h6>Belongs to Groups:</h6>
      <ul class="list-group">
    `;
    
    geneInfo.groups.forEach(group => {
      html += `
        <li class="list-group-item d-flex justify-content-between align-items-center">
          ${group.name}
          <span class="badge bg-primary rounded-pill">Level ${group.level}</span>
        </li>
      `;
    });
    
    html += '</ul>';
    modalBody.innerHTML = html;
    
    // Show modal
    const modal = new bootstrap.Modal(document.getElementById('geneInfoModal'));
    modal.show();
  }

  toggleGeneLabels() {
    labelsVisible = !labelsVisible;
    const button = document.getElementById('toggle-labels');
    
    if (labelsVisible) {
      button.classList.remove('btn-outline-primary');
      button.classList.add('btn-primary');
      // Use the correct labeling function from plotManager that only shows selected genes
      if (typeof plot !== 'undefined' && plot.addGeneLabelsForSelected) {
        plot.addGeneLabelsForSelected();
      } else {
        console.warn('PlotManager not available, fallback not implemented - labels will not be shown');
      }
    } else {
      button.classList.remove('btn-primary');
      button.classList.add('btn-outline-primary');
      this.removeGeneLabels();
    }
  }

  toggleThresholdLabels() {
    thresholdLabelsVisible = !thresholdLabelsVisible;
    const button = document.getElementById('toggle-threshold-labels');
    
    if (thresholdLabelsVisible) {
      button.classList.remove('btn-outline-info');
      button.classList.add('btn-info');
      this.addThresholdLabels();
    } else {
      button.classList.remove('btn-info');
      button.classList.add('btn-outline-info');
      this.removeThresholdLabels();
    }
  }



  removeGeneLabels() {
    // Use the correct label removal function from plotManager
    if (typeof plot !== 'undefined' && plot.removeGeneLabels) {
      plot.removeGeneLabels();
    } else {
      console.warn('PlotManager not available, using fallback removeGeneLabels');
      // Fallback: clear annotations from all plot instances
      this.plotInstances.forEach(plotDiv => {
        Plotly.relayout(plotDiv, { annotations: [] });
      });
    }
  }

  addThresholdLabels() {
    this.plotInstances.forEach((plotDiv, index) => {
      const plot = plotData[index];
      const annotations = [];
      
      // Get thresholds from plot data
      const fcThresh = plot.fc_thresh || 1; // Default fold change threshold
      const pThresh = plot.p_thresh_log || 1.3; // Default p-value threshold (-log10)
      
      // Add labels for genes above/below thresholds
      plot.x.forEach((x, i) => {
        const geneName = plot.text[i];
        const logFC = plot.x[i];
        const negLogP = plot.y[i];
        
        // Check if gene meets threshold criteria
        const meetsThreshold = (
          (Math.abs(logFC) >= fcThresh) && 
          (negLogP >= pThresh)
        );
        
        if (meetsThreshold) {
          annotations.push({
            x: logFC,
            y: negLogP,
            text: geneName,
            showarrow: true,
            arrowhead: 2,
            arrowsize: 1,
            arrowwidth: 1,
            arrowcolor: '#666',
            font: { size: 9, color: '#333' },
            bgcolor: 'rgba(255,255,255,0.8)',
            bordercolor: '#999',
            borderwidth: 1,
            borderradius: 2,
            captureevents: true,
            // Store original position for reset capability
            _originalX: logFC,
            _originalY: negLogP,
            _geneName: geneName,
            _isThresholdLabel: true
          });
        }
      });
      
      Plotly.relayout(plotDiv, { annotations: annotations });
      this.setupLabelDragging(plotDiv);
    });
  }

  removeThresholdLabels() {
    this.plotInstances.forEach(plotDiv => {
      // Remove only threshold labels, keep gene labels if visible
      if (labelsVisible) {
        // Use the correct labeling function from plotManager that only shows selected genes
        if (typeof plot !== 'undefined' && plot.addGeneLabelsForSelected) {
          plot.addGeneLabelsForSelected();
        } else {
          console.warn('PlotManager not available, fallback not implemented - labels will not be shown');
        }
      } else {
        Plotly.relayout(plotDiv, { annotations: [] });
      }
    });
  }

  setupLabelDragging(plotDiv) {
    // Enable draggable annotations
    plotDiv.on('plotly_relayout', (eventData) => {
      // This event fires when annotations are moved
      if (eventData.annotations) {
        // Store the updated annotation positions
        console.log('Labels moved:', eventData.annotations);
      }
    });
    
    // Add double-click to reset label position
    plotDiv.on('plotly_doubleclick', () => {
      // Reset all labels to original positions
      this.resetLabelPositions();
      return false; // Prevent default zoom reset
    });
  }

  resetLabelPositions() {
    this.plotInstances.forEach(plotDiv => {
      const currentLayout = plotDiv.layout;
      if (currentLayout.annotations) {
        const resetAnnotations = currentLayout.annotations.map(ann => {
          if (ann._originalX !== undefined && ann._originalY !== undefined) {
            return {
              ...ann,
              x: ann._originalX,
              y: ann._originalY
            };
          }
          return ann;
        });
        
        Plotly.relayout(plotDiv, { annotations: resetAnnotations });
      }
    });
  }

  toggleLegend() {
    legendVisible = !legendVisible;
    const button = document.getElementById('toggle-legend');
    
    if (legendVisible) {
      button.classList.remove('btn-outline-secondary');
      button.classList.add('btn-secondary');
    } else {
      button.classList.remove('btn-secondary');
      button.classList.add('btn-outline-secondary');
    }

    // Re-render plots with updated legend visibility
    this.updatePlotColorsWithDynamicLegend();
  }

  exportSelectedGenes() {
    if (selectedGroups.size === 0 && selectedGenes.size === 0) {
      alert('No genes selected to export');
      return;
    }

    // Get current selection
    const selectionData = {
      groups: Array.from(selectedGroups),
      genes: Array.from(selectedGenes)
    };

    fetch('/api/selection/genes', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(selectionData)
    })
                .then(response => response.json())
                .then(data => {
      if (data.success) {
        this.downloadGeneList(data.result);
      } else {
        alert('Error exporting genes: ' + data.error);
      }
    })
    .catch(error => {
      console.error('Export error:', error);
      alert('Error exporting genes');
    });
  }

  downloadGeneList(geneData) {
    const csvContent = this.generateCSV(geneData);
    const blob = new Blob([csvContent], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    
    const a = document.createElement('a');
    a.href = url;
    a.download = `selected_genes_${new Date().toISOString().slice(0,10)}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
  }

  generateCSV(geneData) {
    let csv = 'Gene Name,Source Groups\n';
    
    geneData.genes.forEach(gene => {
      const sources = geneData.sources[gene] || [];
      csv += `"${gene}","${sources.join('; ')}"\n`;
    });
    
    return csv;
  }

  async runSynonymCheck(selectionData) {
    // Automatically check for synonyms (no UI button needed)
    try {
      // In a real implementation, this would come from the uploaded file
      // For now, simulate some common gene names that might need synonym resolution
      const simulatedUploadedGenes = [
        'yvgX', 'copZ', 'dapA', 'dapF', 'artP', 'cueR', 'unknownGene1', 'unknownGene2'
      ];

      const response = await fetch('/api/synonyms/resolve', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          groups: selectionData.groups,
          uploaded_genes: simulatedUploadedGenes
        })
      });

      const data = await response.json();
      
      if (data.success) {
        return data.synonym_analysis;
      } else {
        console.warn('Synonym check failed:', data.error);
        return null;
      }
    } catch (error) {
      console.warn('Synonym check error:', error);
      return null;
    }
  }

  async checkSynonyms() {
    if (selectedGroups.size === 0) {
      alert('Please select at least one group to check for synonyms');
      return;
    }

    const selectionData = {
      groups: Array.from(selectedGroups),
      genes: Array.from(selectedGenes)
    };

    const button = document.getElementById('check-synonyms');
    const originalHTML = button.innerHTML;
    
    try {
      button.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Checking...';
      button.disabled = true;

      const synonymAnalysis = await this.runSynonymCheck(selectionData);
      
      if (synonymAnalysis) {
        this.displaySynonymResults(synonymAnalysis);
      } else {
        alert('No synonym data available');
      }
    } catch (error) {
      console.error('Synonym check error:', error);
      alert('Error checking synonyms: ' + error.message);
    } finally {
      button.innerHTML = originalHTML;
      button.disabled = false;
    }
  }

  displaySynonymResults(analysis) {
    const resultsDiv = document.getElementById('synonym-results');
    const summaryDiv = document.getElementById('synonym-summary');
    const mappingsDiv = document.getElementById('synonym-mappings');

    // Show summary
    summaryDiv.innerHTML = `
      <div class="row text-center">
        <div class="col-3">
          <div class="text-success"><strong>${analysis.found_direct}</strong></div>
          <div class="small">Direct</div>
        </div>
        <div class="col-3">
          <div class="text-info"><strong>${analysis.found_via_synonyms}</strong></div>
          <div class="small">Synonyms</div>
        </div>
        <div class="col-3">
          <div class="text-warning"><strong>${analysis.still_missing}</strong></div>
          <div class="small">Missing</div>
        </div>
        <div class="col-3">
          <div class="text-primary"><strong>${analysis.total_expected}</strong></div>
          <div class="small">Total</div>
        </div>
      </div>
    `;

    // Show synonym mappings if any
    if (Object.keys(analysis.synonym_mappings).length > 0) {
      let mappingsHTML = '<h6 class="small mb-2">Synonym Mappings Found:</h6>';
      for (const [synonym, canonical] of Object.entries(analysis.synonym_mappings)) {
        mappingsHTML += `
          <div class="small mb-1">
            <span class="badge bg-light text-dark">${synonym}</span> 
            <i class="fas fa-arrow-right mx-1"></i> 
            <span class="badge bg-info">${canonical}</span>
          </div>
        `;
      }
      mappingsDiv.innerHTML = mappingsHTML;
    } else {
      mappingsDiv.innerHTML = '<div class="small text-muted">No synonym mappings found</div>';
    }

    // Show missing genes if any
    if (analysis.missing_genes.length > 0) {
      mappingsDiv.innerHTML += `
        <h6 class="small mb-2 mt-3">Still Missing:</h6>
        <div class="small">
          ${analysis.missing_genes.slice(0, 10).map(gene => 
            `<span class="badge bg-light text-dark me-1">${gene}</span>`
          ).join('')}
          ${analysis.missing_genes.length > 10 ? `<span class="text-muted">... and ${analysis.missing_genes.length - 10} more</span>` : ''}
        </div>
      `;
    }

    resultsDiv.style.display = 'block';
  }
}

// Initialize app
const app = new VolcanoPlotApp();

// Make app globally available for onclick handlers
window.app = app;