/**
 * Plot Manager - Handles all Plotly.js interactions and plot rendering
 */

import { 
    plotData, plotInstances, currentColorMap, labelsVisible, 
    legendVisible, pinnedGene, isTooltipPinned, FILTER_TYPE_COLORS,
    INTERSECTION_COLOR, INTERSECTION_MARKER_SIZE,
    selectedGroups, selectedRegulons, selectedGenes, currentSelectionResult,
    coloringMode, selectedPathways, geneToPathwayMap, pathwayInfo,
    selectedComplexes, proteinComplexMap
} from './state.js';

import GeneOverlapDebugger from './debugManager.js';

// Initialize comprehensive debugger
const geneDebugger = new GeneOverlapDebugger({
    debugMode: true,
    logLevel: 'verbose',
    enableVisualDebugging: true
});

// ===========================================================================
// PLOT COLORING AND HIGHLIGHTING
// ===========================================================================

export function updatePlotWithSelection(result) {
    console.log('🎨 DYNAMIC LEGEND UPDATE: Executing comprehensive plot update with:', result);
    console.log('🎨 Categories in result:', Object.keys(result.categories || {}));
    console.log('🎨 Category names in result:', result.category_names || {});
    
    // STORE CURRENT SELECTION RESULT FOR LABELS
    currentSelectionResult.set(result);
    
    // DEFENSIVE CHECK: Ensure plots are initialized first
    const plotContainers = document.querySelectorAll('[id^="plot-"]:not(#plot-container):not(#plot-data):not(#plot-loading)');
    console.log('🔍 Found plot containers:', Array.from(plotContainers).map(c => c.id));
    
    if (!plotContainers.length) {
        console.warn('⚠️ No plot containers found. Initializing plots first...');
        initializePlots();
        return;
    }

    const categories = result.categories || result;  // Support both old and new format
    const categoryNames = result.category_names || {};  // Get display names

    // PERFORMANCE OPTIMIZATION: Batch DOM updates
    const startTime = performance.now();
    
    plotContainers.forEach(container => {
        const plotData = window[container.id + '_data'];
        if (!plotData) {
            console.error('❌ No plot data for', container.id);
            return;
        }

        // STEP 1: ANALYZE GENE INTERSECTIONS (OPTIMIZED)
        const geneIntersections = analyzeGeneIntersections(categories);
        console.log('🔍 Gene intersection analysis:', geneIntersections);
        console.log('🔍 Categories received:', Object.keys(categories));
        Object.entries(categories).forEach(([key, genes]) => {
            console.log(`🔍 Category ${key}: ${genes ? genes.length : 0} genes`);
        });

        // STEP 2: CREATE DYNAMIC TRACES WITH PROPER COLORING (OPTIMIZED)
        const traces = buildDynamicTraces(plotData, categories, geneIntersections, categoryNames);
        
        // STEP 3: UPDATE PLOT WITH NEW TRACES AND LEGEND (BATCH UPDATE)
        // Get the current plot element to access its layout
        const plotElement = document.getElementById(container.id);
        const currentLayout = plotElement ? plotElement.layout : {};
        
        Plotly.react(container.id, traces, {
            ...currentLayout,  // Preserve existing layout (including title, axes, etc.)
            showlegend: true,
            legend: {
                x: 1.02,
                y: 1,
                xanchor: 'left',
                yanchor: 'top',
                bgcolor: 'rgba(255,255,255,0.9)',
                bordercolor: '#333',
                borderwidth: 1
            }
        });
        
        console.log(`✅ Plot ${container.id} updated with ${traces.length} dynamic traces and legend.`);
    });
    
    const endTime = performance.now();
    console.log(`⚡ Plot update completed in ${(endTime - startTime).toFixed(2)}ms`);
}

function analyzeGeneIntersections(categories) {
    // Use comprehensive debugging framework for gene overlap detection
    geneDebugger.log("=== GENE INTERSECTION ANALYSIS ===", 'critical');
    
    // PERFORMANCE OPTIMIZATION: Use Maps for faster lookups
    const intersections = new Map();
    const geneToCategories = new Map();
    
    // OPTIMIZED: Single pass through categories
    for (const [categoryName, genes] of Object.entries(categories)) {
        if (genes && Array.isArray(genes) && genes.length > 0) {
            geneDebugger.log(`Processing category: ${categoryName} with ${genes.length} genes`, 'critical');
            for (const gene of genes) {
                if (!geneToCategories.has(gene)) {
                    geneToCategories.set(gene, []);
                }
                geneToCategories.get(gene).push(categoryName);
            }
        }
    }
    
    // OPTIMIZED: Find intersection genes - ONLY GENES IN ALL CATEGORIES
    const intersectionGroups = new Map();
    const allIntersectionGenes = [];
    const totalCategories = Object.keys(categories).length;
    
    geneDebugger.log('Analyzing gene-category relationships:', 'critical');
    geneDebugger.log(`Total selected categories: ${totalCategories}`, 'critical');
    
    for (const [gene, cats] of geneToCategories.entries()) {
        // Only include genes that appear in ALL selected categories
        if (cats.length === totalCategories && totalCategories > 1) {
            allIntersectionGenes.push(gene);
            geneDebugger.log(`Multi-category gene: ${gene} appears in ALL ${cats.length} categories: ${cats.join(', ')}`, 'critical');
        }
    }
    
    // Create single intersection group if there are any multi-category genes
    if (allIntersectionGenes.length > 0) {
        intersectionGroups.set('Multi-Category Genes', allIntersectionGenes);
    }
    
    // Convert Maps back to objects for compatibility
    const intersectionObj = {};
    const geneToCategoriesObj = {};
    
    intersectionGroups.forEach((genes, key) => {
        intersectionObj[key] = genes;
    });
    
    geneToCategories.forEach((cats, gene) => {
        geneToCategoriesObj[gene] = cats;
    });
    
    geneDebugger.log('Found intersection groups:', Object.keys(intersectionObj), 'critical');
    Object.entries(intersectionObj).forEach(([key, genes]) => {
        geneDebugger.log(`Intersection "${key}": ${genes.length} genes`, 'critical');
    });
    
    // Use debugging framework to validate results
    const debugResult = geneDebugger.detectOverlaps(
        Object.keys(categories), 
        Object.entries(categories).map(([name, genes]) => ({ id: name, name: name, genes: genes }))
    );
    
    return { intersections: intersectionObj, geneToCategories: geneToCategoriesObj };
}

function buildDynamicTraces(plotData, categories, intersectionData, categoryNames = {}) {
    // Dispatch based on coloring mode
    if (coloringMode.get && coloringMode.get() === 'pathway') {
        return buildPathwayTraces(plotData);
    }

    const traces = [];
    const colorPalette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'];
    let colorIndex = 0;

    // === NEW LOGIC: Correctly identify intersections and update group counts ===
    const geneToCategoriesMap = new Map();
    Object.entries(categories).forEach(([catKey, genes]) => {
        if (!genes) return;
        genes.forEach(gene => {
            if (!geneToCategoriesMap.has(gene)) {
                geneToCategoriesMap.set(gene, []);
            }
            geneToCategoriesMap.get(gene).push(catKey);
        });
    });

    const intersectionGenes = new Set();
    const totalSelectedCategories = Object.keys(categories).length;
    geneToCategoriesMap.forEach((cats, gene) => {
        // Only include genes that appear in ALL selected categories (true intersection)
        if (cats.length === totalSelectedCategories && totalSelectedCategories > 1) {
            intersectionGenes.add(gene);
        }
    });

    // Create a mutable copy of categories to subtract intersections
    const adjustedCategories = {};
    Object.entries(categories).forEach(([catKey, genes]) => {
        if (genes) {
            adjustedCategories[catKey] = genes.filter(gene => !intersectionGenes.has(gene));
        }
    });
    // === END NEW LOGIC ===

    // BASE TRACE: Non-selected genes (light grey)
    const baseTrace = {
        x: [], y: [], text: [],
        mode: 'markers', type: 'scatter',
        name: 'Non-selected Genes',
        marker: { color: '#cccccc', size: 6, opacity: 0.7 },
        hoverinfo: 'skip',
        showlegend: true
    };

    // CATEGORY TRACES: Use the ADJUSTED categories
    const categoryTraces = {};
    Object.entries(adjustedCategories).forEach(([categoryKey, genes]) => {
        if (genes && genes.length > 0) {
            const displayName = categoryNames[categoryKey] || categoryKey;
            categoryTraces[categoryKey] = {
                x: [], y: [], text: [],
                mode: 'markers', type: 'scatter',
                name: `${displayName} (${genes.length})`, // Correct, updated count
                marker: { 
                    color: colorPalette[colorIndex % colorPalette.length], 
                    size: 8, 
                    opacity: 0.9 
                },
                hovertemplate: '<b>%{text}</b><br>Log2FC: %{x}<br>-log10(p): %{y}<br>Category: ' + displayName + '<extra></extra>',
                showlegend: true
            };
            colorIndex++;
        }
    });

    // INTERSECTION TRACE: Genes in multiple categories (PURPLE WITH RED CIRCLE)
    const intersectionTrace = {
        x: [], y: [], text: [],
        mode: 'markers', type: 'scatter',
        name: `Intersection (${intersectionGenes.size})`, // Correct intersection count
        marker: {
            color: '#9467bd', // Neutral intersection color (purple)
            size: 10,
            symbol: 'circle',
            line: {
                color: 'red', // Red circle outline
                width: 2
            }
        },
        hovertemplate: '<b>%{text}</b><br>Log2FC: %{x}<br>-log10(p): %{y}<br>Intersection Gene<extra></extra>',
        showlegend: intersectionGenes.size > 0 // Only show if intersections exist
    };

    // ASSIGN GENES TO TRACES
    const allSelectedGenes = new Set(Object.keys(Object.fromEntries(geneToCategoriesMap)));

    plotData.text.forEach((geneName, index) => {
        const point = {
            x: plotData.x[index],
            y: plotData.y[index],
            text: geneName
        };

        if (intersectionGenes.has(geneName)) {
            intersectionTrace.x.push(point.x);
            intersectionTrace.y.push(point.y);
            intersectionTrace.text.push(point.text);
        } else if (allSelectedGenes.has(geneName)) {
            const categoryKey = geneToCategoriesMap.get(geneName)[0];
            if (categoryTraces[categoryKey]) {
                categoryTraces[categoryKey].x.push(point.x);
                categoryTraces[categoryKey].y.push(point.y);
                categoryTraces[categoryKey].text.push(point.text);
            }
        } else {
            baseTrace.x.push(point.x);
            baseTrace.y.push(point.y);
            baseTrace.text.push(point.text);
        }
    });

    // BUILD FINAL TRACE LIST
    traces.push(baseTrace);
    Object.values(categoryTraces).forEach(trace => {
        if (trace.x.length > 0) traces.push(trace);
    });
    if (intersectionTrace.x.length > 0) {
        traces.push(intersectionTrace);
    }
    
    console.log(`🎨 Built ${traces.length} dynamic traces with intersections highlighted`);
    return traces;
}

export function highlightIntersections(intersections) {
    const instances = plotInstances.get();
    if (!instances || instances.length === 0) return;
    
    const plots = plotData.get();
    if (!plots || plots.length === 0) return;
    
    // Create a set of intersection genes
    const intersectionGenes = new Set(intersections.map(item => item.gene));
    
    // Update each plot
    instances.forEach((plotDiv, plotIndex) => {
        const plot = plots[plotIndex];
        if (!plot || !plot.data) return;
        
        // Add red circles around intersection genes
        intersectionGenes.forEach(gene => {
            const geneIndex = plot.gene_names.indexOf(gene);
            if (geneIndex !== -1) {
                // Make the point larger and add red border
                Plotly.restyle(plotDiv, {
                    'marker.size': INTERSECTION_MARKER_SIZE,
                    'marker.line.color': INTERSECTION_COLOR,
                    'marker.line.width': 3
                }, [geneIndex]);
            }
        });
    });
}

function getCategoryType(categoryName) {
    if (categoryName.startsWith('group_')) return 'group';
    if (categoryName.startsWith('regulation_')) return 'regulation';
    if (categoryName.startsWith('regulon_')) return 'regulon';
    if (categoryName.startsWith('reaction_')) return 'reaction';
    if (categoryName.startsWith('material_')) return 'material';
    if (categoryName.startsWith('pathway_')) return 'pathway';
    if (categoryName.startsWith('individual')) return 'individual';
    return 'unknown';
}



export function resetPlotColors() {
    const instances = plotInstances.get();
    if (!instances || instances.length === 0) return;
    // Reset colors and autorange on each existing plot without recreating DOM
    instances.forEach((plotDiv, index) => {
        const plots = plotData.get();
        const pd = plots[index];
        // Build a single grey trace replacing current data
        const greyTrace = {
            x: pd.x,
            y: pd.y,
            text: pd.text,
            mode: 'markers',
            type: 'scatter',
            name: 'Non-selected Genes',
            marker: { color: '#cccccc', size: 6, opacity: 0.7 },
            hoverinfo: 'skip'
        };
        Plotly.react(plotDiv, [greyTrace], {
            title: pd.title || 'Volcano Plot',
            xaxis: { title: (pd.x_title || 'Log2 Fold Change'), zeroline: true, zerolinecolor: '#000000', zerolinewidth: 1, autorange: true },
            yaxis: { title: (pd.y_title || '-log10(p-value)'), zeroline: false, autorange: true },
            hovermode: 'closest', // Changed from 'x unified' to 'closest' for more precise tooltip behavior
            showlegend: true,
            margin: { t: 50, r: 50, b: 50, l: 60 }
        });
    });
}

// === ADD: Pathway coloring (with multi-pathway intersection) ===
const PATHWAY_PALETTE = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'];

function buildPathwayTraces(plotPoints) {
    const baseTrace = {
        x: [], y: [], text: [],
        mode: 'markers', type: 'scatter',
        name: 'Non-selected Genes',
        marker: { color: '#cccccc', size: 6, opacity: 0.7 },
        hoverinfo: 'skip',
        showlegend: true
    };

    const selected = selectedPathways.get ? selectedPathways.get() : new Set();
    const g2p = geneToPathwayMap.get ? geneToPathwayMap.get() : {};
    const pinfo = pathwayInfo.get ? pathwayInfo.get() : {};

    const pTraces = {};
    Array.from(selected).forEach((pid, idx) => {
        pTraces[pid] = {
            x: [], y: [], text: [], mode: 'markers', type: 'scatter',
            name: pinfo[pid]?.name || `Pathway ${pid}`,
            marker: { color: PATHWAY_PALETTE[idx % PATHWAY_PALETTE.length], size: 8, opacity: 0.9 },
            hovertemplate: '<b>%{text}</b><br>Pathway: ' + (pinfo[pid]?.name || pid) + '<extra></extra>',
            showlegend: true
        };
    });

    const multiTrace = {
        x: [], y: [], text: [], mode: 'markers', type: 'scatter',
        name: '🔴 Multi-Pathway Genes',
        marker: { color: 'red', size: 12, opacity: 0.85, symbol: 'circle', line: { color: 'red', width: 2 } },
        hovertemplate: '<b>%{text}</b><br>Found in multiple selected pathways<extra></extra>',
        showlegend: true
    };

    for (let i = 0; i < plotPoints.text.length; i++) {
        const gene = plotPoints.text[i];
        const x = plotPoints.x[i];
        const y = plotPoints.y[i];
        const pids = (g2p[gene] || []).filter(id => selected.has(String(id)));
        if (pids.length === 0) {
            baseTrace.x.push(x); baseTrace.y.push(y); baseTrace.text.push(gene);
        } else if (pids.length === 1) {
            const pid = String(pids[0]);
            const t = pTraces[pid];
            if (t) { t.x.push(x); t.y.push(y); t.text.push(gene); } else { baseTrace.x.push(x); baseTrace.y.push(y); baseTrace.text.push(gene); }
        } else {
            multiTrace.x.push(x); multiTrace.y.push(y); multiTrace.text.push(gene);
        }
    }

    const out = [baseTrace, ...Object.values(pTraces)];
    if (multiTrace.x.length > 0) out.push(multiTrace);
    return out;
}

function createComplexOverlayTrace(plotPoints, plotId) {
    const selected = selectedComplexes.get ? selectedComplexes.get() : new Set();
    const cmap = proteinComplexMap.get ? proteinComplexMap.get() : {};
    const out = {
        x: [], y: [], text: [], mode: 'markers', type: 'scatter',
        name: 'Protein Complex',
        marker: { symbol: 'star', color: 'gold', size: 16, line: { color: 'black', width: 1 } },
        hovertemplate: '<b>%{text}</b><br>(Protein complex member)<extra></extra>',
        showlegend: selected.size > 0
    };
    if (!selected || selected.size === 0) return out;

    const geneToIndex = window[plotId + '_geneToIndex'];
    if (!geneToIndex) return out;

    const genes = new Set();
    selected.forEach(name => (cmap[name] || []).forEach(g => genes.add(g)));
    genes.forEach(g => {
        const i = geneToIndex.get(g);
        if (i != null) {
            out.x.push(plotPoints.x[i]);
            out.y.push(plotPoints.y[i]);
            out.text.push(g);
        }
    });
    return out;
}

export function updateAllPlots() {
    const instances = plotInstances.get();
    const plots = plotData.get();
    if (!instances || !plots) return;
    instances.forEach((plotDiv) => {
        const plot = window[plotDiv.id + '_data'];
        if (!plot) return;

        let traces = [];
        if (coloringMode.get && coloringMode.get() === 'pathway') {
            traces = buildPathwayTraces(plot);
        } else {
            const result = currentSelectionResult.get();
            if (result && result.categories) {
                const categories = result.categories;
                const categoryNames = result.category_names || {};
                const inter = analyzeGeneIntersections(categories);
                traces = buildDynamicTraces(plot, categories, inter, categoryNames);
            } else {
                traces = [{ x: plot.x, y: plot.y, text: plot.text, mode: 'markers', type: 'scatter', name: 'Non-selected Genes', marker: { color: '#cccccc', size: 6, opacity: 0.7 } }];
            }
        }

        const overlay = createComplexOverlayTrace(plot, plotDiv.id);
        if (overlay.x.length > 0) traces.push(overlay);
        Plotly.react(plotDiv, traces, plotDiv.layout);
    });
}

// === Neighborhood coloring (BFS layers) ===
export function applyNeighborhoodColoring(centerGene, layers) {
    // layers: array of arrays, layers[0] = distance 1, layers[1] = distance 2, ...
    const instances = plotInstances.get();
    const plots = plotData.get();
    if (!instances || !plots) return;

    const palette = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf']; // red, blue, green, purple, orange, brown, pink

    instances.forEach((plotDiv, idx) => {
        const plot = window[plotDiv.id + '_data'];
        if (!plot) return;

        // Base grey trace
        const baseTrace = {
            x: [], y: [], text: [], mode: 'markers', type: 'scatter',
            name: 'Non-selected Genes', marker: { color: '#cccccc', size: 6, opacity: 0.7 },
            hovertemplate: '<b>%{text}</b><br>Log2FC: %{x}<br>-log10(p): %{y}<extra></extra>'
        };

        // Center gene trace (layer 0)
        const centerTrace = {
            x: [], y: [], text: [], mode: 'markers', type: 'scatter',
            name: `${centerGene} (center)`, marker: { color: 'black', size: 11, opacity: 1.0, symbol: 'diamond' }
        };

        const layerTraces = layers.map((_, i) => ({
            x: [], y: [], text: [], mode: 'markers', type: 'scatter',
            name: `Layer ${i+1}`, marker: { color: palette[i % palette.length], size: 9, opacity: 0.9 }
        }));

        // Build lookup sets for fast membership (convert to lowercase for case-insensitive matching)
        const layerSets = layers.map(arr => new Set((arr || []).map(g => g.toLowerCase())));

        for (let i = 0; i < plot.text.length; i++) {
            const g = plot.text[i];
            const x = plot.x[i];
            const y = plot.y[i];
            // CRITICAL FIX: Use case-insensitive comparison for gene matching
            if (g.toLowerCase() === centerGene.toLowerCase()) {
                centerTrace.x.push(x); centerTrace.y.push(y); centerTrace.text.push(g);
                continue;
            }
            let assigned = false;
            for (let li = 0; li < layerSets.length; li++) {
                if (layerSets[li].has(g.toLowerCase())) {
                    layerTraces[li].x.push(x); layerTraces[li].y.push(y); layerTraces[li].text.push(g);
                    assigned = true; break;
                }
            }
            if (!assigned) {
                baseTrace.x.push(x); baseTrace.y.push(y); baseTrace.text.push(g);
            }
        }

        const traces = [baseTrace];
        // Add only non-empty layers to help legend reflect content
        layerTraces.forEach(t => { if (t.x.length > 0) traces.push(t); });
        // Center last so it sits on top visually
        if (centerTrace.x.length > 0) traces.push(centerTrace);
        Plotly.react(plotDiv, traces, plotDiv.layout);
    });
}

// ===========================================================================
// LABEL MANAGEMENT
// ===========================================================================

// ==========================================================================================
// == NOTE FOR AI: THIS IS THE LABEL TOGGLE FUNCTION. IT MANAGES THE LABEL BUTTON STATE    ==
// == AND CALLS THE APPROPRIATE LABEL FUNCTIONS. DO NOT SIMPLIFY THE LOGIC.                ==
// ==========================================================================================
export function toggleGeneLabels() {
    console.log('🔍 Toggle gene labels called');
    
    const labelsButton = document.getElementById('toggle-labels');
    if (!labelsButton) {
        console.error('❌ Labels button not found');
        return;
    }
    
    const isActive = labelsButton.classList.contains('btn-primary');
    console.log('🔍 Labels currently active:', isActive);
    
    // === CRITICAL: TOGGLE LOGIC ===
    // NOTE FOR AI: This toggle logic is deliberate. We check the current state and
    // switch between active (btn-primary) and inactive (btn-outline-primary) states.
    // We also update the global state and call the appropriate label functions.
    if (isActive) {
        // === DEACTIVATE LABELS ===
        labelsButton.classList.remove('btn-primary');
        labelsButton.classList.add('btn-outline-primary');
        labelsVisible.set(false);  // Update global state
        removeGeneLabels();        // Remove all labels from plots
        console.log('🏷️ Labels deactivated');
    } else {
        // === ACTIVATE LABELS ===
        labelsButton.classList.remove('btn-outline-primary');
        labelsButton.classList.add('btn-primary');
        labelsVisible.set(true);   // Update global state
        addGeneLabelsForSelected(); // Add labels for selected genes
        console.log('🏷️ Labels activated');
    }
}

// addGeneLabels function removed - using addGeneLabelsForSelected instead

// Helper function for updating colors after selection
export function updatePlotColors(containerId, selectedGenes, categories) {
    console.log("🔍 Updating plot colors for:", containerId);
    console.log("🔍 Selected genes:", selectedGenes);
    console.log("🔍 Categories:", categories);
    
    const plotData = window[`${containerId}_data`];
    if (!plotData) {
        console.error(`❌ No plot data found for ${containerId}`);
        return;
    }
    
    // ALL GENES START AS LIGHT GREY - ONLY COLOR SELECTED ONES
    const greyPoints = { x: [], y: [], text: [] };
    const selectedPoints = { x: [], y: [], text: [] };
    
    plotData.text.forEach((geneName, index) => {
        const point = {
            x: plotData.x[index],
            y: plotData.y[index],
            text: plotData.text[index]
        };
        
        if (selectedGenes.includes(geneName)) {
            selectedPoints.x.push(point.x);
            selectedPoints.y.push(point.y);
            selectedPoints.text.push(point.text);
        } else {
            greyPoints.x.push(point.x);
            greyPoints.y.push(point.y);
            greyPoints.text.push(point.text);
        }
    });
    
    // Update the plot with new traces - ONLY GREY AND SELECTED
    const newTraces = [
        {
            x: greyPoints.x,
            y: greyPoints.y,
            text: greyPoints.text,
            mode: 'markers',
            type: 'scatter',
            name: 'Non-selected Genes',
            marker: { color: '#cccccc', size: 6, opacity: 0.7 }
        }
    ];
    
    // Add selected genes trace if there are any
    if (selectedPoints.x.length > 0) {
        newTraces.push({
            x: selectedPoints.x,
            y: selectedPoints.y,
            text: selectedPoints.text,
            mode: 'markers',
            type: 'scatter',
            name: 'Selected',
            marker: { color: 'blue', size: 8, opacity: 0.9 }
        });
    }
    
    Plotly.deleteTraces(containerId, [0, 1, 2]); // Remove old traces
    Plotly.addTraces(containerId, newTraces); // Add new traces
    
    console.log("✅ Plot colors updated successfully");
}

// removeGeneLabels function is exported later in the file

// ===========================================================================
// LEGEND MANAGEMENT
// ===========================================================================

export function toggleLegend() {
    const currentVisibility = legendVisible.get();
    legendVisible.set(!currentVisibility);
    
    const instances = plotInstances.get();
    if (!instances || instances.length === 0) return;
    
    instances.forEach(plotDiv => {
        Plotly.relayout(plotDiv, { showlegend: !currentVisibility });
    });
    
    // Update button state
    const button = document.getElementById('toggle-legend');
    if (button) {
        if (!currentVisibility) {
            button.classList.remove('btn-outline-secondary');
            button.classList.add('btn-secondary');
        } else {
            button.classList.remove('btn-secondary');
            button.classList.add('btn-outline-secondary');
        }
    }
}

// ===========================================================================
// TOOLTIP MANAGEMENT
// ===========================================================================

// Helper function to check if a gene is currently part of selected groups
function isGeneCurrentlySelected(geneName, traceName) {
    // CRITICAL: Only show tooltips for genes that are part of selected categories
    // The trace name indicates which category/group the gene belongs to
    
    // Always skip traces that represent non-selected or background genes
    const excludedTraceNames = [
        'Non-selected Genes',
        'Non-selected genes', 
        'Genes',
        'Background',
        'Unselected',
        'All Genes',
        'Base Genes'
    ];
    
    if (excludedTraceNames.includes(traceName)) {
        return false;
    }
    
    // Skip traces that are just generic gene collections
    if (traceName && traceName.includes('Genes') && !traceName.includes('(')) {
        return false;
    }
    
    // Allow tooltips for genes in colored traces (selected categories)
    // These traces represent genes that belong to selected groups/categories
    // Examples: "Group A (15)", "Pathway B (8)", "Regulon C (12)"
    if (traceName && !excludedTraceNames.includes(traceName)) {
        return true;
    }
    
    // Conservative approach: if we can't determine, don't show tooltip
    return false;
}

export function setupTooltipHandlers() {
    const instances = plotInstances.get();
    if (!instances || instances.length === 0) return;
    
    instances.forEach(instance => {
        // Resolve stored instance to a DOM element and ensure Plotly has initialized it
        const el = typeof instance === 'string' ? document.getElementById(instance) : instance;
        if (!el) return;
        if (typeof el.on !== 'function') {
            // Plotly has not augmented this element yet; skip to avoid runtime error
            console.warn('Tooltip handlers skipped: Plotly element not initialized for', el.id);
            return;
        }

        // Handle hover events with comprehensive debugging
        el.on('plotly_hover', (data) => {
            if (!isTooltipPinned.get()) {
                // PRE-VALIDATION: Check if this is a valid hover event with actual data points
                if (!data || !data.points || !data.points[0] || !data.points[0].trace) {
                    console.log('🔍 TOOLTIP: Invalid hover event, skipping tooltip');
                    return; // Skip invalid hover events
                }
                
                // Use comprehensive debugging framework
                geneDebugger.debugTooltipEvent(data);
                showTooltip(data, el);
            }
        });
        
        // Handle click events for pinning
        el.on('plotly_click', (data) => {
            if (isTooltipPinned.get() && pinnedGene.get() === data.points[0].text) {
                // Unpin if clicking the same gene
                unpinTooltip();
            } else {
                // Pin the tooltip
                pinTooltip(data.points[0].text);
                showTooltip(data, el);
            }
        });
        
        // Handle unhover events
        el.on('plotly_unhover', () => {
            if (!isTooltipPinned.get()) {
                hideTooltip();
            }
        });
    });
}

function showTooltip(data, plotDiv) {
    // COMPREHENSIVE VALIDATION with debugging framework
    if (!data || !data.points || !data.points[0]) {
        geneDebugger.log("🔍 TOOLTIP: No valid data points", 'critical');
        return; // Don't show tooltip for invalid data
    }
    
    const point = data.points[0];
    
    // CRITICAL: Validate that we have a gene name
    if (!point.text || typeof point.text !== 'string') {
        geneDebugger.log("🔍 TOOLTIP: No valid gene name", 'critical');
        return; // Don't show tooltip without gene name
    }
    
    const geneName = point.text;
    const x = point.x;
    const y = point.y;
    
    // ENHANCED COORDINATE VALIDATION
    if (typeof x !== 'number' || typeof y !== 'number' || isNaN(x) || isNaN(y)) {
        geneDebugger.log("🔍 TOOLTIP: Invalid coordinates", { x, y }, 'critical');
        return; // Don't show tooltip if coordinates are invalid
    }
    
    // ENHANCED TRACE VALIDATION: ensure we're hovering over a scatter trace, not a shape
    if (!point.trace || point.trace.type !== 'scatter') {
        geneDebugger.log("🔍 TOOLTIP: Not a scatter trace", point.trace?.type, 'critical');
        return; // Don't show tooltip for non-scatter elements (like threshold lines)
    }
    
    // ADDITIONAL VALIDATION: Check if we're hovering over actual data, not empty space
    if (point.trace.hoverinfo === 'skip' || point.trace.hoverinfo === 'none') {
        geneDebugger.log("🔍 TOOLTIP: Trace has hoverinfo set to skip/none", point.trace.hoverinfo, 'critical');
        return; // Don't show tooltip for elements that should not have hover
    }
    
    // VALIDATION: Ensure the trace has proper hover configuration
    if (!point.trace.hovertemplate && !point.trace.hoverinfo) {
        geneDebugger.log("🔍 TOOLTIP: Trace has no hover configuration", 'critical');
        return; // Don't show tooltip for traces without proper hover setup
    }
    
    // DYNAMIC VALIDATION: Check if the trace is currently visible
    // This respects user actions like hiding traces via legend clicks
    if (point.trace.visible === false || point.trace.visible === 'legendonly') {
        geneDebugger.log("🔍 TOOLTIP: Trace is hidden by user", { trace: point.trace.name, visible: point.trace.visible }, 'critical');
        return; // Don't show tooltip for hidden traces
    }
    
    // VALIDATION: Ensure we're hovering over a real data point with valid coordinates
    if (point.x === undefined || point.y === undefined || 
        point.x === null || point.y === null || 
        isNaN(point.x) || isNaN(point.y)) {
        geneDebugger.log("🔍 TOOLTIP: Invalid point coordinates", { x: point.x, y: point.y }, 'critical');
        return; // Don't show tooltip for invalid coordinates
    }
    
    // CRITICAL: Only show tooltips for colored (selected) genes
    // Check if this trace is the base grey trace (non-selected genes)
    const isGreyTrace = point.trace.name === 'Non-selected Genes';
    
    if (isGreyTrace) {
        geneDebugger.log("🔍 TOOLTIP: Skipping grey genes", point.trace.name, 'critical');
        return; // Don't show tooltip for grey (non-selected) genes
    }
    
    // DYNAMIC VALIDATION: Check if this gene is currently part of selected groups
    // This ensures tooltips only appear for genes that are actually selected/visible
    const shouldShowTooltip = isGeneCurrentlySelected(geneName, point.trace.name);
    if (!shouldShowTooltip) {
        geneDebugger.log("🔍 TOOLTIP: Gene not currently selected", { gene: geneName, trace: point.trace.name }, 'critical');
        return; // Don't show tooltip for genes that are not currently selected
    }
    
    geneDebugger.log("🔍 TOOLTIP: Gene is currently selected", { gene: geneName, trace: point.trace.name }, 'info');
    
    // ADDITIONAL VALIDATION: Check if the point is actually visible (not hidden by zoom/pan)
    if (point.curveNumber === undefined || point.pointNumber === undefined) {
        geneDebugger.log("🔍 TOOLTIP: Point not properly indexed", 'critical');
        return;
    }
    
    // EXTRA VALIDATION: Ensure we're not over empty space
    if (!point.trace.x || !point.trace.y || !point.trace.text) {
        geneDebugger.log("🔍 TOOLTIP: Trace data missing", 'critical');
        return;
    }
    
    // FINAL VALIDATION: Check if this is actually a real data point
    if (point.pointNumber >= point.trace.x.length || 
        point.pointNumber >= point.trace.y.length || 
        point.pointNumber >= point.trace.text.length) {
        geneDebugger.log("🔍 TOOLTIP: Point number out of bounds", 'critical');
        return;
    }
    
    // Verify the point data matches the trace data
    if (point.trace.x[point.pointNumber] !== x || 
        point.trace.y[point.pointNumber] !== y || 
        point.trace.text[point.pointNumber] !== geneName) {
        geneDebugger.log("🔍 TOOLTIP: Point data mismatch with trace data", 'critical');
        return;
    }
    
    // EXTRA VALIDATION: Ensure the point is not just a hover event from empty space
    // Check if the point coordinates are within reasonable bounds of the actual data
    const traceX = point.trace.x[point.pointNumber];
    const traceY = point.trace.y[point.pointNumber];
    const tolerance = 0.001; // Small tolerance for floating point precision
    
    if (Math.abs(traceX - x) > tolerance || Math.abs(traceY - y) > tolerance) {
        geneDebugger.log("🔍 TOOLTIP: Point coordinates don't match trace data within tolerance", 
            { pointX: x, pointY: y, traceX, traceY, tolerance }, 'critical');
        return; // Don't show tooltip if coordinates don't match
    }
    
    // FINAL VALIDATION: Check if the point is actually within the plot's data range
    const traceXArray = point.trace.x;
    const traceYArray = point.trace.y;
    if (traceXArray && traceYArray) {
        const minX = Math.min(...traceXArray);
        const maxX = Math.max(...traceXArray);
        const minY = Math.min(...traceYArray);
        const maxY = Math.max(...traceYArray);
        
        // Add a small buffer to the range check
        const buffer = 0.1;
        if (x < minX - buffer || x > maxX + buffer || y < minY - buffer || y > maxY + buffer) {
            geneDebugger.log("🔍 TOOLTIP: Point outside trace data range", 
                { x, y, minX, maxX, minY, maxY, buffer }, 'critical');
            return; // Don't show tooltip for points outside the data range
        }
    }
    
    geneDebugger.log("✅ TOOLTIP: Showing tooltip for", geneName, "in trace", point.trace.name, 'critical');
    
    // Create tooltip content
    const tooltipContent = `
        <div class="tooltip-content">
            <strong>${geneName}</strong><br>
            Log2 FC: ${x.toFixed(3)}<br>
            -log10(p): ${y.toFixed(3)}
        </div>
    `;
    
    // Show tooltip
    const tooltip = document.createElement('div');
    tooltip.id = 'gene-tooltip';
    tooltip.className = 'tooltip';
    tooltip.innerHTML = tooltipContent;
    tooltip.style.cssText = `
        position: absolute;
        background: rgba(0,0,0,0.8);
        color: white;
        padding: 8px;
        border-radius: 4px;
        font-size: 12px;
        pointer-events: none;
        z-index: 1000;
        left: ${data.event.clientX + 10}px;
        top: ${data.event.clientY - 10}px;
    `;
    
    document.body.appendChild(tooltip);
}

function hideTooltip() {
    const tooltip = document.getElementById('gene-tooltip');
    if (tooltip) {
        tooltip.remove();
    }
}

function pinTooltip(geneName) {
    pinnedGene.set(geneName);
    isTooltipPinned.set(true);
}

function unpinTooltip() {
    pinnedGene.set(null);
    isTooltipPinned.set(false);
    hideTooltip();
}

// ===========================================================================
// PLOT CONTROLS
// ===========================================================================

export function downloadPlotAsPNG() {
    const instances = plotInstances.get();
    if (!instances || instances.length === 0) return;
    
    // Download each plot as PNG
    instances.forEach((plotDiv, index) => {
        const plots = plotData.get();
        const plot = plots[index];
        const filename = `${plot.title || `plot_${index + 1}`}.png`;
        
        Plotly.downloadImage(plotDiv, {
            format: 'png',
            filename: filename,
            height: 600,
            width: 800
        });
    });
}

export function resetPlotScale() {
    const instances = plotInstances.get();
    if (!instances || instances.length === 0) return;
    
    instances.forEach(plotDiv => {
        Plotly.relayout(plotDiv, {
            'xaxis.autorange': true,
            'yaxis.autorange': true
        });
    });
}

// ===========================================================================
// PLOT INITIALIZATION
// ===========================================================================

export function initializePlots() {
    console.log('🚀 INITIALIZING PLOTS...');
    
    // GUARD: Prevent multiple initializations
    const existingInstances = plotInstances.get();
    if (existingInstances && existingInstances.length > 0) {
        console.log('⚠️ Plots already initialized, skipping duplicate initialization');
        return;
    }
    
    const plotContainer = document.getElementById('plot-container');
    if (!plotContainer) {
        console.error('❌ plot-container element not found');
        return;
    }
    
    // Get plot data from window object (set by Flask template)
    const plots = window.plotData || [];
    console.log('🔍 Plot data available:', plots.length, 'plots');
    console.log('🔍 Plot data:', plots);
    
    if (plots.length === 0) {
        console.warn('⚠️ No plot data available. Creating dummy plot for testing.');
        // Create a dummy plot if no data (for development/testing)
        plotContainer.innerHTML = '<div class="alert alert-warning">No plot data available. Please upload data first.</div>';
        return;
    }
    
    // Store plot data
    plotData.set(plots);
    
    // Render each plot
    plots.forEach((plotData, index) => {
        const plotDiv = document.createElement('div');
        plotDiv.id = `plot-${index}`;
        plotDiv.style.width = '100%';
        // Match the layout height to avoid cropping/misalignment
        plotDiv.style.height = '600px';
        plotDiv.style.marginBottom = '20px';
        plotDiv.style.position = 'relative'; // Ensure proper positioning for toolbar
        
        plotContainer.appendChild(plotDiv);
        
        // Convert raw data to Plotly format
        console.log('🔍 Creating volcano plot with data:', plotData);
        console.log('🔍 Plot title from data:', plotData.title);
        console.log('🔍 Colors in data:', plotData.colors);
        
        // Count color distribution for debugging
        const colorCounts = {};
        plotData.colors.forEach(color => {
            colorCounts[color] = (colorCounts[color] || 0) + 1;
        });
        console.log('🔍 Color distribution:', colorCounts);
        
        // ALL GENES ARE LIGHT GREY BY DEFAULT - NO RED GENES
        const greyPoints = { 
            x: plotData.x, 
            y: plotData.y, 
            text: plotData.text, 
            mode: 'markers', 
            type: 'scatter', 
            name: 'Genes',
            marker: { color: '#cccccc', size: 6, opacity: 0.7 },
            hovertemplate: '<b>%{text}</b><br>Log2FC: %{x}<br>-log10(p): %{y}<extra></extra>'
        };
        
        const plotlyData = [greyPoints];
        
        console.log('🔍 All genes are light grey:', greyPoints.x.length);
        
        // ADD PROMINENT DOTTED THRESHOLD LINES
        const thresholdLines = [];
        const maxY = Math.max(...plotData.y) * 1.15;
        const minX = Math.min(...plotData.x) * 1.15;
        const maxX = Math.max(...plotData.x) * 1.15;
        
        // VERTICAL DOTTED LINES for fold change thresholds
        if (plotData.fc_thresh) {
            // Positive fold change threshold
            thresholdLines.push({
                type: 'line',
                x0: plotData.fc_thresh,
                x1: plotData.fc_thresh,
                y0: 0,
                y1: maxY,
                line: { 
                    color: '#FF6B6B', 
                    width: 3, 
                    dash: 'dot' 
                },
                name: `Log2FC ≥ ${plotData.fc_thresh}`,
                layer: 'below',
                hoverinfo: 'skip',
                hoveron: 'points', // Only hover on actual data points, not on the line
                hovermode: 'closest' // Use closest mode for this specific shape
            });
            
            // Negative fold change threshold
            thresholdLines.push({
                type: 'line',
                x0: -plotData.fc_thresh,
                x1: -plotData.fc_thresh,
                y0: 0,
                y1: maxY,
                line: { 
                    color: '#FF6B6B', 
                    width: 3, 
                    dash: 'dot' 
                },
                name: `Log2FC ≤ -${plotData.fc_thresh}`,
                layer: 'below',
                hoverinfo: 'skip',
                hoveron: 'points', // Only hover on actual data points, not on the line
                hovermode: 'closest' // Use closest mode for this specific shape
            });
        }
        
        // HORIZONTAL DOTTED LINE for p-value threshold
        if (plotData.p_thresh_log) {
            thresholdLines.push({
                type: 'line',
                x0: minX,
                x1: maxX,
                y0: plotData.p_thresh_log,
                y1: plotData.p_thresh_log,
                line: { 
                    color: '#4ECDC4', 
                    width: 3, 
                    dash: 'dot' 
                },
                name: `p-value ≤ ${Math.pow(10, -plotData.p_thresh_log).toExponential(1)}`,
                layer: 'below',
                hoverinfo: 'skip',
                hoveron: 'points', // Only hover on actual data points, not on the line
                hovermode: 'closest' // Use closest mode for this specific shape
            });
        }
        
        // Ensure we have a proper title
        const plotTitle = plotData.title && plotData.title.trim() !== '' ? plotData.title : 'Volcano Plot';
        console.log('🔧 Setting plot title:', plotTitle);
        
        const layout = {
            title: { text: plotTitle },
            xaxis: { title: { text: (plotData.x_title || 'Log2 Fold Change') }, zeroline: true, zerolinecolor: '#000000', zerolinewidth: 1 },
            yaxis: { title: { text: (plotData.y_title || '-log10(p-value)') }, zeroline: false },
            hovermode: 'closest', // Changed from 'x unified' to 'closest' for more precise tooltip behavior
            showlegend: true,
            width: null,
            height: 600,
            margin: { t: 50, r: 50, b: 50, l: 60 },
            shapes: thresholdLines
        };
        
        // Create a safe filename for downloads
        const safeTitle = plotTitle.replace(/[^a-zA-Z0-9\s-_]/g, '').replace(/\s+/g, '_');
        const downloadFilename = `${safeTitle}_volcano_plot`;
        
        const config = { 
            responsive: true, 
            displayModeBar: true, 
            displaylogo: false,
            editable: true,  // Enable global editing for annotations
            toImageButtonOptions: {
                format: 'png',
                filename: downloadFilename,
                height: 600,
                width: 800,
                scale: 2
            },
            modeBarButtonsToAdd: [
                'pan2d',
                'select2d',
                'lasso2d',
                'resetScale2d',
                'zoomIn2d',
                'zoomOut2d',
                'autoScale2d',
                
                'toggleSpikelines'
            ],
            modeBarButtonsToRemove: [
                'sendDataToCloud',
                'editInChartStudio'
            ]
        };
        
                            // Create the plot
                            console.log('🔧 Creating plot with config:', config);
                            return Plotly.newPlot(plotDiv.id, plotlyData, layout, config)
                        .then(() => {
                            console.log('✅ Plot created successfully in', plotDiv.id);
                            
                            // Ensure toolbar is properly initialized
                            const plotElement = document.getElementById(plotDiv.id);
                            if (plotElement && plotElement._fullLayout) {
                                console.log('🔧 Plot layout initialized, toolbar should be available');
                            }
                            
                            // Store plot data for labels and other functions
                            window[plotDiv.id + '_data'] = plotData;
                            console.log(`🔍 STORED PLOT DATA for ${plotDiv.id}:`, {
                                geneCount: plotData.text.length,
                                firstFewGenes: plotData.text.slice(0, 5),
                                dataKeys: Object.keys(plotData)
                            });
                            
                            // Per-plot gene index map for O(1) lookups
                            const geneToIndex = new Map();
                            plotData.text.forEach((g, i) => geneToIndex.set(g, i));
                            window[plotDiv.id + '_geneToIndex'] = geneToIndex;

                            // Store plot instance
                            const instances = plotInstances.get();
                            instances.push(plotDiv);
                            plotInstances.set(instances);

                                    // Setup tooltip handlers
        setupTooltipHandlers();
        
        // Setup legend click handlers to update labels when visibility changes
        setupLegendClickHandlers(plotDiv.id);

        // Setup draggable labels
        setupDraggableLabels(plotDiv.id);

                            // Hide loading message after all plots are created
                            if (index === plots.length - 1) {
                                const loadingElement = document.getElementById('plot-loading');
                                if (loadingElement) {
                                    loadingElement.style.display = 'none';
                                }
                            }
                        })
                        .catch(error => {
                            console.error('❌ Error creating plot:', error);
                        });
    });
}

// PLOT CONTROL FUNCTIONS (toggleGeneLabels is declared earlier in the file)

// ==========================================================================================
// == NOTE FOR AI: THIS IS THE GENE LABELING LOGIC. IT IS COMPLEX AND DELIBERATE.          ==
// == DO NOT SIMPLIFY OR CHANGE THE CORE LOGIC. The goal is to label ONLY the genes        ==
// == that belong to CATEGORIES the user has selected AND that are currently VISIBLE       ==
// == on the plot (i.e., not toggled off in the legend).                                   ==
// ==========================================================================================
export function addGeneLabelsForSelected() {
    // Get all plot containers (e.g., plot-0, plot-1)
    const plotContainers = document.querySelectorAll('[id^="plot-"]:not(#plot-container):not(#plot-loading):not(#plot-data)');
    if (!plotContainers.length) {
        console.error('LABEL ERROR: No plot containers found.');
        return;
    }

    // === CRITICAL STEP 1: Get Genes from VISIBLE Traces ONLY ===
    // NOTE FOR AI: This section is the most important. We build a Set of genes to label
    // by iterating through each plot's traces and checking two conditions:
    // 1. The trace MUST NOT be the base trace of non-selected genes.
    // 2. The trace MUST be visible (i.e., trace.visible is not false or 'legendonly').
    // This respects the user's actions of hiding specific groups via the legend.
    
    plotContainers.forEach(container => {
        const plotId = container.id;
        const plotElement = document.getElementById(plotId);
        const plotData = window[plotId + '_data']; // The original full dataset for coordinates

        if (!plotElement || !plotElement.data || !plotData) {
            console.error(`LABEL ERROR: Plot data or element missing for ${plotId}`);
            return;
        }

        const genesToLabel = new Set();
        
        plotElement.data.forEach(trace => {
            // === AI GUARD: DO NOT REMOVE THIS CHECK ===
            // We explicitly skip the main grey trace of background genes.
            if (trace.name === 'Non-selected Genes' || trace.name === 'Non-selected genes' || trace.name === 'Genes') {
                console.log(`🏷️ LABEL: Skipping base trace "${trace.name}"`);
                return; // Continue to the next trace
            }
            
            // === AI GUARD: DO NOT REMOVE THIS CHECK ===
            // We only consider traces that are currently visible. This respects legend clicks.
            if (trace.visible !== false && trace.visible !== 'legendonly') {
                if (trace.text && Array.isArray(trace.text)) {
                    trace.text.forEach(gene => genesToLabel.add(gene));
                    console.log(`🏷️ LABEL: Adding ${trace.text.length} genes from visible trace "${trace.name}"`);
                }
            } else {
                console.log(`🏷️ LABEL: Skipping hidden trace "${trace.name}" (visible: ${trace.visible})`);
            }
        });

        if (genesToLabel.size === 0) {
            console.warn(`LABEL INFO: No selected or visible genes to label on ${plotId}. This is normal if no categories are selected.`);
            // It's not an error if no genes are selected, so we just stop here for this plot.
            return; 
        }

        // === CRITICAL STEP 2: Create Annotations for Matched Genes ===
        // NOTE FOR AI: Now we iterate through the ORIGINAL plot data to find the coordinates
        // of the genes we collected in the 'genesToLabel' set.
        
        const annotations = [];
        // Use a Set for the original plot genes for faster lookups if needed, though linear scan is fine.
        const allPlotGenes = new Map(plotData.text.map((gene, index) => [gene, index]));

        genesToLabel.forEach(geneName => {
            if (allPlotGenes.has(geneName)) {
                const index = allPlotGenes.get(geneName);
                annotations.push({
                    x: plotData.x[index],
                    y: plotData.y[index],
                    text: geneName,
                    showarrow: true,
                    arrowhead: 2,
                    arrowcolor: '#333',
                    font: { size: 10 },
                    bgcolor: 'rgba(255, 255, 255, 0.8)',
                    borderwidth: 1,
                    // NOTE FOR AI: The 'ax' and 'ay' pixel offsets are crucial for preventing
                    // the label from perfectly overlapping the data point. Do not remove.
                    ax: 0,
                    ay: -25,
                    axref: 'pixel',
                    ayref: 'pixel',
                });
            }
        });

        // === CRITICAL STEP 3: Apply Annotations via Plotly.relayout ===
        // NOTE FOR AI: We use Plotly.relayout to add the annotations without redrawing the
        // entire plot. The 'editable: true' property is ESSENTIAL for making the labels
        // draggable by the user. DO NOT set it to false or remove it.
        
        // Use debugging framework to analyze label creation
        geneDebugger.debugLabelCreation(plotId, genesToLabel, annotations);
        
        Plotly.relayout(plotId, {
            annotations: annotations,
            editable: true,
        });
        
        console.log(`🏷️ LABELS: Added ${annotations.length} draggable labels to ${plotId}`);
    });
}

// Setup legend click handlers to update labels when trace visibility changes
function setupLegendClickHandlers(plotId) {
    const plotElement = document.getElementById(plotId);
    if (!plotElement) return;
    
    // Listen for legend clicks to update labels
    plotElement.on('plotly_legendclick', function(data) {
        console.log(`🏷️ LEGEND: Legend item clicked for "${data.data[data.curveNumber].name}"`);
        
        // Update labels after a short delay to allow the trace visibility to change
        setTimeout(() => {
            if (labelsVisible.get()) {
                console.log('🏷️ LEGEND: Updating labels after legend click');
                addGeneLabelsForSelected();
            }
        }, 100);
    });
}

function setupDraggableLabels(plotId) {
    // Use comprehensive debugging framework for draggable labels
    geneDebugger.debugDraggableLabels(plotId);
    
    const plotDiv = document.getElementById(plotId);
    if (!plotDiv) {
        console.error(`❌ Plot element ${plotId} not found for draggable labels setup`);
        return;
    }
    
    // Remove existing event handlers to prevent duplicates
    if (plotDiv.removeAllListeners) {
        plotDiv.removeAllListeners('plotly_relayout');
        // REMOVED: plotDiv.removeAllListeners('plotly_hover');  // ← This was too aggressive and broke tooltips
        plotDiv.removeAllListeners('plotly_click');
    }
    
    // Set up event handlers for draggable annotations with debugging
    plotDiv.on('plotly_relayout', function(eventData) {
        // Use debugging framework to analyze relayout events
        geneDebugger.debugLabelDragEvent(eventData, plotId);
        
        console.log('📍 DRAGGABLE: Relayout event detected:', Object.keys(eventData));
        
        // Check if this is an annotation drag event
        let annotationMoved = false;
        Object.keys(eventData).forEach(key => {
            if (key.startsWith('annotations[') && (key.endsWith('].x') || key.endsWith('].y') || key.endsWith('].ax') || key.endsWith('].ay'))) {
                const annotationIndex = key.match(/\[(\d+)\]/)[1];
                annotationMoved = true;
                console.log(`📍 DRAGGABLE: Label ${annotationIndex} moved in ${plotId}`);
                
                // Get current annotation position
                const annotations = plotDiv.layout.annotations;
                if (annotations && annotations[annotationIndex]) {
                    const annotation = annotations[annotationIndex];
                    console.log(`🎯 DRAGGABLE: ${annotation.text} moved to (${annotation.x}, ${annotation.y})`);
                }
            }
        });
        
        if (annotationMoved) {
            console.log('✅ DRAGGABLE: Annotation successfully moved!');
        }
    });
    
    // REMOVED: Conflicting hover handler that was causing tooltip issues
    // The main tooltip handler in setupTooltipHandlers() should be the only one
    
    plotDiv.on('plotly_click', function(data) {
        if (data && data.points && data.points[0]) {
            const clickedGene = data.points[0].text;
            console.log(`🖱️ DRAGGABLE: Clicked on gene: ${clickedGene}`);
        }
    });
    
    // Additional event for annotation editing specifically
    plotDiv.on('plotly_beforeplot', function() {
        console.log(`🎮 DRAGGABLE: Plot ${plotId} ready for annotation editing`);
    });
    
    console.log(`🎮 DRAGGABLE: Enhanced draggable label controls enabled for ${plotId}`);
}

// ==========================================================================================
// == NOTE FOR AI: THIS IS THE LABEL REMOVAL FUNCTION. IT CLEARS ALL ANNOTATIONS FROM      ==
// == ALL PLOTS. DO NOT MODIFY THE CORE LOGIC.                                            ==
// ==========================================================================================
export function removeGeneLabels() {
    const plotContainers = document.querySelectorAll('[id^="plot-"]:not(#plot-container):not(#plot-data):not(#plot-loading)');
    console.log('🚫 LABELS: Removing labels from plots:', Array.from(plotContainers).map(c => c.id));
    
    plotContainers.forEach(container => {
        try {
            // === CRITICAL: SAFE LABEL REMOVAL ===
            // NOTE FOR AI: We check if the plot exists before trying to modify it.
            // We use Plotly.relayout to clear annotations without redrawing the entire plot.
            // We keep editable: true to maintain the plot's interactive state.
            if (container && document.getElementById(container.id)) {
                // Use debugging framework to analyze label removal
                geneDebugger.debugLabelRemoval(container.id);
                
                Plotly.relayout(container.id, { 
                    annotations: [],  // Clear all annotations
                    editable: true   // Keep editable mode for future label additions
                });
                console.log(`🚫 LABELS: Removed labels from ${container.id}`);
            }
        } catch (error) {
            console.warn(`⚠️ LABELS: Could not remove labels from ${container.id}:`, error.message);
        }
    });
    console.log('🚫 LABELS: Gene labels removal completed');
}