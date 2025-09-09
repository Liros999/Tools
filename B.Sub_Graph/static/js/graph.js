import cytoscape from 'https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.28.1/cytoscape.esm.min.js';
import cola from 'https://unpkg.com/cytoscape-cola/cytoscape-cola.js';
import centrality from 'https://unpkg.com/cytoscape-centrality@^1.0.0/cytoscape-centrality.js';

cytoscape.use(cola);
cytoscape.use(centrality);

/**
 * Manages all graph rendering and interaction using Cytoscape.js.
 */
export class CytoscapeManager {
    constructor(containerId) {
        this.container = document.getElementById(containerId);
        if (!this.container) throw new Error(`Container with id "${containerId}" not found.`);
        
        this.cy = cytoscape({
            container: this.container,
            style: this._getStylesheet(),
            ...this._optimizeForLargeGraphs()
        });
    }

    _getStylesheet() {
        return [
            {
                selector: 'node',
                style: {
                    'background-color': '#ffffff',
                    'border-color': '#666',
                    'border-width': 2,
                    'label': 'data(name)',
                    'width': 'mapData(degree, 0, 50, 20, 60)',
                    'height': 'mapData(degree, 0, 50, 20, 60)',
                    'font-size': '10px',
                    'color': '#333',
                    'text-valign': 'bottom',
                    'text-halign': 'center',
                    'text-margin-y': '5px'
                }
            },
            {
                selector: 'edge',
                style: {
                    'width': 2,
                    'line-color': '#ccc',
                    'target-arrow-color': '#ccc',
                    'target-arrow-shape': 'triangle',
                    'curve-style': 'bezier'
                }
            },
            {
                selector: 'edge[interaction_type="regulation"]',
                style: { 'line-color': '#e74c3c', 'target-arrow-color': '#e74c3c' }
            },
            {
                selector: 'edge[interaction_type="binding"]',
                style: { 'line-color': '#3498db', 'target-arrow-color': '#3498db', 'line-style': 'dashed' }
            },
            {
                selector: 'edge[interaction_type="metabolic"]',
                style: { 'line-color': '#2ecc71', 'target-arrow-color': '#2ecc71' }
            },
            {
                selector: '.search-highlight',
                style: {
                    'background-color': '#f1c40f',
                    'border-color': '#f39c12',
                    'border-width': 4
                }
            }
        ];
    }

    _optimizeForLargeGraphs() {
        return { wheelSensitivity: 0.2, textureOnViewport: true, motionBlur: true };
    }

    getOptimalLayout(nodeCount) {
        if (nodeCount < 1000) {
            return { name: 'cose', idealEdgeLength: 100, nodeRepulsion: 400000, edgeElasticity: 100, numIter: 1000 };
        } else {
            return { name: 'cola', edgeLength: 70, animate: false, randomize: false, maxSimulationTime: 3000 };
        }
    }

    render(molecules, interactions) {
        const elements = [];
        const moleculeIds = new Set();
        molecules.forEach(mol => {
            if (mol.type === 'Protein') {
                elements.push({ group: 'nodes', data: { id: mol.id, name: mol.name } });
                moleculeIds.add(mol.id);
            }
        });
        interactions.forEach(intr => {
            if (intr.molecule_ids.length === 2) {
                const [source, target] = intr.molecule_ids;
                if (moleculeIds.has(source) && moleculeIds.has(target)) {
                    elements.push({ group: 'edges', data: { id: intr.id, source, target, interaction_type: (intr.type || 'unknown').toLowerCase() } });
                }
            }
        });
        this.cy.elements().remove();
        this.cy.add(elements);

        // Add degree data to nodes for styling
        this.cy.nodes().forEach(node => {
            node.data('degree', node.degree());
        });

        this.runLayout();
    }

    runLayout() {
        const nodeCount = this.cy.nodes().length;
        const layoutOptions = this.getOptimalLayout(nodeCount);
        const layout = this.cy.layout(layoutOptions);
        layout.run();
        this.cy.fit(50);
    }

    updateNodeColors(dataMap, colorScale) {
        this.cy.nodes().forEach(node => {
            const geneName = node.data('name');
            let color = '#ffffff';
            if (dataMap.has(geneName)) {
                color = colorScale(dataMap.get(geneName));
            }
            node.style('background-color', color);
        });
    }
    
highlightNodesByName(geneNames) {
        this.cy.elements().removeClass('search-highlight');
        if (!geneNames || geneNames.length === 0) return;

        const selector = geneNames.map(name => `node[name = "${name}"]`).join(',');
        const highlightedNodes = this.cy.nodes(selector);
        
        highlightedNodes.addClass('search-highlight');
        this.cy.animate({ fit: { eles: highlightedNodes, padding: 100 }, duration: 500 });
    }

    // --- Network Metrics --- //
    calculateCentrality(type) {
        let centralityResult;
        switch(type) {
            case 'degree':
                // Degree is already calculated and stored in node.data('degree')
                return;
            case 'betweenness':
                centralityResult = this.cy.elements().betweennessCentrality();
                break;
            case 'closeness':
                centralityResult = this.cy.elements().closenessCentrality();
                break;
            default:
                console.warn(`Unknown centrality type: ${type}`);
                return;
        }
        
        this.cy.nodes().forEach(node => {
            node.data(type, centralityResult.node(node));
        });
    }

    getCentralityValue(node, type) {
        return node.data(type);
    }

    // Tooltip and event binding methods would go here...

    export(format, options = {}) {
        const exportOptions = { 
            bg: options.backgroundColor || 'white', 
            full: options.includeHidden || false, 
            scale: options.scale || 2 
        };

        let content;
        let filename;
        let mimeType;

        switch(format) {
            case 'png':
                content = this.cy.png(exportOptions);
                filename = `subtiwiki_network_${Date.now()}.png`;
                mimeType = 'image/png';
                break;
            case 'svg':
                content = this.cy.svg(exportOptions);
                filename = `subtiwiki_network_${Date.now()}.svg`;
                mimeType = 'image/svg+xml';
                break;
            case 'cytoscape':
                content = JSON.stringify(this.cy.json(), null, 2);
                filename = `subtiwiki_network_${Date.now()}.cyjs`;
                mimeType = 'application/json';
                break;
            case 'sif':
                content = this.exportSIFFormat();
                filename = `subtiwiki_network_${Date.now()}.sif`;
                mimeType = 'text/plain';
                break;
            default:
                console.error('Unsupported export format:', format);
                return;
        }
        this._downloadFile(content, filename, mimeType);
    }

    exportSIFFormat() {
        const edges = this.cy.edges();
        const sifLines = [];
        edges.forEach(edge => {
            const source = edge.source().data('name');
            const target = edge.target().data('name');
            const interaction = edge.data('interaction_type') || 'interacts_with';
            sifLines.push(`${source}\t${interaction}\t${target}`);
        });
        return sifLines.join('\n');
    }

    _downloadFile(content, filename, mimeType) {
        const blob = new Blob([content], { type: mimeType });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        link.click();
        URL.revokeObjectURL(url);
    }
}

function createTooltipHtml(context) {
    let html = `<div class="tooltip-content"><strong>${context.name}</strong>`;
    if (context.function) html += `<br><em>${context.function}</em>`;
    html += '<hr>';
    if (context.synonyms && context.synonyms.length > 0) {
        html += `<strong>Synonyms:</strong> ${context.synonyms.join(', ')}<br>`;
    }
    if (context.categories && context.categories.length > 0) {
        html += `<strong>Categories:</strong> ${context.categories.map(c => c.name).join(', ')}`;
    }
    html += `</div>`;
    return html;
}
