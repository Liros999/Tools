import cytoscape from 'https://unpkg.com/cytoscape@3.28.1/dist/cytoscape.esm.min.js';
import dagre from 'https://unpkg.com/dagre@0.8.5/dist/dagre.min.js';
import cyDagre from 'https://unpkg.com/cytoscape-dagre@2.5.0/cytoscape-dagre.esm.js';

cytoscape.use(cyDagre);

export class CytoscapeManager {
    constructor(containerId = 'cy') {
        this.cy = cytoscape({
            container: document.getElementById(containerId),
            elements: [],
            layout: { name: 'dagre', rankSep: 60, nodeSep: 40, edgeSep: 20 },
            style: [
                { selector: 'node', style: { 'label': 'data(label)', 'font-size': 10, 'text-valign': 'center', 'color': '#000', 'background-color': '#fff', 'border-width': 1, 'border-color': '#777' } },
                { selector: 'edge', style: { 'width': 1.5, 'line-color': '#999', 'curve-style': 'bezier', 'target-arrow-shape': 'triangle', 'target-arrow-color': '#999' } }
            ]
        });
    }

    draw(molecules, interactions) {
        const proteinSet = new Set(molecules.filter(m => m.type === 'Protein').map(m => m.name));
        const nodes = [...proteinSet].map(name => ({ data: { id: name, label: name } }));
        const edges = [];
        const idToName = new Map(molecules.map(m => [m.id, m.name]));
        for (const it of interactions) {
            if (Array.isArray(it.molecule_ids) && it.molecule_ids.length === 2) {
                const a = idToName.get(it.molecule_ids[0]);
                const b = idToName.get(it.molecule_ids[1]);
                if (a && b && proteinSet.has(a) && proteinSet.has(b)) {
                    edges.push({ data: { id: `${a}__${b}__${it.id}`, source: a, target: b } });
                }
            }
        }
        this.cy.elements().remove();
        this.cy.add(nodes).add(edges);
        this.cy.layout({ name: 'dagre', rankSep: 60, nodeSep: 40, edgeSep: 20 }).run();
    }

    applyColorsFromMap(dataMap) {
        const values = Object.values(dataMap || {});
        if (!values.length) return this.resetColors();
        const maxAbs = Math.max(...values.map(v => Math.abs(v)));
        if (maxAbs === 0) return this.resetColors();
        const min = -maxAbs, max = maxAbs;
        const scale = v => {
            if (v <= 0) {
                const t = (v - min) / (0 - min); // [0..1]
                const c1 = [33, 102, 172]; // blue
                const c2 = [247, 247, 247]; // grey
                const c = c1.map((c1i, i) => Math.round(c1i + t * (c2[i] - c1i)));
                return `rgb(${c[0]},${c[1]},${c[2]})`;
            } else {
                const t = v / max;
                const c1 = [247, 247, 247]; // grey
                const c2 = [178, 24, 43]; // red
                const c = c1.map((c1i, i) => Math.round(c1i + t * (c2[i] - c1i)));
                return `rgb(${c[0]},${c[1]},${c[2]})`;
            }
        };
        this.cy.nodes().forEach(n => {
            const name = n.data('id');
            if (name in dataMap) {
                const v = dataMap[name];
                const color = scale(v);
                n.style('background-color', color);
                n.style('color', Math.abs(v) > maxAbs * 0.6 ? '#fff' : '#000');
            } else {
                n.style('background-color', '#E0E0E0');
                n.style('color', '#000');
            }
        });
        const minEl = document.getElementById('legendMin');
        const maxEl = document.getElementById('legendMax');
        if (minEl && maxEl) {
            minEl.textContent = (-maxAbs).toFixed(2);
            maxEl.textContent = `+${maxAbs.toFixed(2)}`;
        }
    }

    resetColors() {
        this.cy.nodes().style({ 'background-color': '#fff', 'color': '#000' });
        const minEl = document.getElementById('legendMin');
        const maxEl = document.getElementById('legendMax');
        if (minEl && maxEl) {
            minEl.textContent = 'No data';
            maxEl.textContent = 'No data';
        }
    }
}



