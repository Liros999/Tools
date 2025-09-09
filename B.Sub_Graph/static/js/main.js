(function() {
    // Helper functions for UI updates
    function updateStatus(message, isError = false) {
        const statusEl = document.getElementById('status');
        const messageDiv = `<div class="${isError ? 'error' : 'ok'}">${escapeHtml(message)}</div>`;
        statusEl.innerHTML = messageDiv;
    }

    function escapeHtml(str) {
        return String(str).replace(/[&<"']/g, c => ({
            '&': '&amp;',
            '<': '&lt;',
            '>': '&gt;',
            '"': '&quot;',
            "'": '&#39;'
        } [c]));
    }

    // API interaction
    async function fetchFullGraph() {
        const response = await fetch('/api/full_graph');
        const json = await response.json();
        if (!response.ok || !json.is_success) {
            throw (json.message || 'API error');
        }
        return json.data;
    }

    async function uploadAndProcess(file, geneCol, log2fcCol, pvalCol) {
        const formData = new FormData();
        formData.append('file', file);
        if (geneCol) formData.append('gene_col', geneCol);
        if (log2fcCol) formData.append('log2fc_col', log2fcCol);
        if (pvalCol) formData.append('pval_col', pvalCol);

        const response = await fetch('/api/process_upload', {
            method: 'POST',
            body: formData
        });
        const json = await response.json();
        if (!response.ok || !json.is_success) {
            throw (json.message || 'Upload error');
        }
        return json.data;
    }

    // Cytoscape graph handling
    let cy = null;

    function initializeCytoscape() {
        cy = window.cytoscape({
            container: document.getElementById('cy'),
            elements: [],
            layout: {
                name: 'dagre',
                rankSep: 60,
                nodeSep: 40,
                edgeSep: 20
            },
            style: [{
                selector: 'node',
                style: {
                    'label': 'data(label)',
                    'font-size': 10,
                    'text-valign': 'center',
                    'color': '#000',
                    'background-color': '#fff',
                    'border-width': 1,
                    'border-color': '#777'
                }
            }, {
                selector: 'edge',
                style: {
                    'width': 1.5,
                    'line-color': '#999',
                    'curve-style': 'bezier',
                    'target-arrow-shape': 'triangle',
                    'target-arrow-color': '#999'
                }
            }]
        });
    }

    function drawGraph(molecules, interactions) {
        const proteinNames = new Set(molecules.filter(m => m.type === 'Protein').map(m => m.name));
        const idToNameMap = new Map(molecules.map(m => [m.id, m.name]));

        const nodes = [...proteinNames].map(name => ({
            data: {
                id: name,
                label: name
            }
        }));

        const edges = [];
        for (const it of interactions) {
            if (Array.isArray(it.molecule_ids) && it.molecule_ids.length === 2) {
                const sourceName = idToNameMap.get(it.molecule_ids[0]);
                const targetName = idToNameMap.get(it.molecule_ids[1]);
                if (sourceName && targetName && proteinNames.has(sourceName) && proteinNames.has(targetName)) {
                    edges.push({
                        data: {
                            id: `${sourceName}__${targetName}__${it.id}`,
                            source: sourceName,
                            target: targetName
                        }
                    });
                }
            }
        }

        cy.elements().remove();
        cy.add(nodes).add(edges);
        cy.layout({
            name: 'dagre',
            rankSep: 60,
            nodeSep: 40,
            edgeSep: 20
        }).run();
    }

    function applyDataMapToGraph(dataMap) {
        const values = Object.values(dataMap).filter(v => v !== null);
        if (values.length === 0) return;

        const min = Math.min(...values);
        const max = Math.max(...values);
        const absMax = Math.max(Math.abs(min), Math.abs(max));

        const blue = [33, 102, 172];
        const white = [247, 247, 247];
        const red = [178, 24, 43];

        cy.nodes().forEach(node => {
            const name = node.data('id');
            const val = dataMap[name];

            if (val === undefined || val === null) {
                node.style({
                    'background-color': '#fff',
                    'border-color': '#777'
                });
            } else {
                const t = val / absMax;
                const color = t < 0 ? interpolateColor(white, blue, -t) : interpolateColor(white, red, t);
                node.style({
                    'background-color': color,
                    'border-color': '#000'
                });
            }
        });

        updateLegend(absMax, white, blue, red);
    }

    function interpolateColor(color1, color2, factor) {
        const result = color1.slice();
        for (let i = 0; i < 3; i++) {
            result[i] = Math.round(result[i] + factor * (color2[i] - result[i]));
        }
        return `rgb(${result[0]}, ${result[1]}, ${result[2]})`;
    }

    function updateLegend(absMax, white, blue, red) {
        document.getElementById('legendMin').textContent = -absMax.toFixed(2);
        document.getElementById('legendMid').textContent = '0';
        document.getElementById('legendMax').textContent = absMax.toFixed(2);
        document.getElementById('legendGradient').style.background =
            `linear-gradient(to right, ${interpolateColor(white, blue, 1)}, ${interpolateColor(white, red, 1)})`;
    }


    // File parsing
    function parseFileHeaders(file) {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            reader.onload = () => {
                try {
                    const text = reader.result;
                    const firstLine = text.slice(0, text.indexOf('\n'));
                    const delimiter = /\t/.test(firstLine) ? '\t' : ',';
                    const headers = firstLine.split(delimiter).map(h => h.trim().replace(/^"|"$/g, ''));
                    resolve(headers);
                } catch (e) {
                    reject(new Error('Could not parse file headers.'));
                }
            };
            reader.onerror = () => reject(new Error('Error reading file.'));
            reader.readAsText(file);
        });
    }


    // Event Listeners and DOM Initialization
    document.addEventListener('DOMContentLoaded', () => {
        initializeCytoscape();

        const fileInput = document.getElementById('fileInput');
        const columnsBlock = document.getElementById('columnsBlock');
        const processBtn = document.getElementById('processBtn');
        const geneColSelect = document.getElementById('geneCol');
        const log2fcColSelect = document.getElementById('log2fcCol');
        const pvalColSelect = document.getElementById('pvalCol');

        // Load initial full graph
        fetchFullGraph().then(graph => {
            drawGraph(graph.molecules, graph.interactions);
            updateStatus('Full graph loaded. Upload a file to see your data.');
        }).catch(error => {
            updateStatus(error.message || 'Could not load full graph.', true);
            console.error(error);
        });

        // File input change handler
        fileInput.addEventListener('change', async () => {
            columnsBlock.style.display = 'none';
            if (!fileInput.files || fileInput.files.length === 0) return;

            updateStatus('Reading file headers...');
            try {
                const headers = await parseFileHeaders(fileInput.files[0]);
                [geneColSelect, log2fcColSelect, pvalColSelect].forEach(select => select.innerHTML = '');

                headers.forEach(header => {
                    const option = document.createElement('option');
                    option.value = header;
                    option.textContent = header;
                    geneColSelect.appendChild(option.cloneNode(true));
                    log2fcColSelect.appendChild(option.cloneNode(true));
                    pvalColSelect.appendChild(option.cloneNode(true));
                });

                // Auto-select columns based on common names
                geneColSelect.value = headers.find(h => /gene|name|id/i.test(h)) || headers[0];
                log2fcColSelect.value = headers.find(h => /log2fc|log2/i.test(h)) || (headers.length > 1 ? headers[1] : headers[0]);
                pvalColSelect.value = headers.find(h => /pval|p-value|padj/i.test(h)) || (headers.length > 2 ? headers[2] : headers[0]);

                columnsBlock.style.display = 'flex';
                updateStatus('Ready. Select columns and click "Process File".');
            } catch (error) {
                updateStatus(error.message, true);
            }
        });

        // Process button click handler
        processBtn.addEventListener('click', async () => {
            if (!fileInput.files || fileInput.files.length === 0) {
                updateStatus('No file selected.', true);
                return;
            }
            updateStatus('Processing...');
            try {
                const data = await uploadAndProcess(fileInput.files[0], geneColSelect.value, log2fcColSelect.value, pvalColSelect.value);
                applyDataMapToGraph(data.dataMap);
                updateStatus(`Processed ${Object.keys(data.dataMap).length} genes with data.`);
            } catch (error) {
                updateStatus(error.message || 'Processing failed.', true);
                console.error(error);
            }
        });
    });
})();