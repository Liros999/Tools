import { fetchFullGraph, uploadAndProcess } from './api.js';
import { CytoscapeManager } from './graph.js';
import { setStatusOk, setStatusError, populateColumnSelectors } from './ui.js';

let cyManager;
let lastHeaders = null;
let lastFile = null;

async function init() {
    cyManager = new CytoscapeManager('cy');
    setStatusOk('Fetching SubtiWiki interaction graph...');
    try {
        const { molecules, interactions } = await fetchFullGraph();
        cyManager.draw(molecules, interactions);
        setStatusOk(`Loaded ${molecules.length} molecules and ${interactions.length} interactions`);
    } catch (e) {
        setStatusError(describeError(e));
    }

    document.getElementById('fileInput').addEventListener('change', onFileChosen);
    document.getElementById('processBtn').addEventListener('click', onProcess);
}

function describeError(e) {
    if (!e) return 'Unknown error';
    if (e.message && e.type) return `${e.type}: ${e.message}`;
    if (e.type && e.message) return `${e.type}: ${e.message}`;
    return typeof e === 'string' ? e : (e.message || JSON.stringify(e));
}

async function onFileChosen(ev) {
    const f = ev.target.files[0];
    if (!f) return;
    lastFile = f;
    // Read headers client-side for selector population
    try {
        const text = await f.text();
        const firstLine = text.split(/\r?\n/)[0] || '';
        const delim = /\t/.test(firstLine) ? '\t' : (firstLine.includes(';') ? ';' : ',');
        lastHeaders = firstLine.split(new RegExp(delim)).map(h => h.trim().replace(/^"|"$/g, ''));
        populateColumnSelectors(lastHeaders);
        document.getElementById('columnsBlock').style.display = '';
        setStatusOk('Select columns and click Process File');
    } catch (e) {
        setStatusError('Failed reading file headers. Please ensure CSV/TSV format.');
    }
}

async function onProcess() {
    if (!lastFile) {
        setStatusError('No file selected');
        return;
    }
    const geneCol = document.getElementById('geneCol').value;
    const log2fcCol = document.getElementById('log2fcCol').value;
    const pvalCol = document.getElementById('pvalCol').value;
    setStatusOk('Processing file and resolving genes via SubtiWiki...');
    try {
        const data = await uploadAndProcess(lastFile, geneCol, log2fcCol, pvalCol);
        // Apply colors dynamically from actual data
        cyManager.applyColorsFromMap(data.dataMap);
        if (data.notFoundGenes && data.notFoundGenes.length) {
            setStatusOk(`Colored ${Object.keys(data.dataMap).length} genes. ${data.notFoundGenes.length} names not found (see console).`);
            console.log('Not found genes:', data.notFoundGenes);
        } else {
            setStatusOk(`Colored ${Object.keys(data.dataMap).length} genes.`);
        }
    } catch (e) {
        setStatusError(describeError(e));
    }
}

document.addEventListener('DOMContentLoaded', init);



