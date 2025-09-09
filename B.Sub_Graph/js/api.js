export async function fetchFullGraph() {
    const r = await fetch('http://localhost:5001/api/full_graph');
    const j = await r.json();
    if (!r.ok || !j.is_success) throw j.message || j;
    return j.data;
}

export async function uploadAndProcess(file, geneCol, log2fcCol, pvalCol) {
    const form = new FormData();
    form.append('file', file);
    if (geneCol) form.append('gene_col', geneCol);
    if (log2fcCol) form.append('log2fc_col', log2fcCol);
    if (pvalCol) form.append('pval_col', pvalCol);

    const r = await fetch('http://localhost:5001/api/process_upload', { method: 'POST', body: form });
    const j = await r.json();
    if (!r.ok || !j.is_success) throw j.message || j;
    return j.data;
}



