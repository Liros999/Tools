export function setStatusOk(msg) {
    const el = document.getElementById('status');
    el.innerHTML = `<div class="ok">${escapeHtml(msg)}</div>`;
}

export function setStatusError(msg) {
    const el = document.getElementById('status');
    el.innerHTML = `<div class="error">${escapeHtml(msg)}</div>`;
}

export function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;','\'':'&#39;'}[c]));
}

export function populateColumnSelectors(headers) {
    const geneSel = document.getElementById('geneCol');
    const l2Sel = document.getElementById('log2fcCol');
    const pvSel = document.getElementById('pvalCol');
    for (const h of headers) {
        const o1 = document.createElement('option'); o1.value = h; o1.textContent = h; geneSel.appendChild(o1);
        const o2 = document.createElement('option'); o2.value = h; o2.textContent = h; l2Sel.appendChild(o2);
        const o3 = document.createElement('option'); o3.value = h; o3.textContent = h; pvSel.appendChild(o3);
    }
    // Best-guess selection (no hardcoded thresholds)
    const lower = headers.map(h => h.toLowerCase());
    geneSel.value = headers[lower.findIndex(h => h.includes('gene') || h.includes('name') || h.includes('id'))] || headers[0];
    l2Sel.value = headers[lower.findIndex(h => h.includes('log2fc') || h.includes('log2') || h.includes('fold'))] || headers[1] || headers[0];
    const pidx = lower.findIndex(h => h.includes('pval') || h.includes('p_value') || h.includes('padj') || h.includes('fdr'));
    if (pidx >= 0) pvSel.value = headers[pidx];
}



