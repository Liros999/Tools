/**
 * User Interface and DOM Manipulation Module
 */

/**
 * Sets up the file input to parse the CSV header and enable the UI.
 * @param {Function} onFileReady - Callback executed with the selected file.
 */
export function setupFileInput(onFileReady) {
    const fileInput = document.getElementById('fileInput');
    fileInput.addEventListener('change', async (event) => {
        const file = event.target.files[0];
        if (!file) return;

        try {
            const headers = await parseCSVHeader(file);
            populateColumnDropdowns(headers);
            document.getElementById('column-selection').style.display = 'block';
            document.getElementById('filter-options').style.display = 'block';
            if (onFileReady) onFileReady(file);
        } catch (error) {
            showErrorMessage(`Error reading file: ${error.message}`);
        }
    });
}

function parseCSVHeader(file) {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (e) => {
            const firstLine = e.target.result.slice(0, e.target.result.indexOf('\n')).trim();
            if (!firstLine) return reject(new Error('Cannot find header row.'));
            const headers = firstLine.split(/[\t,]/).map(h => h.trim().replace(/"/g, ''));
            resolve(headers);
        };
        reader.onerror = () => reject(new Error('Failed to read file.'));
        reader.readAsText(file);
    });
}

function populateColumnDropdowns(headers) {
    const selectors = {
        gene: document.getElementById('gene-column'),
        log2fc: document.getElementById('log2fc-column'),
        pval: document.getElementById('pval-column'),
    };
    Object.values(selectors).forEach(sel => sel.innerHTML = '');
    headers.forEach(header => {
        Object.values(selectors).forEach(sel => sel.add(new Option(header, header)));
    });
    selectors.gene.value = findBestMatch(headers, ['gene', 'name', 'id']);
    selectors.log2fc.value = findBestMatch(headers, ['log2fc', 'log2']);
    selectors.pval.value = findBestMatch(headers, ['pval', 'padj']);
}

function findBestMatch(headers, keywords) {
    for (const keyword of keywords) {
        const match = headers.find(h => h.toLowerCase().includes(keyword));
        if (match) return match;
    }
    return headers[0] || null;
}

/**
 * Dynamically generates and updates the color legend.
 * @param {d3.ScaleLinear<number, string>} scale - The D3 color scale.
 */
export function updateLegend(scale) {
    const legendContainer = d3.select('#color-legend');
    legendContainer.html(''); // Clear previous legend

    const legendWidth = 250, legendHeight = 20;

    const canvas = legendContainer.append('canvas')
        .attr('width', legendWidth)
        .attr('height', 1)
        .style('width', `${legendWidth}px`)
        .style('height', `${legendHeight}px`)
        .style('border', '1px solid #ccc')
        .style('border-radius', '4px')
        .node();

    const context = canvas.getContext('2d');
    const image = context.createImageData(legendWidth, 1);
    const domain = scale.domain();
    const range = d3.range(domain[0], domain[domain.length - 1], (domain[domain.length - 1] - domain[0]) / legendWidth);

    for (let i = 0; i < legendWidth; ++i) {
        const c = d3.rgb(scale(range[i]));
        image.data[4 * i] = c.r;
        image.data[4 * i + 1] = c.g;
        image.data[4 * i + 2] = c.b;
        image.data[4 * i + 3] = 255;
    }
    context.putImageData(image, 0, 0);

    // Add axis labels
    const legendAxis = d3.scaleLinear()
        .domain(scale.domain())
        .range([0, legendWidth]);

    const svg = legendContainer.append('svg')
        .attr('width', legendWidth + 20)
        .attr('height', 30)
        .style('margin-left', '-10px');

    const axis = d3.axisBottom(legendAxis)
        .ticks(5)
        .tickFormat(d3.format(".1f"));

    svg.append('g')
        .attr('transform', `translate(10, 0)`)
        .call(axis);
}

/**
 * Creates the HTML content for a gene tooltip.
 * @param {object} context - The gene context data from the backend.
 * @returns {string} HTML string.
 */
export function createTooltipHtml(context) {
    let html = `<div class=\"tooltip-content\"><strong>${context.name}</strong>`;
    if (context.function) {
        html += `<br><em>${context.function}</em>`;
    }
    html += '<hr>';
    if (context.synonyms && context.synonyms.length > 0) {
        html += `<strong>Synonyms:</strong> ${context.synonyms.join(', ')}<br>`;
    }
    if (context.categories && context.categories.length > 0) {
        html += `<strong>Categories:</strong> ${context.categories.map(c => c.name).join(', ')}`;
    }
    // Add more fields as needed, e.g., outlinks
    html += `</div>`;
    return html;
}

// --- Messaging --- //
export function showLoadingMessage(message) {
    const loader = document.getElementById('loader');
    if (loader) {
        loader.textContent = message;
        loader.style.display = 'flex';
    }
}

export function hideLoadingMessage() {
    const loader = document.getElementById('loader');
    if (loader) {
        loader.style.display = 'none';
    }
}
