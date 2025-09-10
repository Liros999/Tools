/**
 * JavaScript for the Configure Plots page
 * Handles dynamic addition of plot configurations
 */

document.addEventListener('DOMContentLoaded', function() {
    console.log('🔧 Configure plots page loaded');
    
    // Track plot count for each file
    const plotCounters = new Map();
    
    // Initialize plot counters
    document.querySelectorAll('.plot-configs-container').forEach(container => {
        const fileIndex = container.id.replace('plot-configs-', '');
        plotCounters.set(fileIndex, 1); // Start with 1 plot per file
    });
    
    // Handle "Add Another Plot" button clicks (robust to inner elements)
    document.addEventListener('click', function(event) {
        const addBtn = event.target.closest && event.target.closest('.add-plot-btn');
        if (addBtn) {
            // Prevent any default form behavior or bubbling that could refresh/reset the DOM
            event.preventDefault();
            event.stopPropagation();
            addNewPlotConfig(addBtn);
        }
    });
    
    function addNewPlotConfig(button) {
        console.log('➕ Adding new plot configuration');
        
        const containerSelector = button.dataset.container;
        const container = document.getElementById(containerSelector);
        
        if (!container) {
            console.error('❌ Container not found:', containerSelector);
            return;
        }
        
        const fileIndex = containerSelector.replace('plot-configs-', '');
        const currentPlotCount = plotCounters.get(fileIndex) || 1;
        const newPlotNumber = currentPlotCount + 1;
        
        // Get the file data to access columns
        const fileConfigSection = container.closest('.file-config-section');
        const fileName = fileConfigSection.querySelector('h4 .text-muted').textContent;
        
        // Get columns from the first plot config in this file
        const firstPlotConfig = container.querySelector('.plot-config-card');
        const columnOptions = Array.from(firstPlotConfig.querySelectorAll('select option')).map(option => ({
            value: option.value,
            text: option.textContent
        }));
        
        // Create new plot configuration card
        const newPlotCard = createPlotConfigCard(fileIndex, newPlotNumber, columnOptions);
        
        // Add the new card to the container
        container.appendChild(newPlotCard);
        
        // Update plot counter
        plotCounters.set(fileIndex, newPlotNumber);
        
        // Ensure card stays in view (no animation to avoid flicker)
        try { newPlotCard.scrollIntoView({ behavior: 'auto', block: 'center' }); } catch (_) {}
        
        console.log(`✅ Added plot ${newPlotNumber} for file ${fileIndex}`);
    }
    
    function createPlotConfigCard(fileIndex, plotNumber, columnOptions) {
        const cardDiv = document.createElement('div');
        cardDiv.className = 'plot-config-card card mb-3';
        
        // Create column options HTML
        const optionsHtml = columnOptions.map(option => 
            `<option value="${option.value}">${option.text}</option>`
        ).join('');
        
        cardDiv.innerHTML = `
            <div class="card-header d-flex justify-content-between align-items-center">
                <span>Plot ${plotNumber}</span>
                <button type="button" class="btn btn-sm btn-outline-danger remove-plot-btn" title="Remove this plot">
                    <i class="fas fa-trash"></i>
                </button>
            </div>
            <div class="card-body">
                <div class="row mb-3">
                    <div class="col-md-4">
                        <label for="gene_col_${fileIndex}_${plotNumber}" class="form-label">Gene Name Column</label>
                        <select class="form-select" id="gene_col_${fileIndex}_${plotNumber}" name="config_${fileIndex}_plot_${plotNumber}_gene_col" required>
                            ${optionsHtml}
                        </select>
                    </div>
                    <div class="col-md-4">
                        <label for="x_col_${fileIndex}_${plotNumber}" class="form-label">X-Axis (Log2 Fold Change)</label>
                        <select class="form-select" id="x_col_${fileIndex}_${plotNumber}" name="config_${fileIndex}_plot_${plotNumber}_x_col" required>
                            ${optionsHtml}
                        </select>
                    </div>
                    <div class="col-md-4">
                        <label for="y_col_${fileIndex}_${plotNumber}" class="form-label">Y-Axis (P-value)</label>
                        <select class="form-select" id="y_col_${fileIndex}_${plotNumber}" name="config_${fileIndex}_plot_${plotNumber}_y_col" required>
                            ${optionsHtml}
                        </select>
                    </div>
                </div>
                <div class="row">
                    <div class="col-md-6">
                        <label for="title_${fileIndex}_${plotNumber}" class="form-label">Plot Title</label>
                        <input type="text" class="form-control" id="title_${fileIndex}_${plotNumber}" name="config_${fileIndex}_plot_${plotNumber}_title" placeholder="e.g., My Volcano Plot ${plotNumber}">
                    </div>
                    <div class="col-md-3">
                        <label for="p_thresh_${fileIndex}_${plotNumber}" class="form-label">P-value Threshold</label>
                        <input type="number" step="any" class="form-control" id="p_thresh_${fileIndex}_${plotNumber}" name="config_${fileIndex}_plot_${plotNumber}_p_thresh" value="0.05">
                    </div>
                    <div class="col-md-3">
                        <label for="fc_thresh_${fileIndex}_${plotNumber}" class="form-label">Log2FC Threshold</label>
                        <input type="number" step="any" class="form-control" id="fc_thresh_${fileIndex}_${plotNumber}" name="config_${fileIndex}_plot_${plotNumber}_fc_thresh" value="1">
                    </div>
                </div>
            </div>
        `;
        
        // Add remove functionality
        const removeBtn = cardDiv.querySelector('.remove-plot-btn');
        removeBtn.addEventListener('click', function(e) {
            // Prevent any unintended form submission or outer handlers
            e.preventDefault();
            e.stopPropagation();
            removePlotConfig(cardDiv, fileIndex);
        });
        
        return cardDiv;
    }
    
    function removePlotConfig(cardElement, fileIndex) {
        console.log('🗑️ Removing plot configuration');
        
        // Immediate removal (no animation to avoid transient disappearance issues)
        cardElement.remove();
        
        // Update plot numbers for remaining cards in this file
        const container = document.getElementById(`plot-configs-${fileIndex}`);
        const plotCards = container.querySelectorAll('.plot-config-card');
        
        plotCards.forEach((card, index) => {
            const plotNumber = index + 1;
            const header = card.querySelector('.card-header span');
            header.textContent = `Plot ${plotNumber}`;
            
            // Update all IDs and names in this card
            updatePlotConfigReferences(card, fileIndex, plotNumber);
        });
        
        // Update counter
        plotCounters.set(fileIndex, plotCards.length);
        
        console.log(`✅ Plot removed, ${plotCards.length} plots remaining for file ${fileIndex}`);
    }
    
    function updatePlotConfigReferences(card, fileIndex, newPlotNumber) {
        // Update all form elements in the card to use the new plot number
        const elements = card.querySelectorAll('[id], [name], [for]');
        
        elements.forEach(element => {
            ['id', 'name', 'for'].forEach(attr => {
                if (element.hasAttribute(attr)) {
                    const value = element.getAttribute(attr);
                    // Replace the plot number in the attribute
                    const newValue = value.replace(
                        new RegExp(`_${fileIndex}_\\d+`), 
                        `_${fileIndex}_${newPlotNumber}`
                    );
                    element.setAttribute(attr, newValue);
                }
            });
        });
        
        // Update placeholder text
        const titleInput = card.querySelector('input[placeholder*="My Volcano Plot"]');
        if (titleInput) {
            titleInput.placeholder = `e.g., My Volcano Plot ${newPlotNumber}`;
        }
    }
});

// Remove animations/styles that could cause flicker/disappearance
