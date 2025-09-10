/**
 * Comprehensive Gene Group Overlap Detection & Debugging Manager
 * Implements systematic debugging framework for tooltip and gene overlap issues
 */

class GeneOverlapDebugger {
    constructor(config = {}) {
        this.config = {
            debugMode: true,
            logLevel: 'verbose', // 'minimal', 'verbose', 'critical'
            enableVisualDebugging: true,
            ...config
        };
        this.debugMode = this.config.debugMode;
        this.overlapHistory = [];
        this.tooltipEvents = [];
        this.labelEvents = [];
        this.initializeDebugConsole();
    }
    
    // ===========================================================================
    // PHASE 1: DATA VALIDATION
    // ===========================================================================
    
    debugGeneData(genes, groups) {
        this.log("=== GENE DATA DEBUG ===", 'critical');
        
        // Check for duplicate gene IDs
        const geneIds = genes.map(g => g.id || g);
        const duplicates = geneIds.filter((id, index) => geneIds.indexOf(id) !== index);
        this.log("Duplicate Gene IDs:", duplicates);
        
        // Validate group memberships
        groups.forEach(group => {
            this.log(`Group ${group.name || group.id}:`, group.genes?.length || 0, "genes");
            this.log("Sample genes:", group.genes?.slice(0, 5));
        });
        
        // Check for genes in multiple groups
        const geneGroupMap = new Map();
        groups.forEach(group => {
            group.genes?.forEach(geneId => {
                if (!geneGroupMap.has(geneId)) {
                    geneGroupMap.set(geneId, []);
                }
                geneGroupMap.get(geneId).push(group.name || group.id);
            });
        });
        
        const overlappingGenes = Array.from(geneGroupMap.entries())
            .filter(([geneId, groupList]) => groupList.length > 1);
        
        this.log("Genes in multiple groups:", overlappingGenes);
        return { duplicates, geneGroupMap, overlappingGenes };
    }
    
    // ===========================================================================
    // PHASE 2: SELECTION LOGIC VERIFICATION
    // ===========================================================================
    
    debugGroupSelection(selectedGroups, allGroups) {
        this.log("=== GROUP SELECTION DEBUG ===", 'critical');
        this.log("Selected groups:", selectedGroups);
        
        if (selectedGroups.length < 2) {
            this.log("WARNING: Less than 2 groups selected - no overlaps possible", 'critical');
            return [];
        }
        
        // Find overlapping genes
        const selectedGroupData = allGroups.filter(g => selectedGroups.includes(g.id));
        const geneOccurrences = new Map();
        
        selectedGroupData.forEach(group => {
            this.log(`Processing group ${group.name || group.id} with ${group.genes?.length || 0} genes`);
            group.genes?.forEach(geneId => {
                geneOccurrences.set(geneId, (geneOccurrences.get(geneId) || 0) + 1);
            });
        });
        
        const overlappingGenes = Array.from(geneOccurrences.entries())
            .filter(([geneId, count]) => count > 1)
            .map(([geneId, count]) => ({ geneId, count }));
        
        this.log("Overlapping genes found:", overlappingGenes);
        this.overlapHistory.push({
            timestamp: Date.now(),
            selectedGroups: selectedGroups,
            overlappingGenes: overlappingGenes
        });
        return overlappingGenes;
    }
    
    // ===========================================================================
    // PHASE 3: VISUAL RENDERING DEBUG
    // ===========================================================================
    
    debugRedCircles(overlappingGenes, renderedElements) {
        this.log("=== RED CIRCLE RENDERING DEBUG ===", 'critical');
        
        overlappingGenes.forEach(({ geneId }) => {
            const geneElement = document.querySelector(`[data-gene-id="${geneId}"]`);
            const circleElement = document.querySelector(`[data-gene-id="${geneId}"] .red-circle`);
            
            this.log(`Gene ${geneId}:`);
            this.log("  - Element exists:", !!geneElement);
            this.log("  - Red circle exists:", !!circleElement);
            this.log("  - Element classes:", geneElement?.className);
            this.log("  - Circle styles:", circleElement?.style.cssText);
        });
    }
    
    // ===========================================================================
    // PHASE 4: LEGEND UPDATE DEBUG
    // ===========================================================================
    
    debugLegend(overlappingGenes, legendElement) {
        this.log("=== LEGEND DEBUG ===", 'critical');
        
        const expectedLegendItems = overlappingGenes.map(({ geneId, count }) => 
            `${geneId} (in ${count} groups)`
        );
        
        this.log("Expected legend items:", expectedLegendItems);
        
        const actualLegendItems = Array.from(legendElement?.children || [])
            .map(child => child.textContent);
        
        this.log("Actual legend items:", actualLegendItems);
        
        const missing = expectedLegendItems.filter(expected => 
            !actualLegendItems.some(actual => actual.includes(expected.split(' ')[0]))
        );
        
        this.log("Missing from legend:", missing);
    }
    
    // ===========================================================================
    // TOOLTIP DEBUGGING FRAMEWORK
    // ===========================================================================
    
    debugTooltipEvent(eventData) {
        this.log("=== TOOLTIP EVENT DEBUG ===", 'critical');
        
        const debugInfo = {
            timestamp: Date.now(),
            hasData: !!eventData,
            hasPoints: !!(eventData && eventData.points),
            pointsLength: eventData?.points?.length,
            eventType: eventData?.event?.type,
            clientX: eventData?.event?.clientX,
            clientY: eventData?.event?.clientY,
            firstPoint: eventData?.points?.[0] ? {
                x: eventData.points[0].x,
                y: eventData.points[0].y,
                text: eventData.points[0].text,
                traceType: eventData.points[0].trace?.type,
                traceName: eventData.points[0].trace?.name,
                curveNumber: eventData.points[0].curveNumber,
                pointNumber: eventData.points[0].pointNumber
            } : null
        };
        
        this.log("Tooltip event data:", debugInfo);
        this.tooltipEvents.push(debugInfo);
        
        // Analyze if this is a problematic event
        this.analyzeTooltipEvent(debugInfo);
        
        return debugInfo;
    }
    
    analyzeTooltipEvent(debugInfo) {
        const issues = [];
        
        if (!debugInfo.hasData) {
            issues.push("No event data");
        }
        
        if (!debugInfo.hasPoints) {
            issues.push("No points in event data");
        }
        
        if (debugInfo.pointsLength === 0) {
            issues.push("Empty points array");
        }
        
        if (debugInfo.firstPoint) {
            const point = debugInfo.firstPoint;
            
            if (!point.text) {
                issues.push("No gene name in point");
            }
            
            if (typeof point.x !== 'number' || typeof point.y !== 'number') {
                issues.push("Invalid coordinates");
            }
            
            if (point.traceName && this.isGreyTrace(point.traceName)) {
                issues.push("Hovering over grey/unselected trace");
            }
            
            if (point.traceType !== 'scatter') {
                issues.push("Not a scatter trace");
            }
        }
        
        if (issues.length > 0) {
            this.log("🚨 TOOLTIP ISSUES DETECTED:", issues, 'critical');
        }
    }
    
    isGreyTrace(traceName) {
        return traceName && (
            traceName.includes('Non-selected') || 
            traceName.includes('Genes') ||
            traceName === 'Genes' ||
            traceName === 'Non-selected genes' ||
            traceName === 'Non-selected Genes'
        );
    }
    
    // ===========================================================================
    // DRAGGABLE LABELS DEBUGGING FRAMEWORK
    // ===========================================================================
    
    debugDraggableLabels(plotId) {
        this.log("=== DRAGGABLE LABELS DEBUG ===", 'critical');
        
        const plotDiv = document.getElementById(plotId);
        if (!plotDiv) {
            this.log(`❌ Plot element ${plotId} not found for draggable labels setup`, 'critical');
            return false;
        }
        
        // Check if Plotly has initialized the element
        if (typeof plotDiv.on !== 'function') {
            this.log(`❌ Plotly element not initialized for ${plotId}`, 'critical');
            return false;
        }
        
        // Debug existing event handlers
        this.debugExistingEventHandlers(plotDiv);
        
        // Debug annotation state
        this.debugAnnotationState(plotDiv);
        
        // Debug draggable functionality
        this.debugDraggableFunctionality(plotDiv);
        
        return true;
    }
    
    debugExistingEventHandlers(plotDiv) {
        this.log("=== EXISTING EVENT HANDLERS DEBUG ===", 'critical');
        
        // Check for existing event listeners
        const hasRelayoutHandler = plotDiv._eventListeners && plotDiv._eventListeners.plotly_relayout;
        const hasClickHandler = plotDiv._eventListeners && plotDiv._eventListeners.plotly_click;
        const hasHoverHandler = plotDiv._eventListeners && plotDiv._eventListeners.plotly_hover;
        
        this.log("Existing event handlers:", {
            relayout: hasRelayoutHandler ? hasRelayoutHandler.length : 0,
            click: hasClickHandler ? hasClickHandler.length : 0,
            hover: hasHoverHandler ? hasHoverHandler.length : 0
        }, 'critical');
        
        // Check for conflicting handlers
        if (hasHoverHandler && hasHoverHandler.length > 1) {
            this.log("🚨 CONFLICTING HOVER HANDLERS DETECTED!", hasHoverHandler.length, 'critical');
        }
    }
    
    debugAnnotationState(plotDiv) {
        this.log("=== ANNOTATION STATE DEBUG ===", 'critical');
        
        const layout = plotDiv.layout || {};
        const annotations = layout.annotations || [];
        
        this.log("Current annotations:", {
            count: annotations.length,
            annotations: annotations.map((ann, index) => ({
                index: index,
                text: ann.text,
                x: ann.x,
                y: ann.y,
                showarrow: ann.showarrow,
                editable: ann.editable
            }))
        }, 'critical');
        
        // Check for editable annotations
        const editableAnnotations = annotations.filter(ann => ann.editable !== false);
        this.log("Editable annotations:", editableAnnotations.length, 'critical');
    }
    
    debugDraggableFunctionality(plotDiv) {
        this.log("=== DRAGGABLE FUNCTIONALITY DEBUG ===", 'critical');
        
        // Test if relayout events are working
        const testRelayout = () => {
            this.log("Testing relayout event functionality...", 'critical');
            
            // Simulate a small relayout to test if events fire
            const currentAnnotations = plotDiv.layout?.annotations || [];
            if (currentAnnotations.length > 0) {
                const testAnnotation = currentAnnotations[0];
                const newX = testAnnotation.x + 0.001; // Tiny movement
                
                Plotly.relayout(plotDiv, {
                    [`annotations[0].x`]: newX
                }).then(() => {
                    this.log("✅ Relayout test successful", 'critical');
                }).catch((error) => {
                    this.log("❌ Relayout test failed:", error, 'critical');
                });
            } else {
                this.log("⚠️ No annotations to test relayout", 'critical');
            }
        };
        
        // Run test after a short delay
        setTimeout(testRelayout, 1000);
    }
    
    debugLabelDragEvent(eventData, plotId) {
        this.log("=== LABEL DRAG EVENT DEBUG ===", 'critical');
        
        const debugInfo = {
            timestamp: Date.now(),
            plotId: plotId,
            eventType: 'label_drag',
            hasEventData: !!eventData,
            eventKeys: eventData ? Object.keys(eventData) : [],
            annotationChanges: this.extractAnnotationChanges(eventData)
        };
        
        this.log("Label drag event data:", debugInfo, 'critical');
        
        // Analyze drag event
        this.analyzeLabelDragEvent(debugInfo);
        
        return debugInfo;
    }
    
    extractAnnotationChanges(eventData) {
        if (!eventData) return [];
        
        const changes = [];
        Object.keys(eventData).forEach(key => {
            if (key.startsWith('annotations[') && (key.endsWith('].x') || key.endsWith('].y') || key.endsWith('].ax') || key.endsWith('].ay'))) {
                const match = key.match(/\[(\d+)\]/);
                if (match) {
                    changes.push({
                        annotationIndex: parseInt(match[1]),
                        property: key.split('.').pop(),
                        value: eventData[key]
                    });
                }
            }
        });
        
        return changes;
    }
    
    analyzeLabelDragEvent(debugInfo) {
        const issues = [];
        
        if (!debugInfo.hasEventData) {
            issues.push("No event data");
        }
        
        if (debugInfo.eventKeys.length === 0) {
            issues.push("No event keys");
        }
        
        if (debugInfo.annotationChanges.length === 0) {
            issues.push("No annotation changes detected");
        }
        
        if (issues.length > 0) {
            this.log("🚨 LABEL DRAG ISSUES DETECTED:", issues, 'critical');
        } else {
            this.log("✅ Label drag event appears normal", 'critical');
        }
    }
    
    debugLabelCreation(plotId, genesToLabel, annotations) {
        this.log("=== LABEL CREATION DEBUG ===", 'critical');
        
        const debugInfo = {
            timestamp: Date.now(),
            plotId: plotId,
            genesToLabelCount: genesToLabel.size,
            genesToLabel: Array.from(genesToLabel),
            annotationsCreated: annotations.length,
            annotations: annotations.map((ann, index) => ({
                index: index,
                text: ann.text,
                x: ann.x,
                y: ann.y,
                showarrow: ann.showarrow,
                editable: ann.editable,
                ax: ann.ax,
                ay: ann.ay
            }))
        };
        
        this.log("Label creation data:", debugInfo, 'critical');
        
        // Validate annotations
        this.validateCreatedAnnotations(debugInfo);
        
        return debugInfo;
    }
    
    validateCreatedAnnotations(debugInfo) {
        const issues = [];
        
        if (debugInfo.annotationsCreated === 0) {
            issues.push("No annotations created");
        }
        
        if (debugInfo.annotationsCreated !== debugInfo.genesToLabelCount) {
            issues.push(`Mismatch: ${debugInfo.genesToLabelCount} genes to label vs ${debugInfo.annotationsCreated} annotations created`);
        }
        
        debugInfo.annotations.forEach((ann, index) => {
            if (!ann.editable) {
                issues.push(`Annotation ${index} (${ann.text}) is not editable`);
            }
            
            if (ann.ax === undefined || ann.ay === undefined) {
                issues.push(`Annotation ${index} (${ann.text}) missing pixel offsets`);
            }
        });
        
        if (issues.length > 0) {
            this.log("🚨 LABEL CREATION ISSUES DETECTED:", issues, 'critical');
        } else {
            this.log("✅ Label creation appears successful", 'critical');
        }
    }
    
    debugLabelRemoval(plotId) {
        this.log("=== LABEL REMOVAL DEBUG ===", 'critical');
        
        const plotDiv = document.getElementById(plotId);
        if (!plotDiv) {
            this.log(`❌ Plot element ${plotId} not found for label removal`, 'critical');
            return false;
        }
        
        const beforeAnnotations = plotDiv.layout?.annotations || [];
        this.log(`Before removal: ${beforeAnnotations.length} annotations`, 'critical');
        
        // Attempt removal
        try {
            Plotly.relayout(plotId, { 
                annotations: [],
                editable: true
            }).then(() => {
                const afterAnnotations = plotDiv.layout?.annotations || [];
                this.log(`After removal: ${afterAnnotations.length} annotations`, 'critical');
                
                if (afterAnnotations.length === 0) {
                    this.log("✅ Label removal successful", 'critical');
                } else {
                    this.log("❌ Label removal failed - annotations still present", 'critical');
                }
            }).catch((error) => {
                this.log("❌ Label removal error:", error, 'critical');
            });
        } catch (error) {
            this.log("❌ Label removal exception:", error, 'critical');
        }
        
        return true;
    }
    
    // ===========================================================================
    // COMPREHENSIVE OVERLAP DETECTION
    // ===========================================================================
    
    detectOverlaps(selectedGroups, allGroups) {
        this.log("Starting overlap detection", 'critical');
        
        // Step 1: Validate inputs
        if (!this.validateInputs(selectedGroups, allGroups)) {
            return [];
        }
        
        // Step 2: Build gene-to-groups mapping
        const geneGroupMap = this.buildGeneGroupMapping(selectedGroups, allGroups);
        
        // Step 3: Find overlaps
        const overlaps = this.findOverlappingGenes(geneGroupMap);
        
        // Step 4: Debug output
        this.debugResults(overlaps);
        
        return overlaps;
    }
    
    validateInputs(selectedGroups, allGroups) {
        const isValid = selectedGroups?.length >= 2 && Array.isArray(allGroups);
        this.log("Input validation:", isValid ? "PASS" : "FAIL", 'critical');
        return isValid;
    }
    
    buildGeneGroupMapping(selectedGroups, allGroups) {
        const mapping = new Map();
        const selectedGroupData = allGroups.filter(g => selectedGroups.includes(g.id));
        
        selectedGroupData.forEach(group => {
            group.genes?.forEach(geneId => {
                if (!mapping.has(geneId)) {
                    mapping.set(geneId, []);
                }
                mapping.get(geneId).push(group);
            });
        });
        
        this.log("Gene-group mapping built:", mapping.size, "unique genes");
        return mapping;
    }
    
    findOverlappingGenes(geneGroupMap) {
        const overlaps = Array.from(geneGroupMap.entries())
            .filter(([geneId, groups]) => groups.length > 1)
            .map(([geneId, groups]) => ({
                geneId,
                groups: groups.map(g => g.name || g.id),
                count: groups.length
            }));
        
        this.log("Overlapping genes found:", overlaps.length);
        return overlaps;
    }
    
    debugResults(overlaps) {
        if (overlaps.length > 0) {
            console.table(overlaps);
        } else {
            this.log("No overlapping genes found");
        }
    }
    
    // ===========================================================================
    // DEBUG CONSOLE MANAGEMENT
    // ===========================================================================
    
    initializeDebugConsole() {
        if (this.config.enableVisualDebugging) {
            this.createDebugPanel();
        }
        
        // Make debugger globally accessible
        window.geneOverlapDebugger = this;
        
        this.log("Gene Overlap Debugger initialized", 'critical');
    }
    
    createDebugPanel() {
        const panel = document.createElement('div');
        panel.id = 'debug-panel';
        panel.style.cssText = `
            position: fixed;
            top: 10px;
            right: 10px;
            width: 300px;
            max-height: 400px;
            background: rgba(0,0,0,0.9);
            color: white;
            padding: 10px;
            border-radius: 5px;
            font-family: monospace;
            font-size: 12px;
            z-index: 10000;
            overflow-y: auto;
            display: none;
        `;
        
        panel.innerHTML = `
            <div style="margin-bottom: 10px;">
                <strong>Gene Overlap Debugger</strong>
                <button onclick="this.parentElement.parentElement.style.display='none'" style="float: right; background: red; color: white; border: none; padding: 2px 5px;">X</button>
            </div>
            <div id="debug-content"></div>
        `;
        
        document.body.appendChild(panel);
        
        // Add toggle button
        const toggleBtn = document.createElement('button');
        toggleBtn.textContent = '🐛 Debug';
        toggleBtn.style.cssText = `
            position: fixed;
            top: 10px;
            right: 10px;
            z-index: 10001;
            background: #007bff;
            color: white;
            border: none;
            padding: 5px 10px;
            border-radius: 3px;
            cursor: pointer;
        `;
        toggleBtn.onclick = () => {
            panel.style.display = panel.style.display === 'none' ? 'block' : 'none';
        };
        document.body.appendChild(toggleBtn);
    }
    
    updateDebugPanel(content) {
        const debugContent = document.getElementById('debug-content');
        if (debugContent) {
            debugContent.innerHTML = content;
        }
    }
    
    // ===========================================================================
    // LOGGING SYSTEM
    // ===========================================================================
    
    log(...args) {
        if (!this.debugMode) return;
        
        const level = args[args.length - 1];
        const isCritical = level === 'critical';
        const isVerbose = this.config.logLevel === 'verbose';
        
        if (isCritical || isVerbose || this.config.logLevel === 'minimal') {
            const prefix = isCritical ? '🚨' : '🔍';
            console.log(prefix, "[GENE OVERLAP DEBUG]", ...args.slice(0, -1));
        }
    }
    
    // ===========================================================================
    // UTILITY METHODS
    // ===========================================================================
    
    getDebugSummary() {
        return {
            overlapHistory: this.overlapHistory,
            tooltipEvents: this.tooltipEvents,
            labelEvents: this.labelEvents,
            totalOverlapsDetected: this.overlapHistory.reduce((sum, entry) => sum + entry.overlappingGenes.length, 0),
            totalTooltipEvents: this.tooltipEvents.length,
            totalLabelEvents: this.labelEvents.length,
            problematicTooltipEvents: this.tooltipEvents.filter(event => !event.hasData || !event.hasPoints || event.pointsLength === 0).length
        };
    }
    
    clearHistory() {
        this.overlapHistory = [];
        this.tooltipEvents = [];
        this.labelEvents = [];
        this.log("Debug history cleared");
    }
    
    exportDebugData() {
        const data = {
            timestamp: new Date().toISOString(),
            summary: this.getDebugSummary(),
            overlapHistory: this.overlapHistory,
            tooltipEvents: this.tooltipEvents,
            labelEvents: this.labelEvents
        };
        
        const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `gene-overlap-debug-${Date.now()}.json`;
        a.click();
        URL.revokeObjectURL(url);
        
        this.log("Debug data exported");
    }
}

// Export for use in other modules
export default GeneOverlapDebugger;
