// state.js
const createStateSet = () => {
    const internalSet = new Set();
    return {
        get: () => internalSet,
        add: (item) => internalSet.add(item),
        delete: (item) => internalSet.delete(item),
        clear: () => internalSet.clear(),
        has: (item) => internalSet.has(item),
    };
};

// Selection state sets
export const selectedGroups = createStateSet();
export const selectedRegulons = createStateSet();
export const selectedGenes = createStateSet();
export const selectedRegulations = createStateSet();
export const selectedExpressionClusters = createStateSet();
export const selectedMetaboliteCategories = createStateSet();
export const selectedMetabolites = createStateSet();

// Plot data and instances
export const plotData = {
    get: () => _plotData,
    set: (value) => { _plotData = value; }
};
let _plotData = [];

export const plotInstances = {
    get: () => _plotInstances,
    set: (value) => { _plotInstances = value; }
};
let _plotInstances = [];

// UI state
export const labelsVisible = {
    get: () => _labelsVisible,
    set: (value) => { _labelsVisible = value; }
};
let _labelsVisible = false;

export const legendVisible = {
    get: () => _legendVisible,
    set: (value) => { _legendVisible = value; }
};
let _legendVisible = true;

export const pinnedGene = {
    get: () => _pinnedGene,
    set: (value) => { _pinnedGene = value; }
};
let _pinnedGene = null;

export const isTooltipPinned = {
    get: () => _isTooltipPinned,
    set: (value) => { _isTooltipPinned = value; }
};
let _isTooltipPinned = false;

// Color management
export const currentColorMap = {
    get: () => _currentColorMap,
    set: (value) => { _currentColorMap = value; },
    clear: () => _currentColorMap.clear(),
    has: (key) => _currentColorMap.has(key),
    getValue: (key) => _currentColorMap.get(key),
    setValue: (key, value) => _currentColorMap.set(key, value)
};
let _currentColorMap = new Map();

// Color constants
export const FILTER_TYPE_COLORS = {
    group: '#1f77b4',       // Muted blue
    regulon: '#ff7f0e',    // Safety orange
    individual: '#2ca02c',  // Cooked asparagus green
    regulation: '#d62728',   // Brick red
    unknown: '#9467bd'       // Faded purple
};

export const INTERSECTION_COLOR = '#FF0000';
export const INTERSECTION_MARKER_SIZE = 12;

// Current selection result storage
export const currentSelectionResult = {
    get: () => _currentSelectionResult,
    set: (value) => { _currentSelectionResult = value; }
};
let _currentSelectionResult = null;

// === ADD: Pathways & Complexes State ===
let _pathwayInfo = {};
let _geneToPathwayMap = {};
let _selectedPathways = new Set();
let _coloringMode = 'category'; // 'category' | 'pathway'

let _proteinComplexMap = {};
let _selectedComplexes = new Set();

export const pathwayInfo = { get: () => _pathwayInfo, set: v => { _pathwayInfo = v || {}; } };
export const geneToPathwayMap = { get: () => _geneToPathwayMap, set: v => { _geneToPathwayMap = v || {}; } };
export const selectedPathways = {
    get: () => _selectedPathways,
    add: id => _selectedPathways.add(String(id)),
    delete: id => _selectedPathways.delete(String(id)),
    has: id => _selectedPathways.has(String(id)),
    clear: () => _selectedPathways.clear()
};
export const coloringMode = { get: () => _coloringMode, set: m => { _coloringMode = m === 'pathway' ? 'pathway' : 'category'; } };

export const proteinComplexMap = { get: () => _proteinComplexMap, set: v => { _proteinComplexMap = v || {}; } };
export const selectedComplexes = {
    get: () => _selectedComplexes,
    add: name => _selectedComplexes.add(String(name)),
    delete: name => _selectedComplexes.delete(String(name)),
    has: name => _selectedComplexes.has(String(name)),
    clear: () => _selectedComplexes.clear()
};
