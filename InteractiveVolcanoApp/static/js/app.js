import * as api from './modules/apiService.js';
import * as ui from './modules/uiManager.js';
import { pathwayInfo, geneToPathwayMap } from './modules/state.js';
import * as plot from './modules/plotManager.js';
import * as state from './modules/state.js';

class VolcanoPlotApp {
  constructor() {
    this.init();
  }

  async init() {
    console.log('Initializing Modular Volcano Plot App');
    // Use a more robust DOM ready check
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', () => this.start());
    } else if (document.readyState === 'interactive') {
      // DOM is parsed but subresources may still be loading
      setTimeout(() => this.start(), 100);
    } else {
      this.start();
    }
  }

  async start() {
    try {
      // Wait for DOM to be fully loaded
      if (document.readyState === 'loading') {
        await new Promise(resolve => {
          document.addEventListener('DOMContentLoaded', resolve);
        });
      } else if (document.readyState === 'interactive') {
        // Give a bit more time for DOM to be fully ready
        await new Promise(resolve => setTimeout(resolve, 100));
      }
      
      // Verify key elements exist before proceeding
      const keyElements = [
        'gene-search',
        'group-selection-tree',
        'apply-selection',
        'clear-selection'
      ];
      
      const missingElements = keyElements.filter(id => !document.getElementById(id));
      if (missingElements.length > 0) {
        console.warn('Some elements not found:', missingElements);
        console.log('Available elements with similar IDs:');
        keyElements.forEach(id => {
          const element = document.getElementById(id);
          if (!element) {
            // Look for similar elements
            const similarElements = Array.from(document.querySelectorAll(`[id*="${id.split('-')[0]}"]`));
            if (similarElements.length > 0) {
              console.log(`  Similar to ${id}:`, similarElements.map(el => el.id));
            }
          }
        });
        
        // Check if Bootstrap is loaded and tabs are initialized
        console.log('Checking Bootstrap and tab initialization...');
        if (typeof bootstrap !== 'undefined') {
          console.log('Bootstrap is loaded');
          // Try to initialize tabs manually if needed
          const tabElements = document.querySelectorAll('[data-bs-toggle="pill"]');
          console.log('Found tab elements:', tabElements.length);
        } else {
          console.warn('Bootstrap not loaded');
        }
        
        // Wait a bit more and try again
        await new Promise(resolve => setTimeout(resolve, 1000));
        
        // Check again after waiting
        const stillMissing = keyElements.filter(id => !document.getElementById(id));
        if (stillMissing.length > 0) {
          console.warn('Elements still missing after retry:', stillMissing);
        }
      }
      
      // Setup event handlers first
      ui.setupEventHandlers();
      ui.setupGeneSearch();
      
      // Load and render hierarchical groups tree
      ui.renderHierarchicalTree();
      // Initialize coloring controls
      ui.initializeColoringUI();

      // Fetch pathway and complex LIST endpoints per final design; show informative UI on 404
      try {
        const t0 = performance.now();
        const r = await fetch('/api/pathways/list');
        if (r.ok) {
          const ttfbMs = performance.now() - t0; // rough client-side TTFB estimate
          const d = await r.json();
          if (d.success && Array.isArray(d.pathways)) {
            // Convert list to minimal info map for UI builder compatibility
            const infoMap = {};
            d.pathways.forEach(p => { if (p && p.id != null) infoMap[String(p.id)] = { name: p.name || `Pathway ${p.id}` }; });
            pathwayInfo.set(infoMap);
            ui.buildPathwayFilter();
            console.log('[net] /api/pathways/list', {
              client_ttfb_ms: Math.round(ttfbMs),
              count: d.pathways.length
            });
          }
        } else if (r.status === 404) {
          const pc = document.getElementById('pathway-filter-container');
          if (pc) { pc.classList.remove('d-none'); pc.innerHTML = '<p class="text-muted small"><em>Pathway coloring is not available. Cache data missing on server.</em></p>'; }
        }
      } catch (e) {
        const pc = document.getElementById('pathway-filter-container');
        if (pc) { pc.classList.remove('d-none'); pc.innerHTML = '<p class="text-muted small"><em>Pathway coloring failed to load.</em></p>'; }
      }

      // Complexes feature removed per updated design
      
      if (window.plotData && Array.isArray(window.plotData)) {
        state.plotData.set(window.plotData);
        plot.initializePlots();
      }
      
      console.log('App initialized successfully');
    } catch (error) {
      console.error('Failed to initialize modular app:', error);
    }
  }

  // Make UI functions globally available for onclick handlers
  addGeneToSelection = ui.addGeneToSelection;
  removeGeneFromSelection = ui.removeGeneFromSelection;
}

const app = new VolcanoPlotApp();
window.app = app;

// Global functions for HTML onclick handlers
window.closeSidebarTooltip = function() {
    const sidebar = document.getElementById('gene-sidebar-tooltip');
    if (sidebar) {
        sidebar.style.display = 'none';
    }
    state.pinnedGene.set(null);
    state.isTooltipPinned.set(false);
};
