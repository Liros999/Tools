import sys
from pathlib import Path

# Ensure project root on sys.path when executed directly
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.api.metabolite_cache import (
    get_gene_to_metabolites,
    get_metabolite_to_genes,
    get_metabolite_categories,
    get_metabolite_to_categories,
)
from src.api.selection_api import (
    list_metabolite_categories,
    list_metabolites_by_category,
    get_genes_for_metabolite_selection,
)


def main():
    g2m = get_gene_to_metabolites() or {}
    m2g = get_metabolite_to_genes() or {}
    cats = get_metabolite_categories() or []
    m2c = get_metabolite_to_categories() or {}

    print("CACHE SUMMARY")
    print("  gene_to_metabolites:", len(g2m))
    print("  metabolite_to_genes:", len(m2g))
    print("  metabolite_categories:", len(cats))
    print("  metabolite_to_categories:", len(m2c))

    if not cats:
        print("No categories found; aborting sample selection test.")
        return

    # Find first non-empty category
    cid = None
    mets = []
    for c in cats:
        try:
            mc = list_metabolites_by_category(int(c['id']))
        except Exception:
            mc = []
        if mc:
            cid = int(c['id'])
            mets = mc
            break
    print(f"SELECTED CATEGORY id={cid} metabolites={len(mets)} (searched {cats.index(c)+1 if cid is not None else 0} categories)")

    if not mets:
        print("No category with metabolites found; aborting.")
        return

    payload = {
        'metabolite_ids': [int(m['id']) for m in mets[:3] if 'id' in m],
        'metabolite_category_ids': [],
        'interaction_types': [],
    }
    sel = get_genes_for_metabolite_selection(payload)
    print("SELECTION RESULT buckets=", len(sel))
    if sel:
        k = next(iter(sel))
        print("  first bucket:", k, "genes=", len(sel[k]))
        # print a few genes for visibility
        print("  sample genes:", (sel[k][:10] if sel[k] else []))


if __name__ == "__main__":
    main()


