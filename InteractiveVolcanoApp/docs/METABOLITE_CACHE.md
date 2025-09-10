### Metabolite Cache Builder

This prebuild step fetches real SubtiWiki data and writes deterministic caches that the app can load at runtime (no placeholders, no fallbacks).

Outputs (written to `src/data/cache/`):
- `gene_to_metabolites.json`: `{ gene_name: [{metabolite_id, metabolite_name, type, transport_type?}], ... }`
- `metabolite_to_genes.json`: `{ metabolite_id: [gene_name, ...], ... }`
- `metabolite_categories.json`: list of category nodes `{id, name, dot_notation, parent_id, children}`
- `metabolite_to_categories.json`: `{ metabolite_id: [{id, name, dot_notation, parent_id}], ... }`

Data sources (see `Tools/Subtiwiki_API.txt`):
- `/interaction/metabolite-protein`
- `/protein/{id}` → `gene_id`
- `/gene/{gene_id}?representation=minimal` → `name`
- `/metabolite-category/` and `/metabolite-category/{id}`

Run (PowerShell, Windows):
```powershell
python Tools/InteractiveVolcanoApp/docs/development/build_metabolite_mappings.py
```

Notes:
- Parallelized resolution of protein→gene names; polite short sleeps on category hydration.
- All responses are accessed via the standard `data` envelope when present.
- If SubtiWiki is temporarily unavailable, rerun the builder to refresh caches.


