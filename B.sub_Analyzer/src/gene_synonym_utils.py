import os
import json
from tqdm import tqdm

def build_synonym_index_from_cache(cache_dir, index_file):
    print(f"Building synonym index from cache: {cache_dir}")
    synonym_to_id = {}
    synonym_to_canonical = {}
    name_to_id = {}
    for fname in tqdm(os.listdir(cache_dir), desc="Cache files", unit="file"):
        if not fname.endswith('.json'):
            continue
        with open(os.path.join(cache_dir, fname), 'r') as f:
            try:
                data = json.load(f)
                canonical = data.get('gene_name')
                gene_id = data.get('gene_id')
                if canonical and gene_id:
                    name_to_id[canonical.lower()] = gene_id
                    synonym_to_id[canonical.lower()] = gene_id
                    synonym_to_canonical[canonical.lower()] = canonical
                for syn in data.get('synonyms', []):
                    if syn:
                        synonym_to_id[syn.lower()] = gene_id
                        synonym_to_canonical[syn.lower()] = canonical
            except Exception as e:
                continue
    index = {
        'name_to_id': name_to_id,
        'synonym_to_id': synonym_to_id,
        'synonym_to_canonical': synonym_to_canonical
    }
    with open(index_file, 'w') as f:
        json.dump(index, f)
    print(f"Synonym index saved to {index_file}")
    return name_to_id, synonym_to_id, synonym_to_canonical

def load_synonym_index(cache_dir):
    index_file = os.path.join(cache_dir, 'gene_synonym_index.json')
    if os.path.exists(index_file):
        with open(index_file, 'r') as f:
            index = json.load(f)
        name_to_id = index['name_to_id']
        synonym_to_id = index['synonym_to_id']
        synonym_to_canonical = index['synonym_to_canonical']
    else:
        name_to_id, synonym_to_id, synonym_to_canonical = build_synonym_index_from_cache(cache_dir, index_file)
    return name_to_id, synonym_to_id, synonym_to_canonical

def resolve_gene_name(name, name_to_id, synonym_to_id, synonym_to_canonical):
    key = name.lower()
    if key in name_to_id:
        return name_to_id[key], key, None
    elif key in synonym_to_id:
        canonical = synonym_to_canonical[key]
        return synonym_to_id[key], canonical, key
    else:
        return None, None, None 