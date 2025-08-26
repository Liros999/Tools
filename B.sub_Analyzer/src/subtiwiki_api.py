import requests

BASE_URL = "https://subtiwiki.uni-goettingen.de/v5/api"

def fetch_gene_data(gene_name):
    """
    Fetch all available data for a gene from SubtiWiki.
    Returns a dictionary with gene info, or None if not found.
    """
    url = f"{BASE_URL}/gene/{gene_name}/all"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error fetching data for {gene_name}: {e}")
        return None

# --- Gene Search ---
def search_gene_id(gene_name):
    """
    Search for a gene by name and return its SubtiWiki gene ID (or None if not found).
    """
    endpoint = f"/search/"
    params = {"q": gene_name, "category": "Gene", "mode": "exact"}
    try:
        response = requests.get(BASE_URL + endpoint, params=params, timeout=20)
        response.raise_for_status()
        data = response.json()
        hits = data.get("data", {}).get("exact_hits", [])
        for hit in hits:
            if hit.get("category") == "Gene" and hit.get("name", "").lower() == gene_name.lower():
                return hit.get("id")
        hits = data.get("data", {}).get("partial_hits", [])
        for hit in hits:
            if hit.get("category") == "Gene" and hit.get("name", "").lower() == gene_name.lower():
                return hit.get("id")
        return None
    except Exception as e:
        print(f"Error searching for gene ID of {gene_name}: {e}")
        return None

# --- Fetch Regulations ---
def fetch_gene_regulations(gene_id):
    """
    Fetch regulatory information for a gene by SubtiWiki gene ID.
    Returns a list of regulators or None if not found.
    """
    endpoint = f"/gene/{gene_id}/regulations"
    try:
        response = requests.get(BASE_URL + endpoint, timeout=20)
        response.raise_for_status()
        data = response.json()
        return data.get("data", [])
    except Exception as e:
        print(f"Error fetching regulations for gene ID {gene_id}: {e}")
        return None

# --- Fetch Gene Metadata ---
def fetch_gene_metadata(gene_id):
    endpoint = f"/gene/{gene_id}"
    try:
        response = requests.get(BASE_URL + endpoint, timeout=20)
        response.raise_for_status()
        data = response.json()
        return data.get("data", {})
    except Exception as e:
        print(f"Error fetching metadata for gene ID {gene_id}: {e}")
        return None

# --- Menu ---
def main_menu():
    print("\nSubtiWiki API Client Menu:")
    print("1. Search for gene metadata")
    print("2. Fetch gene regulations")
    print("3. Fetch gene metadata (full)")
    print("0. Exit")
    choice = input("Enter your choice: ").strip()
    return choice

if __name__ == "__main__":
    while True:
        choice = main_menu()
        if choice == "0":
            print("Exiting.")
            break
        gene_name = input("Enter gene name: ").strip()
        gene_id = search_gene_id(gene_name)
        if not gene_id:
            print(f"Gene '{gene_name}' not found in SubtiWiki.")
            continue
        if choice == "1":
            print(f"Gene ID for {gene_name}: {gene_id}")
        elif choice == "2":
            regulations = fetch_gene_regulations(gene_id)
            print(f"Regulations for {gene_name}:")
            print(regulations)
        elif choice == "3":
            metadata = fetch_gene_metadata(gene_id)
            print(f"Metadata for {gene_name}:")
            print(metadata)
        else:
            print("Invalid choice. Please try again.") 