'''
Author: Gabriel Furtado Noll
Purpose: Change Ensembl IDs to Gene Names
'''

import json
import requests
import pandas as pd

# MODIFY path to data as needed
PATH_TO_DATA = 'SRP181622/SRP181622.tsv'

# Read the TSV file into a DataFrame
df = pd.read_csv(PATH_TO_DATA, sep='\t', index_col=0)
indexes = df.index.tolist()

# Function to split list into chunks of size n
def chunks(lst, n):
    for i in range(0, len(lst), n):
        print(f"Processing chunk starting at index {i}")
        yield lst[i:i+n]

# Function to fetch display names for a list of Ensembl IDs
def fetch_display_names(ensembl_ids):
    url = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    name_map = {}
    for block in chunks(ensembl_ids, 1000):  # Ensembl allows up to 1000 IDs/post
        payload = {"ids": block}
        r = requests.post(url, headers=headers, data=json.dumps(payload), timeout=60)
        r.raise_for_status()
        data = r.json()
        for _id, meta in data.items():
            if isinstance(meta, dict) and "display_name" in meta:
                name_map[_id] = meta["display_name"]
    return name_map

name_map = fetch_display_names(indexes)

# Keep Ensembl IDs as index and add a 'gene' column
df_with_gene = df.copy()
df_with_gene["gene"] = df_with_gene.index.map(name_map.get)
df_with_gene = df_with_gene.dropna(subset=['gene']) # Drop rows where gene name is not found
df_with_gene.to_csv('ensembl_and_gene.tsv', sep='\t', index=True)


# Option B: replace index with gene symbols where available
df_just_gene = df_with_gene.copy()
df_just_gene = df_just_gene.set_index('gene') # Drop rows where gene name is not found
df_just_gene.to_csv('gene.tsv', sep='\t', index=True)
