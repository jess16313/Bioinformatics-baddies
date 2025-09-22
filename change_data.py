import pandas as pd
import ensembl_rest

PATH_TO_DATA = 'SRP181622/SRP181622.tsv'

# Read the TSV file into a DataFrame
df = pd.read_csv(PATH_TO_DATA, sep='\t', index_col=0)

# Get the list of indexes (Ensembl IDs)
indexes = df.index.tolist()

# Initialize an empty list to store gene names
genes = []

for i in indexes:
    # Lookup the gene name using Ensembl REST API and append to list
    try:
        print(f"looking up {i}")
        gene_name = ensembl_rest.lookup(i)['display_name']
        print(f"found {gene_name}")
        genes.append(gene_name)
    # If there's ID not found, drop that row from the DataFrame
    except ensembl_rest.HTTPError as e:
        df.drop(i, inplace=True)

# # If the lengths match, add the gene names as a new column and save to a new TSV file
# if len(genes) == len(df.index.tolist()):
# Ensembl and Genes
df_gene = df.copy()
df_gene['gene'] = genes
df_gene.to_csv('ensembl_and_gene.tsv', sep='\t', index=True)

# Just genes
df.index = genes
df.to_csv('ensembl_and_gene.tsv', sep='\t', index=True)
# # If lengths don't match, print the lengths for debugging
# else:
#     print(len(genes))
#     print(len(df.index.tolist()))
