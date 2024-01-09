""" 
Compare ids and sequences from FOMOnet and Ensembl fasta files to make list
of sequences unique to FOMOnet. Then, check the results files to find which
protein groups contains at leat one of these unique sequences and which contains
only these unique sequences.
"""
import numpy as np
import pandas as pd
import utils

# Infiles
QUANT_FILE = snakemake.input.quant
ENSEMBL = snakemake.params.ensembl
FOMONET = snakemake.params.fomonet
CONDITIONS = snakemake.params.conditions
# Outfiles
UNIQUES_PG = snakemake.output.unique_pg
CONTAINS_PG = snakemake.output.contains_pg
UNIQUES_ID = snakemake.output.unique_id

# # Infiles
# ENSEMBL = '/home/francis/Documents/Université/MabLab/FOMOnet/canonical.fasta'
# FOMONET = '/home/francis/Documents/Université/MabLab/FOMOnet/fomonet.fasta'
# # QUANT_FILE = lambda cond: f'Results/renamed_fomonet_protein_groups_quant_{cond}.tsv'
# QUANT_FILE = '/home/francis/Documents/Université/MabLab/FOMOnet/Results/msstats_quant_protein.tsv'

# # Outfiles
# UNIQUES_PG = '/home/francis/Desktop/TEST_unique_fomonet_pg.txt' # Identified pg with only fomonet id
# CONTAINS_PG = '/home/francis/Desktop/TEST_contains_fomonet_pg.txt' # Identified pg with both fomonet and ensembl id
# UNIQUES_ID = '/home/francis/Desktop/TEST_unique_fomonet_id.txt' # Identified fomonet id which are unique

# Conditions name of the files
CONDITIONS = ['ctr', 'hyp', 'rep']

# Convert fasta files
ensembl_dict = utils.read_fasta(ENSEMBL)
fomonet_dict = utils.read_fasta(FOMONET)

# Find sequences which are predictions from FOMOnet
unique_ids = []
for id_, sequence in fomonet_dict.items():
    base_id = id_.split('-')[0]
    if base_id not in ensembl_dict:
        unique_ids.append(id_)
    else:
        # Make sure the sequences are not the same
        if sequence != ensembl_dict[base_id]:
            unique_ids.append(id_)

# Read quantification file to get identified protein groups
with open(QUANT_FILE, 'r', encoding='utf-8') as f:
    df = pd.read_csv(f, sep='\t').pivot(index='Protein',
                                        columns='originalRUN',
                                        values='LogIntensities')

protein_groups = df.index.to_numpy()

# Bool array to splice proteins groups which contains only ids unique to fomonet
is_only_fomonet = np.array([np.all([prot in unique_ids for prot in pg.split(';')])
                            for pg in protein_groups])

# Bool array to splice proteins groups which contains AT LEAST 1 id unique to fomonet
with_fomonet = np.array([np.any([prot in unique_ids for prot in pg.split(';')])
                                    for pg in protein_groups])

with open(UNIQUES_PG, 'w', encoding='UTF-8') as f:
    for pg in np.unique(protein_groups[is_only_fomonet]):
        f.write(f'{pg}\n')
with open(CONTAINS_PG, 'w', encoding='UTF-8') as f:
    for pg in np.unique(protein_groups[with_fomonet]):
        f.write(f'{pg}\n')
with open(UNIQUES_ID, 'w', encoding='UTF-8') as f:
    for id_ in np.unique(unique_ids):
        f.write(f'{id_}\n')
