""" 
Find (1) identified protein groups with only FOMOnet predictions and (2)
identified protein groups with at least one FOMOnet prediction
"""
import numpy as np
import pandas as pd

# Infiles
QUANT_FILE = snakemake.input.quant
UNIQUE_IDS = snakemake.input.unique_ids

# Outfiles
UNIQUES_PG = snakemake.output.unique_pg
CONTAINS_PG = snakemake.output.contains_pg

with open(UNIQUE_IDS, 'r', encoding='UTF-8') as f:
    fomonet_ids = [line.strip() for line in f]

# Read quantification file to get identified protein groups
with open(QUANT_FILE, 'r', encoding='utf-8') as f:
    df = pd.read_csv(f, sep='\t').pivot(index='Protein',
                                        columns='originalRUN',
                                        values='LogIntensities')

protein_groups = df.index.to_numpy()

# Bool array to splice proteins groups which contains only ids unique to fomonet
is_only_fomonet = np.array([np.all([prot in fomonet_ids for prot in pg.split(';')])
                            for pg in protein_groups])

# Bool array to splice proteins groups which contains AT LEAST 1 id unique to fomonet
with_fomonet = np.array([np.any([prot in fomonet_ids for prot in pg.split(';')])
                                    for pg in protein_groups])

with open(UNIQUES_PG, 'w', encoding='UTF-8') as f:
    for pg in np.unique(protein_groups[is_only_fomonet]):
        f.write(f'{pg}\n')

with open(CONTAINS_PG, 'w', encoding='UTF-8') as f:
    for pg in np.unique(protein_groups[with_fomonet]):
        f.write(f'{pg}\n')
