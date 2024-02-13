"""
Compare Ensembl and FOMOnet fasta to extract (1) protein ids which are
unique to FOMOnet and (2) peptides sequences of these proteins 
"""

import numpy as np

from Kiwi import reader, digestion


ENSEMBL = snakemake.input['ensembl']
FOMONET = snakemake.input['fomonet']
UNIQUE_IDS = snakemake.output['unique_ids']
PEPTIDE_SEQUENCES = snakemake.output['fomonet_peptides']

# Convert fasta files
ensembl_dict = reader.fasta_to_dict(ENSEMBL)
fomonet_dict = reader.fasta_to_dict(FOMONET)

# Find sequences which are predictions from FOMOnet
fomonet_ids = []
for id_, sequence in fomonet_dict.items():
    base_id = id_.split('-')[0]
    if base_id not in ensembl_dict:
        fomonet_ids.append(id_)
    else:
        # Make sure the sequences are not the same
        if sequence != ensembl_dict[base_id]:
            fomonet_ids.append(id_)

# Make dict {protein: sequence} with only fomonet proteins
only_fomonet_dict = {id_: fomonet_dict[id_] for id_ in fomonet_ids}

# Make a list of all peptide sequences obtained from fomonet cleaved proteins
dig = digestion.Digestion(only_fomonet_dict)
dig.cleave_proteins()
dig.check_peptide_uniqueness()
sequences = np.unique(dig.get_sequences()) # Remove duplicates

# Write fomonet protein ids
with open(UNIQUE_IDS, 'w', encoding='UTF-8') as f:
    for id_ in np.unique(fomonet_ids):
        f.write(f'{id_}\n')

# Write fomonet peptide sequences
with open(PEPTIDE_SEQUENCES, 'w', encoding='UTF-8') as f:
    f.write(', '.join(sequences))
