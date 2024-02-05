""" 
Digest FOMOnet fasta file with trypsin using KIWI. Get all parent 
proteins for each peptide obtained.
"""

from collections import defaultdict
import numpy as np

from Kiwi.digestion import Digestion

FASTA = snakemake.input['fasta']
DIGESTED = snakemake.output['digested']

# Perform digestion, check which peptides are unique and write to file using Kiwi
digestion = Digestion(FASTA)
digestion.cleave_proteins()
digestion.check_peptide_uniqueness()
digestion.write(outdir=DIGESTED)    # Write temporary file (will be overwritten)

# Read file
with open(DIGESTED, 'r', encoding='utf-8') as f:
    peptides, proteins, _, uniqueness = np.genfromtxt(f, delimiter=', ',
                                                      dtype=str, unpack=True)

# Make list of all parent proteins for each peptide
parent_proteins = defaultdict(list)
for pep, prot in zip(peptides, proteins):
    parent_proteins[pep].append(prot)

# Write to file
with open(DIGESTED, 'w', encoding='utf-8') as f:
    for pep, prot, uniq in zip(peptides, proteins, uniqueness):
        if uniq == 'unique':
            pp = prot
        else:
            pp = ','.join(parent_proteins[pep])
        f.write('\t'.join([pep, prot, uniq, pp]) + '\n')
