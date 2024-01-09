import copy
import pickle
from collections import defaultdict
import pandas as pd

import pybiomart

FOMONET_ID = snakemake.input.fomonet_id
QUANT_PG = snakemake.input.quant_pg
IDENTIFIED_XLSX = snakemake.output.identified_xlsx
CONDITIONS = snakemake.params.conditions

# def get_transcript_info(transcript_list):
#     """ 
#     From a txt file listing the transcript ids, make a biomart query to
#     extract the gene name and biotype of each transcript. Return a dictionary.
#     """
#     clean_ids = {id_.split('-')[0]: id_ for id_ in transcript_list}

#     # Extract identified transcript information from Ensembl
#     dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl',
#                                 host='http://www.ensembl.org')

#     # Query everything because transcript_id as filters doesn't work
#     query = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name',
#                                     'gene_biotype', 'ensembl_transcript_id'])

#     mask = query['Transcript stable ID'].isin(clean_ids.keys())
#     filtered = query[mask] # Filter out unrelated transcripts

#     # Create a dictionary of results
#     res = defaultdict(dict)

#     for ensembl_id, biotype, gene_name, gene_id in zip(filtered['Transcript stable ID'],
#                                               filtered['Gene type'],
#                                               filtered['Gene name'],
#                                               filtered['Gene stable ID']):
#         res[clean_ids[ensembl_id]]['Associated gene biotype'] = biotype
#         res[clean_ids[ensembl_id]]['Associated gene name'] = gene_name
#         res[clean_ids[ensembl_id]]['Associated gene id'] = gene_id

#     return res

# # Get transcript information from biomart
# biomart_dict = utils.get_transcript_info(id_request)

# Extract identified transcript information from Ensembl
dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl',
                            host='http://www.ensembl.org')

# Query everything because transcript_id as filters doesn't work
query = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name',
                                'gene_biotype', 'ensembl_transcript_id'])

with open(FOMONET_ID, 'rb') as f:
    id_request = pickle.load(f)

clean_ids = {id_.split('-')[0]: id_ for id_ in id_request}

# Keep only the line which contains identified FOMOnet predictions
mask = query['Transcript stable ID'].isin(clean_ids.keys())
filtered = query[mask]

# Create a dictionary of results
biomart_dict = defaultdict(dict)

for ensembl_id, biotype, gene_name, gene_id in zip(filtered['Transcript stable ID'],
                                            filtered['Gene type'],
                                            filtered['Gene name'],
                                            filtered['Gene stable ID']):
    biomart_dict[clean_ids[ensembl_id]]['Associated gene biotype'] = biotype
    biomart_dict[clean_ids[ensembl_id]]['Associated gene name'] = gene_name
    biomart_dict[clean_ids[ensembl_id]]['Associated gene id'] = gene_id


keys = ['Protein group', 'Ensembl id', 'Associated gene id',
        'Associated gene name', 'Associated gene biotype', 'ctr', 'hyp', 'rep']
tmp_dict = {key: pd.NA for key in keys}

with open(QUANT_PG, 'rb') as f:
    quants_dict = pickle.load(f)

res = defaultdict(lambda: copy.deepcopy(tmp_dict))
for cond in CONDITIONS:
    for pg in quants_dict[cond]:
        for prot in pg.split(';'):
            # res[prot].setdefault(copy.deepcopy(tmp_dict))
            res[prot]['Protein group'] = pg
            res[prot][cond] = quants_dict[cond][pg]
            res[prot]['Ensembl id'] = prot.split('-')[0]
            # Add biomart information, skip if not found
            try:
                res[prot]['Associated gene id'] = biomart_dict[prot]['Associated gene id']
                res[prot]['Associated gene name'] = biomart_dict[prot]['Associated gene name']
                res[prot]['Associated gene biotype'] = biomart_dict[prot]['Associated gene biotype']
            except KeyError:
                pass

df = pd.DataFrame.from_dict(res, orient='index')
df.index.name = 'FOMOnet id'
df.to_excel(IDENTIFIED_XLSX)
