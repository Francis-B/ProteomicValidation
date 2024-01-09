""" 
Set of generic functions used in the analysis of the results. 
"""

from collections import defaultdict
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore") # Filter out np.nanmean RuntimeWarning

# Function to convert fasta files to dictionary
def read_fasta(fasta):
    fasta_dict = defaultdict(str)
    with open(fasta, 'r', encoding='UTF-8') as f:
        for line in f:
            # Get the transcript id
            if line.startswith('>'):
                id_ = line.split(' ')[0].strip('>\n')
            # Add the sequence
            else:
                fasta_dict[id_] += line.strip('\n')
    return fasta_dict


def extract_msstats_quant(file, conditions, allowed_na=1):
    """ Read MSstats quantification file and return dictionary with an array of
        identified pg names and an array of quantification for each condition."""

    with open(file, 'r', encoding='utf-8') as f:
        df = pd.read_csv(f, sep='\t').pivot(index='Protein',
                                            columns='originalRUN',
                                            values='LogIntensities')

    pg = {cond: {'names': [], 'quants': []} for cond in conditions}

    for cond in conditions:

        # Get colnames corresponding to the condition
        colnames = [col for col in df.columns
                    if col.split('_')[2] == cond]

        # Calculate mean quantification
        quant_arrs = [df[col].to_numpy() for col in colnames]
        mean_quants = np.nanmean(quant_arrs, axis=0)

        # Filter out NA values and add log2 transformed mean quantifications to dict
        na_filter = np.sum(np.array([np.isnan(arr) for arr in quant_arrs]), axis=0) <= allowed_na
        pg[cond]['quants'] = mean_quants[na_filter]
        pg[cond]['names'] = df.index.to_numpy()[na_filter]

    return pg

def extract_diann_quant(file, rep_of_interest, conditions, allowed_na=1):
    """ Read DIA-NN quantification file and return dictionary with an array of
        identified pg names and an array of quantification for each condition.
        """
    pg_dict = {cond: {'names': [], 'quants': []} for cond in conditions}
    for cond in conditions:
        # Read the protein groups quantification file
        with open(file(cond), 'r', encoding='UTF-8') as f:
            raw_res = np.rot90(np.loadtxt(f, delimiter='\t', dtype=str))

        protein_groups = np.array([pg.strip('"') for pg in raw_res[-1]][1:])  # Protein groups array
        quant = raw_res[0:-1]   # Quantification values arrays

        # Create a 2d np array of quantifications values for the 3 best replicates
        rep_ids = [int(arr[0].split('\\')[-1].split('_')[-1].split('.')[0]) for arr in quant]
        quants_arr = np.array([res_arr[1:] for rep_id, res_arr in zip(rep_ids, quant)
                    if rep_id in rep_of_interest[cond]])

        # Correct type of NA values and convert all values to float
        quants_arr = np.where(quants_arr == "NA", np.nan, quants_arr).astype(float)

        # Calculate mean quantifications and log2 transform it
        mean_quants = np.log2(np.nanmean(quants_arr, axis=0))

        # Filter out NA values and add log2 transformed mean quantifications to dict
        na_filter = np.sum(np.array([np.isnan(arr) for arr in quants_arr]), axis=0) <= allowed_na
        pg_dict[cond]['quants'] = mean_quants[na_filter]
        pg_dict[cond]['names'] = protein_groups[na_filter]

    return pg_dict


def get_fomonet_pg(pg_dict, unique_ids):
    """ From a dictionary containing all protein groups identified in each
        condition, return a dictionary of bool arrays to splice the pg which
        contains only ids from FOMOnet. """
    is_unique_fomonet = {cond: [] for cond in pg_dict.keys()}
    for cond in pg_dict.keys():
        # Bool array to splice proteins groups which contains only ids unique to fomonet
        pg_list = pg_dict[cond]['names']
        is_unique_fomonet[cond] = np.array([np.all([prot in unique_ids for prot in pg.split(';')])
                                            for pg in pg_list])

    return is_unique_fomonet

def get_condition_unique_pg(pg_dict):
    """ From a dictionary containing all protein groups identified in each
        condition, return a dictionary of bool arrays to splice the pg which
        are unique to each condition. """
    # Get the pg which are unique to each condition
    is_unique_cond = {cond: [] for cond in pg_dict.keys()}
    for cond in pg_dict.keys():
        pg_other_cond = np.unique([np.concatenate([pg_cond['names']
                                                   for key, pg_cond in pg_dict.items()
                                                   if key != cond])])
        is_unique_cond[cond] = np.in1d(pg_dict[cond]['names'], pg_other_cond, invert=True)

    return is_unique_cond
