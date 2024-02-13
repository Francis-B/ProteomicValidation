""" 
This pipeline run a DIA-MS analysis with FOMOnet predictions as protein databases. 
"""

configfile: "config.yaml"

rule all:
    input:
        "results/figures/upsetplot.svg",
        "results/identified_FOMOnet.xlsx",
        "results/figures/volcano_plot.svg",
        "results/fomonet_id.pkl"
        
rule get_fomonet_sequences:
    """ Compare Ensembl and FOMOnet fasta to extract (1) protein ids which are
        unique to FOMOnet and (2) peptides sequences of these proteins """
    input:
        ensembl = config['ensembl'],
        fomonet = config['fomonet']
    output:
        unique_ids = "results/unique_fomonet_id.txt", # All ids unique to FOMOnet
        fomonet_peptides = "results/fomonet_peptides.txt"
    script:
        "scripts/get_fomonet_sequences.py"


rule get_parent_proteins:
    """ Digest FOMOnet fasta file with trypsin using KIWI. Get all parent 
        proteins for each peptide obtained."""
    input:
        fasta = config['fomonet']
    output:
        digested = "results/peptides.tsv"
    script:
        "scripts/get_parent_proteins.py"
        
    # shell:
    #     "./scripts/Kiwi/digestion.py "
    #     "--miscleavages 1 "
    #     "--enzyme trypsin "
    #     "--unique " # Add boolean column to indicate if peptide is unique
    #     "--output {output.digested_fasta} "
    #     "{input.fasta}"

# rule convert_raw_to_mzml:
#     """ Convert raw files to mzML using the ThermoRawFileParser tool to run 
#         DIA-NN in linux """
#     input:
#         "data/raw/"
#     output:
#         expand("data/mzml/COVARIS_50_{{condition}}_{rep}.mzML",
#                 condition=config['conditions'], rep=config['replicates'])
#     params:
#         outdir = "data/mzml",
#         thermorawfileparser = config['thermorawfileparser_path']
#     log:
#         "logs/thermorawfileparser_{condition}.log"
#     shell:
#         "mono {params.thermorawfileparser} "
#         "-d {input} "
#         "-o {params.outdir}"


rule run_DIANN:
    """ Analyse mzML files with DIA-NN. Make a separate run for each condition """
    input:
        mzml = expand('data/mzml/COVARIS_50_{{condition}}_{rep}.mzML',
                        rep=config['replicates']),
        fomonet_peptides = "results/fomonet_peptides.txt"
    output:
        "results/DIA-NN/{condition}/report.tsv"
    params:
        fasta = config['fomonet'],
        flagged = lambda wildcards, input: ' --f '.join(input.mzml),
        lib = "results/DIA-NN/{condition}/report-lib.tsv",
        diann = config['diann_path']
    threads:
        5
    log:
        "logs/DIA-NN_{condition}.log"
    shell:
        "{params.diann} "  # Path to DIA-NN executable
        "--f {params.flagged} " # Path to mzML files
        "--out-lib {params.lib} "   # Path to generated spectral library
        "--out {output} "
        "--fasta {params.fasta} "
        "--qvalue 0.01 "
        "--matrices "
        "--gen-spec-lib "   # Generate spectral library from fasta
        "--predictor "
        "--fasta-search "
        "--min-fr-mz 200 "
        "--max-fr-mz 1800 "
        "--met-excision "
        "--cut K*,R*,!*P "  # Trypsin digestion
        "--missed-cleavages 1 "
        "--min-pep-len 7 "
        "--max-pep-len 30 "
        "--min-pr-mz 300 "
        "--max-pr-mz 1800 "
        "--min-pr-charge 1 "
        "--max-pr-charge 4 "
        "--unimod4 "
        "--int-removal 0 "
        "--report-lib-info "  # Add info about precursors and fragments to report
        "--reanalyse "     # Match between runs
        "--smart-profiling "
        "--pg-level 0 "
        "--peak-center "
        "--no-ifs-removal "
        "--threads {threads} "
        # "--vis 15, $(cat {input.fomonet_peptides}) " # Get XICs for FOMOnet peptides


rule DIANN_to_MSstats:
    """ Transform DIA-NN output for MSstats compatibility. Concretely, merge 
        DIA-NN report files and change Protein.Group for Protein.Names 
        column """
    input:
        reports = expand("results/DIA-NN/{condition}/report.tsv",
                         condition=config['conditions'])
    output:
        pooled = "results/DIA-NN/pooled_report.tsv"
    params:
        header = lambda wildcards, input: input.reports[0], # Extract header from first report file
        reports = lambda wildcards, input: ' '.join(input.reports)
    shell:
        """                                   
        head -n 1 {params.header} |awk -v OFS='\t' '{{$3=$5; print $0}}' |cut -f 5 --complement > {output.pooled}
        for file in {params.reports}
        do
            cat $file |tail -n +2 |cut -f 5 --complement >> {output.pooled}
        done
        """


rule MSstats_normalize_intensities:
    """ With MSstats, normalize intensities and calculate fold changes 
        between conditions """
    input:
        pooled = "results/DIA-NN/pooled_report.tsv"
    output:
        comparison = "results/msstats_comparaison.tsv",
        protein_level = "results/msstats_protein_quant.tsv",
        feature_level = "results/msstats_feature_quant.tsv"
    params:
        summaries = expand("results/DIA-NN/{condition}/report.stats.tsv",
                            condition=config['conditions']),
        annotation = "results/DIA-NN/annotation.csv"
    conda:
        "envs/R.yaml"''
    log:
        "logs/msstats.log"
    shell:
        "Rscript scripts/msstats_process_intensities.R "
        "--pooled {input.pooled} "
        "--comparaison {output.comparison} "
        "--protein_level {output.protein_level} "
        "--feature_level {output.feature_level} "
        "--summaries '{params.summaries}' "
        "--annotation {params.annotation} "
        "--log {log}" 
        

rule get_identified_predictions:
    """ Find (1) identified protein groups with only FOMOnet predictions and (2)
        identified protein groups with at least one FOMOnet prediction """
    input:
        quant = "results/msstats_protein_quant.tsv",
        unique_ids = "results/unique_fomonet_id.txt"
    output:
        unique_pg = "results/unique_fomonet_pg.txt", # Only FOMOnet ids in pg
        contains_pg = "results/contains_fomonet_pg.txt", # At least one FOMOnet id in pg (but not all)
    script:
        "scripts/get_fomonet_pg.py"


rule make_volcano_plot:
    """ Make volcano plots from MSStats comparaison data """
    input:
        comparaison = "results/msstats_comparaison.tsv",
        unique_fomonet_pg = "results/unique_fomonet_pg.txt",
        contains_fomonet_pg = "results/contains_fomonet_pg.txt"
    output:
        plot = "results/figures/volcano_plot.svg"
    conda:
        "envs/R.yaml"
    shell:
        "Rscript scripts/volcano_plot.R "
        "--comparaison {input.comparaison} "
        "--unique-fomonet-pg {input.unique_fomonet_pg} "
        "--contains-fomonet-pg {input.contains_fomonet_pg} "
        "--plot {output.plot}"


rule analyse_quantification:
    """ From normalized intensities with MSstats, plot the intensities
        distributions and make an upsetplot of conditions in which FOMOnet pg 
        are identified """
    input:
        quant = "results/msstats_protein_quant.tsv",
        unique_ids = "results/unique_fomonet_pg.txt"
    output:
        upsetplot = "results/figures/upsetplot.svg",
        quant_plot = "results/figures/quant_plot.svg",
        fomonet_id = "results/fomonet_id.pkl",  # TODO: mark it back as temp
        quant_pg = temp("results/quant_pg.pkl")
    params:
        conditions = config['conditions'],
        allowed_na = config['allowed_na']
    conda:
        "envs/py_plot.yaml"
    script:
        "scripts/analyze_quantification.py"


rule list_identified_fomonet:
    """ List all pg uniquely composed of FOMOnet predictions and use pybiomart
        to query Ensembl for information about them """
    input:
        fomonet_id = "results/fomonet_id.pkl",
        quant_pg = "results/quant_pg.pkl"
    output:
        identified_xlsx = "results/identified_FOMOnet.xlsx"
    params:
        conditions = config['conditions']
    conda:
        "envs/py_biomart.yaml"
    script:
        "scripts/list_identified_fomonet.py"
