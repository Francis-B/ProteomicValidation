""" 
This pipeline run a DIA-MS analysis with FOMOnet predictions as protein databases. 
"""

configfile: "config.yaml"

rule all:
    input:
        "results/Figures/upsetplot.svg",
        "results/identified_FOMOnet.xlsx"


rule convert_raw_to_mzml:
    """ Convert raw files to mzML using the ThermoRawFileParser tool to run 
        DIA-NN in linux """

    input:
        "data/raw/"

    output:
        expand("data/mzml/COVARIS_50_{condition}_{rep}.mzML",
                condition=config['conditions'], rep=config['replicates'])
    params:
        outdir = "data/mzml"

    shell:
        "mono ~/ThermoRawFileParser1.4.2/ThermoRawFileParser.exe "
        "-d {input} "
        "-o {params.outdir}"


rule run_DIANN:
    """ Analyse mzML files with DIA-NN. Make a separate run for each condition """

    input:
        mzml = expand('data/mzml/COVARIS_50_{{condition}}_{rep}.mzML',
                        rep=config['replicates'])

    output:
        "results/DIA-NN/{condition}/report.tsv"

    params:
        fasta = config['fomonet'],
        flagged = lambda wildcards, input: ' --f '.join(input.mzml),
        lib = "results/DIA-NN/{condition}/report-lib.tsv",

    threads:
        8

    shell:
        "/usr/diann/1.8.1/diann-1.8.1 "
        "--f {params.flagged} "
        "--out-lib {params.lib} "
        "--out {output} "
        "--fasta {params.fasta} "
        "--qvalue 0.01 "
        "--matrices "
        "--gen-spec-lib "
        "--predictor "
        "--fasta-search "
        "--min-fr-mz 200 "
        "--max-fr-mz 1800 "
        "--met-excision "
        "--cut K*,R*,!*P "
        "--missed-cleavages 1 "
        "--min-pep-len 7 "
        "--max-pep-len 30 "
        "--min-pr-mz 300 "
        "--max-pr-mz 1800 "
        "--min-pr-charge 1 "
        "--max-pr-charge 4 "
        "--unimod4 "
        "--int-removal 0 "
        "--reanalyse "
        "--smart-profiling "
        "--pg-level 0 "
        "--peak-center "
        "--no-ifs-removal "
        "--threads {threads}"


rule DIANN_to_MSstats:
    """ Merge DIA-NN output files and change Protein.Group for Protein.Names 
        column for downstream analysis with MSstats """

    input:
        reports = expand('results/DIA-NN/{condition}/report.tsv',
                         condition=config['conditions'])

    output:
        pooled = "results/DIA-NN/DIANN_pooled_report.tsv"

    params:
        header = lambda wildcards, input: input.reports[0],
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
        pooled = "DIA-NN/DIANN_pooled_report.tsv"
    
    output:
        comparaison = "results/msstats_comparaison.tsv",
        quant = "results/msstats_quant.tsv"

    params:
        conditions = config['conditions'],
        summaries = expand("results/DIA-NN/{condition}/report.stats.tsv",
                            condition=config['conditions']),
        annotation = temp("DIA-NN/annotation.csv")

    shell:
        "Rscript scripts/msstats_process_intensities.R "
        "--pooled {input.pooled} "
        "--comparaison {output.comparaison} "
        "--quant {output.quant} "
        "--summaries '{params.summaries}' "
        "--annotation {params.annotation} " 
        "--conditions '{params.conditions}' "
        

rule get_identified_predictions:
    """ Find which FOMOnet ids belong to new proteins (not found in Ensembl)
        and list proteins groups that contain them """

    input:
        quant = "results/msstats_quant.tsv"

    output:
        unique_pg = "results/unique_fomonet_pg.txt", # Only FOMOnet ids in pg
        contains_pg = "results/contains_fomonet_pg.txt", # At least one FOMOnet id in pg (but not all)
        unique_id = "results/unique_fomonet_id.txt" # All ids unique to FOMOnet

    params:
        ensembl = config['ensembl'],
        fomonet = config['fomonet'],
        conditions = config['conditions']

    script:
        "scripts/get_fomonet_pg.py"


rule make_volcano_plot:
    """ Make volcano plots from MSStats comparaison data"""

    input:
        comparaison = "results/msstats_comparaison.csv",
        unique_fomonet_pg = "results/unique_fomonet_pg.txt",
        contains_fomonet_pg = "results/contains_fomonet_pg.txt"

    output:
        plot = "results/Figures/volcano_plot.svg"

    shell:
        "Rscript scripts/volcano_plot.R "
        "--comparaison {input.comparaison} "
        "--unique-fomonet-pg {input.unique_fomonet_pg} "
        "--contains-fomonet-pg {input.contains_fomonet_pg} "
        "--plot {output.plot}"


rule analyse_quantification:
    """ From normalized intensities with MSstats, plot the intensities
        distributions and make an upsetplot of conditions in which FOMOnet pg 
        are identifiedextract information about identified proteins """
        
    input:
        quant = "results/msstats_quant.tsv",
        unique_ids = "results/unique_fomonet_pg.txt",

    output:
        upsetplot = "results/Figures/upsetplot.svg",
        quant_plot = "results/Figures/quant_plot.svg",
        fomonet_id = temp("results/fomonet_id.pkl"),
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

