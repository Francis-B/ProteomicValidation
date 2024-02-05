library(MSstats)
library(MSstatsConvert)
library(magrittr)
library(gtools)
library(stringr)
library(argparse)

# ------------------------------------------------------------------------------
# Get snakemake files and params

to_list <- function(str_) {
  return(str_split_1(str_, " "))
}

parser <- ArgumentParser()

parser$add_argument("--pooled", "-p", help = "Pooled DIANN report")
parser$add_argument("--summaries", "-s", help = "DIANN summaries")
parser$add_argument("--annotation", "-a", help = "Annotation file")
parser$add_argument("--protein_level", "-q", help = "Protein level quant file")
parser$add_argument("--feature_level", "-f", help = "Feature level quant")
parser$add_argument("--comparaison", "-m", help = "Comparaison file")
parser$add_argument("--log", "-l", help = "Log file")

xargs <- parser$parse_args()
POOLED_DIANN <- xargs$pooled
SUMMARIES <- to_list(xargs$summaries)
ANNOTATION_PATH <- xargs$annotation
PROTEIN_LEVEL <- xargs$protein_level
FEATURE_LEVEL <- xargs$feature_level
COMPARISON <- xargs$comparaison
LOG <- xargs$log

# POOLED_DIANN <- "Documents/FomoNet/results/DIA-NN/pooled_report.tsv"
# SUMMARIES <- c("/home/francis/Documents/FomoNet/results/DIA-NN/ctr/report.stats.tsv",
#                "/home/francis/Documents/FomoNet/results/DIA-NN/hyp/report.stats.tsv",
#                "/home/francis/Documents/FomoNet/results/DIA-NN/rep/report.stats.tsv")
# ANNOTATION_PATH <- "/home/francis/Documents/FomoNet/results/DIA-NN/annotation.csv"
# PROTEIN_LEVEL <- "/home/francis/Documents/FomoNet/results/msstats_protein_quant.tsv"
# FEATURE_LEVEL <- "/home/francis/Documents/FomoNet/results/msstats_feature_quant.tsv"
# COMPARISON <- "/home/francis/Documents/FomoNet/results/msstats_comparison.tsv"
# LOG <- "/home/francis/Documents//FomoNet/logs/msstats.log"
# ------------------------------------------------------------------------------
# Create annotation file required by MSstats

# Extract run names, conditions and replicates num from DIA-NN report
annotation <- data.frame()
for (summary in SUMMARIES) {
  smry <- read.csv(summary, sep = "\t")
  fullnames <- smry$File.Name
  for (fn in fullnames) {
    # Get run name
    run <- tail(strsplit(fn, .Platform$file.sep, fixed = TRUE)[[1]], 1) %>%
      strsplit(split = ".", fixed = TRUE) %>%
      extract2(1) %>%
      extract(1)
    # Get replicate num
    num <- run %>%
      strsplit(split = "_") %>%
      extract2(1) %>%
      tail(1) %>%
      strsplit(".", fixed = TRUE) %>%
      extract2(1) %>%
      extract(1)
    # Get condition
    cond <- run %>%
      strsplit(split = "_") %>%
      extract2(1) %>%
      extract(3)


    # Add run name, condition and replicate num to annotation df
    annotation <- rbind(annotation,
                        data.frame(Run = run,
                                   Condition = cond,
                                   BioReplicate = num))
  }
}

write.csv(annotation, ANNOTATION_PATH, row.names = TRUE)

diann <- DIANNtoMSstatsFormat(POOLED_DIANN,
                              annotation = ANNOTATION_PATH,
                              qvalue_cutoff = 0.01,
                              pg_qvalue_cutoff = 0.01,
                              MBR = TRUE,
                              log_file_path = LOG)

# ------------------------------------------------------------------------------
# Proprocess results with MSstats

# Read data
data <- unique(as.data.frame(diann))

# use MSstats for protein summarization
summarized <- MSstats::dataProcess(diann,
                                   normalization = "equalizeMedians",
                                   logTrans = 2,
                                   nameStandards = c(""),
                                   featureSubset = "all",
                                   n_top_feature = NULL,
                                   summaryMethod = "TMP",
                                   censoredInt = "NA",
                                   MBimpute = TRUE,
                                   remove50missing = FALSE,
                                   maxQuantileforCensored = 0.999,
                                   use_log_file = FALSE)

# Create the contrast matrix
contrast_matrix <- rbind(c(1, -1, 0), c(1, 0, -1), c(0, 1, -1))
row.names(contrast_matrix) <- c("ctr vs hyp", "ctr vs rep", "hyp vs rep")
colnames(contrast_matrix) <- c("ctr", "hyp", "rep")

# Model-based comparison
model <- MSstats::groupComparison(contrast_matrix,
                                  summarized,
                                  use_log_file = FALSE)

# Get comparison results and save it to tsv
write.table(model$ComparisonResult, file = COMPARISON, sep = "\t")
write.table(summarized$ProteinLevelData,
            file = PROTEIN_LEVEL,
            sep = "\t", row.names = FALSE)
write.table(summarized$FeatureLevelData, 
            file = FEATURE_LEVEL,
            sep = "\t", row.names = FALSE)
