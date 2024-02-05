library(magrittr)
library(gridExtra)
library(ggplot2)
library(argparse)
library(svglite)
library(ggrepel)

# ------------------------------------------------------------------------------
# Get snakemake files and params

# parser <- ArgumentParser()
# 
# parser$add_argument("--comparison", "-c", help = "Comparison file")
# parser$add_argument("--unique-fomonet-pg", "-u", help = "Unique pg to fomonet")
# parser$add_argument("--contains-fomonet-pg", "-f", help = "Contains fomonet pg")
# parser$add_argument("--plot", "-p", help = "Plot file")
# 
# xargs <- parser$parse_args()
# COMPARISON <- xargs$comparison
# UNIQUES_FOMO <- xargs$unique_fomonet_pg
# CONTAINS_FOMO <- xargs$contains_fomonet_pg
# PLOT <- xargs$plot

COMPARISON <- "/home/francis/Documents/FomoNet/results/msstats_comparison.tsv"
UNIQUES_FOMO <- "/home/francis/Documents/FomoNet/results/unique_fomonet_pg.txt"
CONTAINS_FOMO <- "/home/francis/Documents/FomoNet/results/contains_fomonet_pg.txt"
PLOT <- "/home/francis/Documents/FomoNet/results/figures/volcano_plot.svg"


# ------------------------------------------------------------------------------
# Read preprocessed data with MSstats and filter out rows with NA pvalue
raw_res <- read.csv(COMPARISON, sep = "\t")
comp <- raw_res[!is.na(raw_res$pvalue), ]

# Get pg with only ids unique to fomonet
unique_fomonet_pg <- read.table(UNIQUES_FOMO, sep = "\n", header = FALSE)
with_fomonet_pg <- read.table(CONTAINS_FOMO, sep = "\n", header = FALSE)

# Create the contrast matrix
contrast_matrix <- rbind(c(1, -1, 0), c(1, 0, -1), c(0, 1, -1))
row.names(contrast_matrix) <- c("ctr vs hyp", "ctr vs rep", "hyp vs rep")
colnames(contrast_matrix) <- c("ctr", "hyp", "rep")

# Make volcano plots
p_list <- lapply(seq_len(dim(contrast_matrix)[1]), function(i) {
  
  # Get data for comparison
  comparison <- row.names(contrast_matrix)[[i]]
  data <- comp[comp$Label == comparison,]
  data$log10pvalue <- -log10(data$pvalue)

  # Colour of point (highlight pg unique to fomonet and containing fomonet)
  data$colour <- rep("PG with no FOMOnet", nrow(data))
  data$colour[data$Protein %in% with_fomonet_pg$V1] <- "PG with at least one FOMOnet"
  data$colour[data$Protein %in% unique_fomonet_pg$V1] <- "PG with only FOMOnet"

  # Get significant value
  fomonet_pg <- c("PG with only FOMOnet", "PG with at least one FOMOnet")
  sign_data <- data[data$pvalue < 0.05 &
                      (data$log2FC > log2(1.2) | data$log2FC < -log2(1.2)) &
                      (data$colour %in% fomonet_pg), ] %>%
    subset(select = c("Protein", "log2FC", "pvalue", "log10pvalue", "colour"))

  sign_fomonet <- sign_data[sign_data$colour == "PG with only FOMOnet", ]
  p <- ggplot(data,
              mapping = aes(x = log2FC, y = log10pvalue, colour = colour)) +
    xlab("log2FC") +
    ylab("-log10(pvalue)") +
    ylim(c(0, 3.5)) +
    geom_point(show.legend = FALSE) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed",
               colour = "#7E7973") +
    geom_vline(xintercept = log2(1.2) * c(-1, 1),
               linetype = "dashed",
               colour = "#7E7973") +
    ggtitle(comparison) +
    scale_colour_manual("",
                        values = c("PG with only FOMOnet" = adjustcolor("#245AB3", alpha = 0.7),
                                    "PG with at least one FOMOnet" = adjustcolor("#B37D24", alpha = 0.7),
                                    "PG with no FOMOnet" = adjustcolor("grey65", alpha.f = 0.5))) +
    geom_text_repel(data = sign_fomonet,
                    mapping = aes(x = sign_fomonet[, "log2FC"],
                                  y = sign_fomonet[,"log10pvalue"]),
                    size = 2,
                    colour = "black",
                    show.legend = FALSE,
                    label = sign_fomonet[, "Protein"]) +
    theme(legend.position = "right",
          panel.background = element_rect(fill = "#F9F5EE"),
          plot.title = element_text(hjust = 0.5))

  return(p)
})

g <- arrangeGrob(grobs = p_list, ncol = dim(contrast_matrix)[1])
ggsave(PLOT, g, width = 30, height = 13, units = "cm", dpi = 300)

