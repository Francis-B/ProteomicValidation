# From the main report output of DIA-NN, count number of protein groups
# identified and create a quantification table for each database.

library(diann)

setwd("~/Documents/Université/MabLab/FOMOnet/")
out_dir <- file.path("Analysis", "Results")

database <- c("ensembl", "uniprot", "fomonet") # Used for path to files
conditions <- c("ctr", "hyp", "rep") # Used for path to files
replicatas <- seq(1, 8, 1)

# Loop through conditions and replicatas to extract number of proteins
# identified in each respective raw file
for (db in database){

  # Empty dataframe to save identified protein groups with db
  res <- data.frame(row.names = replicatas)

  for (cond in conditions){
    # Load DIA-NN results and get protein groups quantification df
    df <- diann_load(file.path("DIA-NN", paste0(db, "_", cond), "report.tsv"))
    protein_groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,],
                                   group.header = "Protein.Group",
                                   id.header = "Precursor.Id",
                                   quantity.header = "Precursor.Normalised") %>%
                      as.data.frame()

    # Create column for actual condition"s results
    res[cond] <- cbind(rep(NA, length(replicatas)))

    for (i in replicatas) {
      # From file name, get number of replicata
      colname <- colnames(protein.groups[i])
      file_name <- tail(strsplit(colname, "\\", fixed = TRUE)[[1]], 1)
      num_rep <- as.double(strsplit(tail(strsplit(file_name, "_")[[1]], 1),
                                    ".",
                                    fixed = TRUE)[[1]][1])
      # Add the number of identified proteins to the res df
      res[num_rep, cond] <- sum(!is.na(protein_groups[i]))
    }
  }

  # Write total protein groups identified to file
  write.table(res,
              file = file.path(out_dir, paste0(db, "_protein_groups.tsv")),
              sep = "\t")
  # Write protein.groups df to file for further analysis
  write.table(protein_groups,
              file = file.path(out_dir,
                               paste0(db, "_protein_groups_quant.tsv")),
              sep = "\t")
}
