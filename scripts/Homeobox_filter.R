# Set the working directory
setwd("path/to/directory")

# Load necessary libraries
library(TCGAbiolinks) # For accessing and analyzing TCGA data
library(sesame)       # For handling Illumina methylation data
library(ExperimentHub)
library(SummarizedExperiment)
library(BiocGenerics)
library(AnnotationHub)
library(BiocFileCache)
library(sesameData)   # Provides genome annotations for sesame
library(data.table)   # Data manipulation
library(readr)        # Reading and writing data

# Load the preprocessed beta values (methylation matrix)
beta_values <- read_rds("file.rds")

# Read the list of Homeobox genes from a CSV file
unedited_list <- read.csv("HoxGenes.csv", header = FALSE)

# Remove every other row (odd-indexed rows)
todelete <- seq(0, nrow(unedited_list), 2)
shortened_list <- unedited_list[-todelete, ]

# Extract gene names by removing the first three characters
HoxGenes <- substring(shortened_list, 4)

# Optionally, save the cleaned list of gene names
# write.csv(HoxGenes, "HomeoboxList.csv")

# Initialize data structures for storing results
total_list <- beta_values
final_list <- data.frame()
name_list <- list()

# Prepare empty dataframe for later use
name_df <- t(as.data.frame(name_list))
final_list_2 <- cbind(name_df, final_list)

# Convert Homeobox gene names to a list
HoxGenes <- as.list(HoxGenes)

# Exclude specific genes (based on indices) from the list
HoxGenes1 <- HoxGenes[-c(57, 77:88, 91:93, 109:112, 119, 224, 225, 226:251, 257:354)]

# Initialize variables for looping and storing probe lists
Probe_List <- list()

# Loop through each gene to retrieve probes associated with it
for (p in seq_along(HoxGenes1)) {
  gene <- HoxGenes1[[p]]
  message(paste0("Working on ", gene, " now."))
  
  # Retrieve probes within 10kb upstream and downstream of the gene
  probes <- sesameData_getProbesByGene(
    gene,
    platform = "HM450", # Specify Illumina platform
    upstream = 10000,
    downstream = 10000,
    genome = "hg38"     # Specify genome version
  )
  
  # Extract probe names and add to the list
  probes1 <- list(probes@ranges@NAMES)
  Probe_List[[p]] <- probes1
}

# Remove specific entries from the list (optional cleanup)
Homeobox_list <- Probe_List[-c(122)]

# Save the list of probes for each gene
write_rds(Homeobox_list, "Probes_by_Gene.rds")

# Match probes in Homeobox_list with the methylation data (total_list)
for (d in seq_along(Homeobox_list)) {
  for (e in seq_along(Homeobox_list[[d]][[1]])) {
    if (length(grep(Homeobox_list[[d]][[1]][e], total_list[, 8])) != 0) {
      final_list <- rbind(final_list, total_list[grep(Homeobox_list[[d]][[1]][e], total_list[, 8]), ])
    }
  }
}

# Optionally save the final list of matched data
# write.csv(final_list, "DMR_final_list_BRCA450b.csv")

# Load Homeobox gene names for naming the results
name_genes <- read.csv("Homeobox_names_all.csv")
names(Homeobox_list) <- name_genes[, 2]

# Extract probes from final_list
probes5 <- final_list$Probe_ID
names1 <- list()
probes6 <- list()

# Assign gene names to probes
for (d in seq_along(Homeobox_list)) {
  for (f in probes5) {
    if (f %in% Homeobox_list[[d]][[1]]) {
      names1 <- append(names1, names(Homeobox_list)[d])
      probes6 <- append(probes6, f)
    }
  }
}

# Combine gene names with the final list
names_1 <- t(as.data.frame(names1))
final_with_names <- cbind(names_1, final_list)
