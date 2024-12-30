# Load required libraries
library(reshape2)
library(ggplot2)
library(stringr)
library(readr)
library(dplyr)

# Function to load and process CSV files
load_and_select_columns <- function(filepath) {
  read.csv(filepath, header = TRUE)[, c(1, 3, 15:17)]
}

# List of cancer types and corresponding file paths
cancer_types <- c("UCEC", "THCA", "BRCA", "COAD", "BLCA", "CESC", "HNSC", 
                  "ESCA", "PAAD", "PRAD", "LUAD", "LUSC", "CHOL", "KIRP", 
                  "KIRC", "LIHC")
csv_files <- paste0(cancer_types, "/Final_list_", cancer_types, "_with_domains.csv")
rds_files <- paste0(cancer_types, "/", cancer_types, "_beta_merged.rds")

# Load all CSV and RDS files into lists
csv_data <- lapply(csv_files, load_and_select_columns)
rds_data <- lapply(rds_files, readRDS)

# Assign column names uniformly
set_column_names <- function(df) {
  colnames(df) <- c("Num", "Seg_ID", "Gene", "Domain", "Tissue")
  return(df)
}
csv_data <- lapply(csv_data, set_column_names)

# Function to calculate summary statistics
calculate_summary <- function(csv_df, rds_df, tissue_name) {
  seg_count <- dim(csv_df %>% count(Seg_ID))[1]
  total_seg_count <- dim(rds_df %>% count(Seg_ID))[1]
  missing_count <- total_seg_count - seg_count
  percentage <- (seg_count / total_seg_count) * 100
  data.frame(`Homeobox DMRs` = seg_count, `DMRs` = missing_count, 
             Percentage = percentage, Tissue = tissue_name)
}

# Generate summary data for all cancer types
cancer_summaries <- mapply(calculate_summary, csv_data, rds_data, cancer_types, SIMPLIFY = FALSE)
summary_df <- do.call(rbind, cancer_summaries)

# Example of filtering and counting for specific cancer types
breast_total <- rds_data[[which(cancer_types == "BRCA")]]
breast_segments <- csv_data[[which(cancer_types == "BRCA")]] %>% count(Seg_ID)
breast_multiple_segments <- breast_segments %>% filter(n > 1)

# Save final output or perform further analysis
print(summary_df)
