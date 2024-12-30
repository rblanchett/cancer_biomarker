# Load necessary libraries for genomic data analysis
library(TCGAbiolinks)        # Access and process TCGA data
library(sesame)              # DNA methylation analysis
library(ExperimentHub)       # Access experiment-related data
library(SummarizedExperiment) # Manage high-throughput assay data
library(BiocGenerics)        # Core infrastructure for Bioconductor
library(AnnotationHub)       # Retrieve annotation data
library(BiocFileCache)       # Cache files for reuse
library(sesameData)          # Provides methylation data resources
library(readr)               # Read and write data efficiently

# Fetch the list of all GDC (Genomic Data Commons) projects
GDCprojects = getGDCprojects()

# Query TCGA for breast cancer (BRCA) methylation data using Illumina Human Methylation 450
query_TCGA_BRCA <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value", 
  platform = "Illumina Human Methylation 450"
) # Data corresponds to hg38 genome build

# Display a summary of the TCGA-BRCA project
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

# Download the queried data to the specified directory
GDCdownload(query = query_TCGA_BRCA, directory = "path/to/directory")

# Prepare the data by formatting it into a SummarizedExperiment object
tcga_data_BRCA <- GDCprepare(query_TCGA_BRCA, directory = "path/to/directory")

# Filter out rows with missing values across all assays
tcga_data_BRCA <- subset(tcga_data_BRCA, subset = (rowSums(is.na(assay(tcga_data_BRCA))) == 0))

# Keep only samples of type "Primary Tumor" or "Solid Tissue Normal"
tcga_data_BRCA <- tcga_data_BRCA[, tcga_data_BRCA$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]

# Set "Solid Tissue Normal" as the reference level for the sample_type factor
tcga_data_BRCA$sample_type <- relevel(factor(tcga_data_BRCA$sample_type), "Solid Tissue Normal")

# Generate summary tables of sample types, tissue/organ of origin, and primary diagnosis
# Save each summary table to a CSV file
da_table <- table(tcga_data_BRCA$sample_type)
write.csv(da_table, "HoxGenes/Sesame/Sesame/main/BRCA/BRCA_samples.csv")

da_table2 <- table(tcga_data_BRCA$tissue_or_organ_of_origin)
write.csv(da_table2, "HoxGenes/Sesame/Sesame/main/BRCA/BRCA_tissues.csv")

da_table3 <- table(tcga_data_BRCA$primary_diagnosis)
write.csv(da_table3, "HoxGenes/Sesame/Sesame/main/BRCA/BRCA_diagnosis.csv")

# Extract and save all metadata (column data) into a CSV file
whole_table <- colData(tcga_data_BRCA)
whole_table <- apply(whole_table, 2, as.character)
write.csv(whole_table, "HoxGenes/Sesame/Sesame/main/BRCA/All_BRCA.csv")

# Convert column data to a data frame for easier manipulation
col_data <- as.data.frame(colData(tcga_data_BRCA))
rownames(col_data) = NULL

# Redefine "Solid Tissue Normal" as the reference level for consistency
tcga_data_BRCA$sample_type <- relevel(factor(tcga_data_BRCA$sample_type), "Solid Tissue Normal")

# Perform differential methylation analysis (DML) based on sample type
smry <- DML(tcga_data_BRCA, ~sample_type)
write_rds(smry, "path/to/directory")

# Perform differential methylation region (DMR) analysis using the DML results
merged <- DMR(tcga_data_BRCA, smry, contrast = "sample_typePrimary.Tumor")
write_rds(merged, "/path/to/directory")

# Filter DMR results to include regions with adjusted p-value < 0.01
merged_2 <- merged %>% dplyr::filter(Seg_Pval_adj < 0.01)
write_rds(merged_2, "path/to/directory")
