# Set the working directory
setwd("path/to/directory")

# Load required libraries
library(rtracklayer)
library(dplyr)
library(stringr)
library(GenomicRanges)

# Import GENCODE gene annotation file and filter for genes
gencode <- import("gencode.v44.annotation.gtf.gz")
genes <- gencode[gencode$type == "gene"]

# Load the list of Homeobox genes from a file
homeo_genes <- read.csv("HomeoboxList.csv")
homeogenes_gr <- genes[genes$gene_name %in% homeo_genes[[2]]]

# Load BRCA-specific DMR data
BRCA_beta <- readRDS("path/to/file")

# Filter DMRs to include only regions defined by multiple CpG sites
filt <- BRCA_beta %>%
  group_by(Seg_ID) %>%
  filter(n() > 1)

# Prepare DMR data for overlap analysis
filt_df <- filt %>%
  select(Seg_ID, Seg_Chrm, Seg_Start, Seg_End) %>%
  rename(Chr = Seg_Chrm, start = Seg_Start, end = Seg_End)
filt_gr <- as(filt_df, "GenomicRanges")

# Identify overlaps between Homeobox genes and DMRs
homeobox_final <- subsetByOverlaps(filt_gr, homeogenes_gr)
homeobox_final_df <- as.data.frame(homeobox_final)
h_distinct <- n_distinct(homeobox_final_df[[6]])

# Filter for hypermethylated DMRs
filt_hyper <- merged.2 %>%
  group_by(Seg_ID) %>%
  filter(n() > 1, Seg_Est > 0)
filt_df_hyper <- filt_hyper %>%
  select(Seg_ID, Seg_Chrm, Seg_Start, Seg_End) %>%
  rename(Chr = Seg_Chrm, start = Seg_Start, end = Seg_End)
filt_gr_hyper <- as(filt_df_hyper, "GenomicRanges")

# Overlap hypermethylated DMRs with Homeobox genes
homeobox_final_hyper <- subsetByOverlaps(filt_gr_hyper, homeogenes_gr)
homeobox_final_df_hyper <- as.data.frame(homeobox_final_hyper)
hh_distinct <- n_distinct(homeobox_final_df_hyper[[6]])

# Load Polycomb data and prepare for overlap analysis
full_bed <- read.csv("Polycomb/hg38_lift_polycomb_full.csv")
bed_df <- full_bed %>%
  select(Column1, Column2, Column3) %>%
  rename(Chr = Column1, start = Column2, end = Column3)
bed_gr <- as(bed_df, "GenomicRanges")

# Identify overlaps between DMRs and Polycomb regions
polycomb_final <- subsetByOverlaps(filt_gr, bed_gr)
filt_poly <- polycomb_final %>%
  group_by(Seg_ID) %>%
  filter(n() > 1)
p_distinct <- n_distinct(filt_poly[[6]])

# Overlap hypermethylated DMRs with Polycomb regions
polycomb_final_hyper <- subsetByOverlaps(filt_gr_hyper, bed_gr)
filt_poly_hyper <- polycomb_final_hyper %>%
  group_by(Seg_ID) %>%
  filter(n() > 1)
ph_distinct <- n_distinct(filt_poly_hyper[[6]])

# Identify overlaps between Homeobox genes and Polycomb regions in hypermethylated DMRs
overlap <- subsetByOverlaps(homeobox_final, subsetByOverlaps(filt_gr, bed_gr))
o_distinct <- n_distinct(overlap[[6]])

# Total number of distinct DMRs
total_distinct <- n_distinct(merged.2[[1]])

# Output summary statistics
print(c(h_distinct, hh_distinct, p_distinct, ph_distinct, o_distinct, total_distinct))
