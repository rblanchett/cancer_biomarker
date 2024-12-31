# Load necessary libraries
library(regioneR)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(stringr)

# Set working directory
setwd("C:/Users/reid.blanchett/OneDrive - Van Andel Institute/Sesame/main")

# Define tissue list
tissue_list <- c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "HNSC", "KIRC", "KIRP", 
                 "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "THCA", "BRCA")

# Load Gencode annotation and filter for genes
gencode <- import("gencode.v44.annotation.gtf.gz")
genes <- gencode[gencode$type == "gene"]

# Load homeobox gene list and subset Gencode data
homeo_genes <- read.csv("HomeoboxList.csv")
homeogenes_gr <- genes[genes$gene_name %in% homeo_genes[[2]]]
homeogenes_df <- as.data.frame(homeogenes_gr)[, 1:3]
homeogenes_df$start <- homeogenes_df$start - 10000
homeogenes_df$end <- homeogenes_df$end + 10000
homeogenes_gr <- as(homeogenes_df, "GenomicRanges")

# Load and process H3K27me3 data
H3k27_bed <- read.csv("Polycomb/polycomb_bed_hg38.csv", header = FALSE)
H3k27_gr <- H3k27_bed %>%
  dplyr::select(V1, V2, V3) %>%
  dplyr::rename(Chr = V1, start = V2, end = V3) %>%
  unique() %>%
  as("GenomicRanges")

# Initialize data frames for results
PRC2_df <- data.frame()
homeo_df <- data.frame()

# Loop through tissues for analysis
for (tissue in tissue_list) {
  # Load data
  merged <- readRDS(paste0(tissue, "/", tissue, "_beta_merged_reverse.rds"))
  universe <- readRDS(paste0(tissue, "/", tissue, "_beta_full_reverse.rds"))
  
  # Filter and prepare data
  filt <- merged %>%
    group_by(Seg_ID) %>%
    filter(n() > 1) %>%
    select(Seg_ID, Seg_Chrm, Seg_Start, Seg_End, Seg_Pval_adj) %>%
    rename(Chr = Seg_Chrm, start = Seg_Start, end = Seg_End, pval = Seg_Pval_adj) %>%
    unique()
  
  universe <- universe %>%
    select(Seg_ID, Seg_Chrm, Seg_Start, Seg_End) %>%
    rename(Chr = Seg_Chrm, start = Seg_Start, end = Seg_End) %>%
    unique()
  
  filt_gr <- as(filt, "GenomicRanges")
  universe_gr <- as(universe, "GenomicRanges")
  
  # Identify overlaps and perform permutation tests
  DMR_and_homeobox <- subsetByOverlaps(filt_gr, homeogenes_gr)
  
  pt_PRC2 <- permTest(A = DMR_and_homeobox, B = H3k27_gr, universe = universe_gr,
                      ntime = 1000, randomize.function = resampleRegions,
                      evaluate.function = numOverlaps, alternative = "auto", verbose = TRUE)
  plot(pt_PRC2)
  PRC2_df <- rbind(PRC2_df, data.frame(Tissue = tissue, 
                                       pval = pt_PRC2$numOverlaps$pval, 
                                       Z_score = pt_PRC2$numOverlaps$zscore, 
                                       Alternative = pt_PRC2$numOverlaps$alternative))
  
  pt_homeo <- permTest(A = filt_gr, B = homeogenes_gr, universe = universe_gr,
                       ntime = 1000, randomize.function = resampleRegions,
                       evaluate.function = numOverlaps, alternative = "auto", verbose = TRUE)
  plot(pt_homeo)
  homeo_df <- rbind(homeo_df, data.frame(Tissue = tissue, 
                                         pval = pt_homeo$numOverlaps$pval, 
                                         Z_score = pt_homeo$numOverlaps$zscore, 
                                         Alternative = pt_homeo$numOverlaps$alternative))
}

# Save results
write.csv(PRC2_df, "Permutation_PRC2.csv", row.names = FALSE)
write.csv(homeo_df, "Permutation_homeobox.csv", row.names = FALSE)
