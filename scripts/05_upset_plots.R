# Set working directory
setwd("C:/Users/reid.blanchett/OneDrive - Van Andel Institute/Sesame/main")

# Required Libraries
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(stringr)

# Import GTF and filter for genes
gencode <- import("gencode.v44.annotation.gtf.gz")
genes <- gencode[gencode$type == "gene"]

# Load and filter homeobox genes
homeo_genes <- read.csv("HomeoboxList.csv")
homeogenes_gr <- genes[genes$gene_name %in% homeo_genes[[2]]]
homeogenes_df <- as.data.frame(homeogenes_gr)[, 1:3]
homeogenes_df$start <- homeogenes_df$start - 10000
homeogenes_df$end <- homeogenes_df$end + 10000
homeogenes_gr2 <- as(homeogenes_df, "GenomicRanges")

# Load DMR data
merged.2 <- readRDS("UCEC/UCEC_beta_merged_reverse.rds")

# Filter DMRs
filt <- merged.2 %>%
  group_by(Seg_ID) %>%
  filter(n() > 1)

filt_extra <- filt %>%
  filter(abs(Estimate) > 0.2)

filt_df <- filt_extra %>%
  select(Seg_ID, Seg_Chrm, Seg_Start, Seg_End) %>%
  rename(Chr = Seg_Chrm, start = Seg_Start, end = Seg_End)

filt_gr <- as(filt_df, "GenomicRanges")

# Subset by overlaps with homeobox genes
homeobox_final <- subsetByOverlaps(filt_gr, homeogenes_gr2)
homeobox_final_df <- as.data.frame(homeobox_final)

# Count distinct overlaps
distinct_count <- dim(distinct(homeobox_final_df[6]))[1]

# Hypermethylated homeobox analysis
filt_extra_hyper <- filt_extra %>%
  filter(Seg_Est > 0, abs(Estimate) > 0.2)

filt_df_hyper <- filt_extra_hyper %>%
  select(Seg_ID, Seg_Chrm, Seg_Start, Seg_End) %>%
  rename(Chr = Seg_Chrm, start = Seg_Start, end = Seg_End)

filt_gr_hyper <- as(filt_df_hyper, "GenomicRanges")
homeobox_final_hyper <- subsetByOverlaps(filt_gr_hyper, homeogenes_gr2)
homeobox_final_df_hyper <- as.data.frame(homeobox_final_hyper)
hyper_distinct_count <- dim(distinct(homeobox_final_df_hyper[6]))[1]

# Combine hypermethylation results
combo <- subsetByOverlaps(as(homeogenes_gr2, "GenomicRanges"), filt_gr_hyper)
combo_df <- as.data.frame(combo)

# Hypomethylated homeobox analysis
filt_extra_hypo <- filt_extra %>%
  filter(Seg_Est < 0, abs(Estimate) > 0.2)

filt_df_hypo <- filt_extra_hypo %>%
  select(Seg_ID, Seg_Chrm, Seg_Start, Seg_End) %>%
  rename(Chr = Seg_Chrm, start = Seg_Start, end = Seg_End)

filt_gr_hypo <- as(filt_df_hypo, "GenomicRanges")
homeobox_final_hypo <- subsetByOverlaps(filt_gr_hypo, homeogenes_gr2)
homeobox_final_df_hypo <- as.data.frame(homeobox_final_hypo)
hypo_distinct_count <- dim(distinct(homeobox_final_df_hypo[6]))[1]

# Load and process tissue-specific data
tissue_files <- list.files(pattern = "homeobox_list_represented.*\\.csv")
tissues <- lapply(tissue_files, read.csv, header = TRUE)

# Rename and combine tissue data
tissues <- lapply(tissues, function(df) {
  colnames(df) <- c("delete", "Cancer", "Homeobox_Gene")
  df <- df[,-1]
  return(df)
})
combined_tissues <- do.call(rbind, tissues)

# Group by homeobox gene for final results
combined_tissues_summary <- combined_tissues %>%
  group_by(Homeobox_Gene)

# Generate upset plot (optional)
# library(ComplexUpset)
# upset(...)

# Final results for hypomethylation and hypermethylation
list(
  HyperMethylatedDistinctCount = hyper_distinct_count,
  HypoMethylatedDistinctCount = hypo_distinct_count
)
