# Load libraries
library(regioneR)
library(annotatr)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)

# Set working directory
setwd("path/to/directory")

# Define tissue list (these are all 16 used in this study)
tissue_list <- c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "THCA", "UCEC")

# Load data: this is the DMR data produced from the intial step, both the entire set and the set with the cut-off
merged_2 <- readRDS("path/to/rds")
universe <- readRDS("path/to/rds")

# Filter DMRs by segment ID and add methylation status
filt <- merged_2 %>%
  group_by(Seg_ID) %>%
  filter(n() > 1) %>%
  mutate(
    DM_status = Seg_Est,
    Seg_Est = ifelse(Seg_Est > 0, "hyper", "hypo")
  ) %>%
  distinct()

# Convert filtered data to GenomicRanges
filt_gr <- filt %>%
  select(Seg_ID, Seg_Chrm, Seg_Start, Seg_End, Seg_Pval_adj, Seg_Est, DM_status) %>%
  rename(Chr = Seg_Chrm, start = Seg_Start, end = Seg_End, pval = Seg_Pval_adj, Meth_status = Seg_Est) %>%
  unique() %>%
  as("GenomicRanges")

# Load annotations
annots <- c('hg38_cpgs', 'hg38_basicgenes', 'hg38_genes_1to5kb', 'hg38_genes_intergenic', 'hg38_genes_exons', 'hg38_genes_promoters')
annotations1 <- build_annotations(genome = "hg38", annotations = annots)

# Subset hypo- and hypermethylated DMRs
filt_gr_hypo <- subset(filt_gr, Meth_status == "hypo")
filt_gr_hyper <- subset(filt_gr, Meth_status == "hyper")

# Convert universe data to GenomicRanges
universe_gr <- universe %>%
  select(Seg_ID, Seg_Chrm, Seg_Start, Seg_End, Seg_Pval_adj, Seg_Est) %>%
  rename(Chr = Seg_Chrm, start = Seg_Start, end = Seg_End, pval = Seg_Pval_adj) %>%
  unique() %>%
  as("GenomicRanges")

# Load and filter homeobox genes
gencode <- import("../code/gencode.v44.annotation.gtf.gz")
genes <- subset(gencode, type == "gene")
homeo_genes <- read.csv("../code/HomeoboxList.csv")[, 2]
homeogenes_gr <- subset(genes, gene_name %in% homeo_genes)
homeogenes_gr2 <- resize(as(homeogenes_gr, "GenomicRanges"), width = 20000, fix = "center")

# Permutation tests for homeobox genes
perm_test <- function(A, B, universe) {
  permTest(A = A, B = B, universe = universe, ntime = 1000, randomize.function = resampleRegions, 
           evaluate.function = numOverlaps, alternative = "auto", verbose = TRUE)
}

pt_homeobox_total <- perm_test(filt_gr, homeogenes_gr2, universe_gr)
pt_homeobox_hypo <- perm_test(filt_gr_hypo, homeogenes_gr2, universe_gr)
pt_homeobox_hyper <- perm_test(filt_gr_hyper, homeogenes_gr2, universe_gr)

# Plot results
lapply(list(pt_homeobox_total, pt_homeobox_hypo, pt_homeobox_hyper), plot)

# Annotate DMRs
dm_annotated <- annotate_regions(filt_gr, annotations1, ignore.strand = TRUE, quiet = FALSE)
dm_annsum <- summarize_annotations(annotated_regions = dm_annotated)
dm_annsum$percent <- dm_annsum$n / sum(dm_annsum$n)

# Plot annotation summaries
ggplot(dm_annsum, aes(fill = annot.type, y = percent, x = "UCEC")) +
  geom_bar(position = "stack", stat = "identity")
