setwd("path/to/directory")
library(dplyr)

#############
### Homeobox section
#############

### This is the list of DMRs of the cancer type in total
merged_2 <- readRDS("path/to/rds")
### This is the list of probes of the homeobox genes
Probe_List <- readRDS("Probes_by_Gene.rds")
probe_list_u <- unlist(Probe_List)
#Here we are winnowing down the large DMRs from the cancer type to just the ones
# with probes in the genes we are concerned about from the homeobox list
tmp_h <- merged_2[merged_2[[8]] %in% probe_list_u,]

tmp_h <- filter(merged_2, merged_2[[8]] %in% probe_list_u)
#Here we are taking the list and narrowing it down to the DMRs with 2 or more probes in it
final_h <- tmp_h %>% 
  group_by(Seg_ID) %>%
  filter(n() >1)
# This gives us the number of distinct DMRs based on the Segment ID
# It accounts for the duplicates which were throwing off the calculations of how 
# many DMRs there actually were
tmp2_h <- distinct(final_h[1])
# This calculation is the same for both because it comes from the original DMR 
# in the specific cancer
total_DMR <- merged_2 %>% group_by(Seg_ID) %>% filter(n()>1)
total_DMR1 <- distinct(total_DMR[1])

#####################
### Polycomb section
#####################
# This is the list of probes associated with each peak from the ChIP-seq data of H3k27me3 and polycomb complex
peaks_probes <- readRDS("C:\\Users\\reid.blanchett\\OneDrive - Van Andel Institute\\Sesame\\main\\Polycomb\\peak_probes.rds")
peaks_probes_u <- unlist(peaks_probes)
# Same idea as above. Finding the probes that intersect and keeping the ones
# in the larger file that do
tmp_p <- merged_2[merged_2[[8]] %in% peaks_probes_u,]
# grouping by segment ID and then choosing which DMRs have 2 probes or more
final_p <- tmp_p %>% 
  group_by(Seg_ID) %>%
  filter(n() >1)
# determining how many there are in total that are unique
tmp2_p <- distinct(final_p[1])

###################
### Overlap of homeobox genes with polycomb mark
###################

tmp_o <- tmp_h[tmp_p$Probe_ID %in% tmp_h$Probe_ID,]

final_o <- tmp_o %>% group_by(Seg_ID) %>% filter(n() > 1)

tmp2_o <- distinct(final_o[1])

print(dim(tmp2_p)[1])
print(dim(tmp2_h)[1])
print(dim(total_DMR1)[1])
print(dim(tmp2_o)[1])

##################
### find number of homeobox genes that cross over with only the hypermethylated cancer DMRs
##################

tmp_h_hyper <- final_h %>% filter(Seg_Est > 0)
tmp2_h_hyper <- distinct(tmp_h_hyper[1])

print(dim(tmp2_h_hyper)[1])

#################
### find number of polycomb genes that cross over with only the hypermethylated cancer DMRs
#################

tmp_p_hyper <- final_p %>% filter(Seg_Est > 0)
tmp2_p_hyper <- distinct(tmp_p_hyper[1])

print(dim(tmp2_p_hyper)[1])

########################
### find the number of stat. significant DMRs from the cancer tissue that are hypermethylated
########################

total_DMR <- merged_2 %>% group_by(Seg_ID) %>% filter(n()>1)
total_DMR_hyper <- total_DMR %>% filter(Seg_Est > 0)
total_DMR1_hyper <- distinct(total_DMR_hyper[1])

print(dim(total_DMR1_hyper)[1])

########################
### find the number of hypermethylated DMRs in polycomb and homeobox
### that are also in the hypermthylated DMR list for the cancer tissue
########################

#final_hh <- filter(tmp_h_hyper, tmp_h_hyper[[8]] %in% total_DMR_hyper[8,])
final_hh <- total_DMR_hyper[total_DMR_hyper[[8]] %in% tmp_h_hyper[[8]],]
distinct_final_hh <- distinct(final_hh[1])

