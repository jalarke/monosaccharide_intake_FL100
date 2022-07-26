## Title: 00_merge_metadata
## Author: Jules Larke
## Date: 010822
## Purpose: Integrate food recall monosaccharide data with subject metadata, physio and GI markers


library(tidyverse)

# Set file path
getwd()
wd = list()
wd$data = "/Users/jules.larke/work/project/glycan_library/mapping_foods/data/final/"

# Load datasets
mono <- read.csv(file.path(wd$data, 'all_items_cal_adjusted_041122.csv'))
colnames(mono)[1] = 'subject_id'
physio = read.csv(file.path(wd$data, 'meta_clinical_data/physiology_vitals.csv'))
physio <- physio %>% select(subject_id, bodyfat_pcnt, age, bmi, sex.factor)
ethnic = read.csv(file.path(wd$data, 'meta_clinical_data/DEXA_ethnicities04272020.csv'))
gi_markers= read.csv(file.path(wd$data, 'meta_clinical_data/gi_markers.csv'))
gi_markers <- gi_markers %>% select(subject_id, fecal_calprotectin, fecal_mpo, fecal_neopterin, fecal_ph)
# Merge glycan_microbe and ethnicity data
gly_eth = left_join(mono, ethnic, by = "subject_id")

sum(is.na(gly_eth$Ethnicity)) # 271 rows corresponding to 3 subjects have missing ethnicity data

gly_meta = left_join(gly_eth, physio, by = "subject_id")

gly_all = left_join(gly_meta, gi_markers, by = "subject_id")

gly_all$idx <- NULL

write_csv(gly_all, "data/glycan_metadata_041122.csv")
