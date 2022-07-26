## Title: 06_gi_marker_ms_correlation
## Author: Jules Larke
## Date: 041122
## Purpose: Correlation of GI markers and monosaccharide intake/diversity

library(tidyverse)
library(bestNormalize)
set.seed(1)

glycan_data <- read_csv('data/glycan_metadata_041122.csv') #Mono data units = g consumed per 1000 kcal averaged over all recalls
colnames(glycan_data)[1] = 'subject_id'
fecal_data <- read.csv('data/FL100_stool_variables.txt', sep = '\t')
#hei <- hei %>% filter(subject_id %in% glycan_data$subject_id)

glycan <- glycan_data %>% select(subject_id, RecallNo, age, bmi, sex.factor, Glucose, Galactose, Fructose, Xylose, Arabinose, Fucose, Rhamnose, GalA, Mannose, Ribose, fiber_consumed_g)

glycan_total <- glycan %>% group_by(subject_id, RecallNo) %>%
  summarise(
    Glucose = sum(Glucose),
    Galactose = sum(Galactose),
    Fructose = sum(Fructose),
    Arabinose = sum(Arabinose),
    Xylose = sum(Xylose),
    Fucose = sum(Fucose),
    Rhamnose = sum(Rhamnose),
    GalA = sum(GalA),
    Mannose = sum(Mannose),
    Ribose = sum(Ribose),
    fiber_consumed_g = sum(fiber_consumed_g)
  ) %>%
  summarise(
    Glucose = mean(Glucose),
    Galactose = mean(Galactose),
    Fructose = mean(Fructose),
    Arabinose = mean(Arabinose),
    Xylose = mean(Xylose),
    Fucose = mean(Fucose),
    Rhamnose = mean(Rhamnose),
    GalA = mean(GalA),
    Mannose = mean(Mannose),
    Ribose = mean(Ribose),
    fiber_consumed_g = mean(fiber_consumed_g)
  )
glycan_total = glycan_total %>%
  mutate(nonglc_ms = select(., Galactose:Ribose) %>% rowSums(na.rm = TRUE)) %>%
  mutate(total_ms = select(., Glucose:Ribose) %>% rowSums(na.rm = TRUE))
glycan_all <- merge(glycan_total, fecal_data, by='subject_id')

library(vegan)
library(reshape2)
## calculate alpha diveristy for monosaccharides
shannon_div <- vegan::diversity(glycan_all[,2:11], "shannon") 
shannon_noglc <- vegan::diversity(glycan_all[,3:11], "shannon")
glycan_all$shannon_div <- shannon_div
glycan_all$shannon_noglc <- shannon_noglc
# Merge phenotype data with monosaccharide intake
pheno <- glycan %>% select(subject_id, age, bmi, sex.factor)
pheno <- distinct(pheno, subject_id, .keep_all = TRUE)
pheno <- pheno %>% mutate(value = 1) %>% spread(sex.factor, value,  fill = 0) 
mono_pheno <- merge(glycan_all, pheno, by='subject_id')

gi_pheno <- mono_pheno[complete.cases(mono_pheno[,28]),]
gi_pheno <- gi_pheno %>% mutate(value = 1) %>% spread(StoolConsistencyClass, value,  fill = 0, drop = TRUE) 

# Correlation between plasma LBP and monosaccharude intake / diversity
lbp <- mono_pheno %>% select(Glucose:total_ms, shannon_div, shannon_noglc, plasma_lbp_bd1)
MASS::truehist(lbp$plasma_lbp_bd1, nbins = 12)
MASS::truehist(mono_pheno$total_ms, nbins = 12)

bnlbp <- bestNormalize::boxcox(lbp$plasma_lbp_bd1)
bn_total_mono <- bestNormalize::boxcox(lbp$total_ms)
bn_nonglc_mono <- bestNormalize::boxcox(lbp$nonglc_ms)
bn_mono_div <- bestNormalize::boxcox(lbp$shannon_div)
bn_noglc_mono_div <- bestNormalize::boxcox(lbp$shannon_noglc)

lbp_noglc_div_res <- lm(bn_noglc_mono_div$x.t ~ bnlbp$x.t)
plot(lbp_noglc_div_res, which = 2)
shapiro.test(lbp_noglc_div_res$residuals)

lbp_div_res <- lm(bn_mono_div$x.t ~ bnlbp$x.t)
plot(lbp_div_res, which = 2)
shapiro.test(lbp_div_res$residuals)

lbp_total_res <- lm(bn_total_mono$x.t ~ bnlbp$x.t)
plot(lbp_total_res, which = 2)
shapiro.test(lbp_total_res$residuals)

lbp_nonglc_res <- lm(bn_nonglc_mono$x.t ~ bnlbp$x.t)
plot(lbp_nonglc_res, which = 2)
shapiro.test(lbp_nonglc_res$residuals)

mono_pheno$bn_total_mono <- bn_total_mono$x.t
mono_pheno$bn_nonglc_mono <- bn_nonglc_mono$x.t
mono_pheno$bn_mono_div <- bn_mono_div$x.t
mono_pheno$bnlbp <- bnlbp$x.t


ppcor::pcor.test(mono_pheno$bn_total_mono, mono_pheno$bnlbp, mono_pheno[,c("age", "bmi", "Female")])
ppcor::pcor.test(mono_pheno$bn_nonglc_mono, mono_pheno$bnlbp, mono_pheno[,c("age", "bmi", "Female")])
ppcor::pcor.test(mono_pheno$bn_mono_div, mono_pheno$bnlbp, mono_pheno[,c("age", "bmi", "Female")])


# Plot result:
r_label = paste0("italic(r) ==", 0.184)
p_label = paste0("italic(p) ==", 0.014)
lbp_total_plot <- ggplot(mono_pheno, aes(x = total_ms, y = plasma_lbp_bd1)) +
  geom_point(color = "black",
             size = 2,
             shape = 16, alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Monosaccharide intake\n(g/day per 1000 kcal)') +
  labs(y = "Plasma LBP (ug/g)") +
  annotate(
    "text",
    x = 38,
    y = 27.5,
    label = r_label,
    parse = TRUE,
    size=3
  ) +
  annotate(
    "text",
    x = 38,
    y = 25.5,
    label = p_label,
    parse = TRUE,
    size=3
  ) +
  theme_classic() +
  theme(
    axis.line.x = element_line(color = "black", size = .5),
    axis.line.y = element_line(color = "black", size = .5),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    ## legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank()
  )
lbp_total_plot

r_label = paste0("italic(r) ==", 0.079)
p_label = paste0("italic(p) ==", 0.3)
lbp_non_glc_plot <- ggplot(mono_pheno, aes(x = nonglc_ms, y = plasma_lbp_bd1)) +
  geom_point(color = "black",
             size = 2,
             shape = 16, alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Non-glucose monosaccharide\nintake (g/day per 1000 kcal)') +
  labs(y = "Plasma LBP (ug/g)") +
  annotate(
    "text",
    x = 10,
    y = 27.5,
    label = r_label,
    parse = TRUE,
    size =3 
  ) +
  annotate(
    "text",
    x = 10,
    y = 25.5,
    label = p_label,
    parse = TRUE,
    size =3
  ) +
  theme_classic() +
  theme(
    axis.line.x = element_line(color = "black", size = .5),
    axis.line.y = element_line(color = "black", size = .5),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    ## legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank()
  )
lbp_non_glc_plot


r_label = paste0("italic(r) ==", -0.11)
p_label = paste0("italic(p) ==", 0.15)
lbp_div_plot <- ggplot(mono_pheno, aes(x = shannon_div, y = plasma_lbp_bd1)) +
  geom_point(color = "black",
             size = 2,
             shape = 16, alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Monosaccharide intake diversity') +
  labs(y = "Plasma LBP (ug/g)") +
  annotate(
    "text",
    x = 0.65,
    y = 27.5,
    label = r_label,
    parse = TRUE,
    size=3
  ) +
  annotate(
    "text",
    x = 0.65,
    y = 25.5,
    label = p_label,
    parse = TRUE,
    size =3
  ) +
  theme_classic() +
  theme(
    axis.line.x = element_line(color = "black", size = .5),
    axis.line.y = element_line(color = "black", size = .5),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    ## legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank()
  )
lbp_div_plot

lbp_plot <- cowplot::plot_grid(lbp_total_plot, lbp_non_glc_plot, lbp_div_plot, nrow = 1, labels = c("D", "E", "F"), label_size = 8)
ggsave("output/correlation/lbp_correlaiton.tiff",
       plot = lbp_plot,
       width = 10,
       height = 4,
       dpi = 'retina')

## Correlation of gut-health outcomes and monsaccharides
# Fecal calprotectin
mono_fcal <- gi_pheno[complete.cases(gi_pheno[,19]),]
mono_fcal <- mono_fcal %>% filter(!AfterV2 == 1)
mono_fcal <- mono_fcal %>% filter(!diff_time_hrs > 24)
bn_fcal <- bestNormalize(mono_fcal$fecal_calprotectin)
MASS::truehist(bn_fcal$x.t, nbins = 12)
MASS::truehist(mono_fcal$shannon_div, nbins = 12)

ppcor::pcor.test(mono_fcal$total_ms, mono_fcal$fecal_calprotectin, mono_fcal[,c("age", "bmi", "Female", "hard", "normal")], method = 'spearman')
plot(mono_fcal$total_ms, mono_fcal$fecal_calprotectin)
ppcor::pcor.test(mono_fcal$nonglc_ms, mono_fcal$fecal_calprotectin, mono_fcal[,c("age", "bmi", "Female", "hard", "normal")], method = 'spearman')
plot(mono_fcal$nonglc_ms, mono_fcal$fecal_calprotectin)
ppcor::pcor.test(mono_fcal$shannon_div, mono_fcal$fecal_calprotectin, mono_fcal[,c("age", "bmi", "Female", "hard", "normal")], method = 'spearman')
plot(mono_fcal$shannon_div, mono_fcal$fecal_calprotectin)

# Fecal neopterin
mono_neop <- gi_pheno[complete.cases(gi_pheno[,17]),]
mono_neop <- mono_neop %>% filter(!AfterV2 == 1)
mono_neop <- mono_neop %>% filter(!diff_time_hrs > 24)
MASS::truehist(mono_neop$fecal_neopterin)

bn_neop <- bestNormalize::boxcox(mono_neop$fecal_neopterin)
bn_total_mono <- bestNormalize::boxcox(mono_neop$total_ms)
bn_nonglc_mono <- bestNormalize::boxcox(mono_neop$nonglc_ms)
bn_mono_div <- bestNormalize::boxcox(mono_neop$shannon_div)

mono_neop$bn_neop <- bn_neop$x.t
mono_neop$bn_total_mono <- bn_total_mono$x.t
mono_neop$bn_nonglc_mono <- bn_nonglc_mono$x.t
mono_neop$bn_mono_div <- bn_mono_div$x.t

neop_total_res <- lm(mono_neop$bn_total_mono ~ mono_neop$bn_neop)
plot(neop_total_res, which = 2)
shapiro.test(neop_total_res$residuals)

neop_nonglc_res <- lm(mono_neop$bn_nonglc_mono ~ mono_neop$bn_neop)
plot(neop_nonglc_res, which = 2)
shapiro.test(neop_nonglc_res$residuals)

neop_div_res <- lm(mono_neop$bn_mono_div ~ mono_neop$bn_neop)
plot(neop_div_res, which = 2)
shapiro.test(neop_div_res$residuals)

ppcor::pcor.test(mono_neop$bn_total_mono, mono_neop$bn_neop, mono_neop[,c("age", "bmi", "Female", "hard", "normal")])
ppcor::pcor.test(mono_neop$bn_nonglc_mono, mono_neop$bn_neop, mono_neop[,c("age", "bmi", "Female", "hard", "normal")])
ppcor::pcor.test(mono_neop$bn_mono_div, mono_neop$bn_neop, mono_neop[,c("age", "bmi", "Female", "hard", "normal")])

# Plot result:
r_label = paste0("italic(r) ==", 0.147)
p_label = paste0("italic(p) ==", 0.076)
neop_total_plot <- ggplot(mono_neop, aes(x = total_ms, y = fecal_neopterin)) +
  geom_point(color = "black",
             size = 2,
             shape = 16, alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Monosaccharide intake\n(g/day per 1000 kcal)') +
  labs(y = "Fecal neopterin (ug/g)") +
  annotate(
    "text",
    x = 40,
    y = 205,
    label = r_label,
    parse = TRUE,
    size=3
  ) +
  annotate(
    "text",
    x = 40,
    y = 190,
    label = p_label,
    parse = TRUE,
    size=3
  ) +
  theme_classic() +
  theme(
    axis.line.x = element_line(color = "black", size = .5),
    axis.line.y = element_line(color = "black", size = .5),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    ## legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank()
  )
neop_total_plot

r_label = paste0("italic(r) ==", -0.09)
p_label = paste0("italic(p) ==", 0.28)
neop_non_glc_plot <- ggplot(mono_neop, aes(x = nonglc_ms, y = fecal_neopterin)) +
  geom_point(color = "black",
             size = 2,
             shape = 16, alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Non-glucose monosaccharide\nintake (g/day per 1000 kcal)') +
  labs(y = "Fecal neopterin (ug/g)") +
  annotate(
    "text",
    x = 10,
    y = 205,
    label = r_label,
    parse = TRUE,
    size=3
  ) +
  annotate(
    "text",
    x = 10,
    y = 190,
    label = p_label,
    parse = TRUE,
    size=3
  ) +
  theme_classic() +
  theme(
    axis.line.x = element_line(color = "black", size = .5),
    axis.line.y = element_line(color = "black", size = .5),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    ## legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank()
  )
neop_non_glc_plot


r_label = paste0("italic(r) ==", -0.247)
p_label = paste0("italic(p) ==", 0.003)
neop_div_plot <- ggplot(mono_neop, aes(x = shannon_div, y = fecal_neopterin)) +
  geom_point(color = "black",
             size = 2,
             shape = 16, alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Monosaccharide intake diversity') +
  labs(y = "Fecal neopterin (ug/g)") +
  annotate(
    "text",
    x = 0.65,
    y = 205,
    label = r_label,
    parse = TRUE,
    size=3
  ) +
  annotate(
    "text",
    x = 0.65,
    y = 190,
    label = p_label,
    parse = TRUE,
    size=3
  ) +
  theme_classic() +
  theme(
    axis.line.x = element_line(color = "black", size = .5),
    axis.line.y = element_line(color = "black", size = .5),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    ## legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank()
  )
neop_div_plot

neop_plot <- cowplot::plot_grid(neop_total_plot, neop_non_glc_plot, neop_div_plot, nrow = 1, labels = "AUTO", label_size = 8)

figure_8 <- cowplot::plot_grid(neop_plot, lbp_plot, nrow = 2, align = 'hv')

ggsave("output/correlation/figure_8_072122.tiff",
       plot = figure_8,
       width = 5.7,
       height = 5,
       dpi = 1000)

# Fecal MPO
mono_mpo <- gi_pheno[complete.cases(gi_pheno[,20]),]
mono_mpo <- mono_mpo %>% filter(!AfterV2 == 1)
mono_mpo <- mono_mpo %>% filter(!diff_time_hrs > 24)
MASS::truehist(mono_mpo$fecal_mpo)
bn_mpo <- bestNormalize(mono_mpo$fecal_mpo)
MASS::truehist(bn_mpo$x.t)

ppcor::pcor.test(mono_mpo$total_ms, bn_mpo$x.t, mono_mpo[,c("age", "bmi", "Female", "hard", "normal")])
ppcor::pcor.test(mono_mpo$nonglc_ms, bn_mpo$x.t, mono_mpo[,c("age", "bmi", "Female", "hard", "normal")])
ppcor::pcor.test(mono_mpo$shannon_div, bn_mpo$x.t, mono_mpo[,c("age", "bmi", "Female", "hard", "normal")])


# Fecal pH
mono_ph <- gi_pheno[complete.cases(gi_pheno[,23]),]
mono_ph <- mono_ph %>% filter(!AfterV2 == 1)
mono_ph <- mono_ph %>% filter(!diff_time_hrs > 24)
#mono_ph <- mono_ph %>% filter(!fecal_ph > 10000)
MASS::truehist(mono_ph$fecal_ph)
bn_ph <- bestNormalize(mono_ph$fecal_ph)
MASS::truehist(bn_ph$x.t)

ppcor::pcor.test(mono_ph$total_ms, mono_ph$fecal_ph, mono_ph[,c("age", "bmi", "Female", "hard", "normal")])
plot(mono_ph$total_ms, mono_ph$fecal_ph)
ppcor::pcor.test(mono_ph$nonglc_ms, mono_ph$fecal_ph, mono_ph[,c("age", "bmi", "Female", "hard", "normal")])
plot(mono_ph$nonglc_ms, mono_ph$fecal_ph)
ppcor::pcor.test(mono_ph$shannon_div, mono_ph$fecal_ph, mono_ph[,c("age", "bmi", "Female", "hard", "normal")])
plot(mono_ph$shannon_div, mono_ph$fecal_ph)
