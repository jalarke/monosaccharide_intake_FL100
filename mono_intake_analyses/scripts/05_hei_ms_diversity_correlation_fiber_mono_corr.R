## Title: 05_hei_ms_diversity_correlation
## Author: Jules Larke
## Date: 041122
## Purpose: Correlation of HEI and monosaccharide diversity. Also includes correlations of fiber intake with total and non-glucose monosaccharide intake (supplemental figure 2)

library(tidyverse)
set.seed(1)
glycan_data <- read_csv('data/glycan_metadata_041122.csv') #Mono data units = g consumed per 1000 kcal averaged over all recalls
colnames(glycan_data)[1] = 'subject_id'
hei <- read_tsv('data/FL100_HEI_n378.txt')
hei <- hei %>% filter(subject_id %in% glycan_data$subject_id)

glycan <- glycan_data %>% select(subject_id, RecallNo, age, bmi, fecal_calprotectin, fecal_mpo, fecal_neopterin, fecal_ph, sex.factor, Glucose, Galactose, Fructose, Xylose, Arabinose, Fucose, Rhamnose, GalA, Mannose, Ribose, fiber_consumed_g)

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

library(vegan)
## calculate alpha diveristy
shannon_div <- vegan::diversity(glycan_total[,2:11], "shannon") 
nonglc_shannon <- vegan::diversity(glycan_total[,3:11], "shannon") 
glycan_total$shannon_div <- shannon_div
glycan_total$nonglc_shannon <- nonglc_shannon

mono_hei <- merge(glycan_total, hei, by = 'subject_id')
pheno <- glycan %>% select(subject_id, age, bmi, sex.factor)
pheno <- distinct(pheno, subject_id, .keep_all = TRUE)
mono_hei <- merge(mono_hei, pheno, by = 'subject_id')
#mono_hei <- mono_hei %>% mutate(value = 1) %>% spread(sex.factor, value,  fill = 0) 

mono_hei = mono_hei %>%
  mutate(nonglc_ms = select(., Galactose:Ribose) %>% rowSums(na.rm = TRUE)) %>%
  mutate(total_ms = select(., Glucose:Ribose) %>% rowSums(na.rm = TRUE))


mono_div_hei <- lm(shannon_div ~ hei_asa24_totalscore + age + bmi + sex.factor, data = mono_hei)
shapiro.test(mono_div_hei$residuals)

# Use bestNormalize
bn_hei_score <- bestNormalize::bestNormalize(mono_hei$hei_asa24_totalscore)
bn_shannon <- bestNormalize::bestNormalize(mono_hei$shannon_div)

mono_hei$bn_shannon <- bn_shannon$x.t
mono_hei$bn_hei_score <- bn_hei_score$x.t

mono_div_hei <- lm(bn_shannon ~ bn_hei_score + age + bmi + sex.factor, data = mono_hei)
summary(mono_div_hei)
plot(mono_div_hei, which = 2)
shapiro.test(mono_div_hei$residuals)

# Use log10 transformation
mono_div_hei <- lm(log10(shannon_div) ~ log10(hei_asa24_totalscore) + age + bmi + sex.factor, data = mono_hei)
summary(mono_div_hei)
plot(mono_div_hei, which = 2)
shapiro.test(mono_div_hei$residuals)

mono_hei$log_shann <- log10(mono_hei$shannon_div)
mono_hei$log_hei <- log10(mono_hei$hei_asa24_totalscore)

mono_hei <- mono_hei %>% mutate(value = 1) %>% spread(sex.factor, value,  fill = 0) 

pp <- ppcor::pcor.test(mono_hei$log_shann, mono_hei$log_hei, mono_hei[,c('age', 'bmi', 'Female')])
pp
round(pp$estimate, 3)
round(pp$p.value, 10)
plot(mono_hei$log_shann, mono_hei$log_hei)


r_label = paste0("italic(r) ==", 0.52)
p_label = paste0("italic(p) ==", 1.2e-13)

mono_hei_plot <- ggplot(mono_hei, aes(x = shannon_div, y = hei_asa24_totalscore)) +
  geom_point(color = "black",
             size = 2,
             shape = 16) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Shannon diversity of total monosaccharide intake') +
  labs(y = "Total HEI score") +
  scale_x_continuous(breaks = c(0.40, 0.80, 1.20)) +
  annotate(
    "text",
    x = 0.5,
    y = 93,
    label = r_label,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = 0.5,
    y = 88,
    label = p_label,
    parse = TRUE
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
mono_hei_plot

ggsave("output/correlation/mono_diversity_HEI_pearson_072122.tiff",
       plot = mono_hei_plot,
       width = 4,
       height = 4,
       dpi = 'retina')


#Correlations between fiber, total monosaccharide and non-glucose monosaccharide intake
nonglc_fiber <- cor.test(mono_hei$fiber_consumed_g, mono_hei$nonglc_ms)
totalms_fiber <- cor.test(mono_hei$fiber_consumed_g, mono_hei$total_ms)

r_label = paste0("italic(r) ==", 0.702)
p_label = paste0("italic(p) ==", 2.2e-16)
nonglc_fiber_plot <- ggplot(mono_hei, aes(x = fiber_consumed_g, y = nonglc_ms)) +
  geom_point(color = "black",
             size = 2,
             shape = 16) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Average fiber intake (g/day) per 1000 kcal') +
  labs(y = 'Average non-glucose monosaccharide intake\n(g/day) per 1000 kcal') +
  annotate(
    "text",
    x = 7,
    y = 24.5,
    label = r_label,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = 7,
    y = 23,
    label = p_label,
    parse = TRUE
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
nonglc_fiber_plot

r_label = paste0("italic(r) ==", 0.314)
p_label = paste0("italic(p) ==", 1.8e-5)
total_ms_fiber_plot <- ggplot(mono_hei, aes(x = fiber_consumed_g, y = total_ms)) +
  geom_point(color = "black",
             size = 2,
             shape = 16) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Average fiber intake (g/day) per 1000 kcal') +
  labs(y = 'Average total monosaccharide intake\n(g/day) per 1000 kcal') +
  annotate(
    "text",
    x = 25,
    y = 105,
    label = r_label,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = 25,
    y = 100,
    label = p_label,
    parse = TRUE
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
total_ms_fiber_plot

ms_fiber_corr <- cowplot::plot_grid(total_ms_fiber_plot, nonglc_fiber_plot, nrow = 1, align = 'hv', labels = 'AUTO')

ggsave("output/correlation/supp_fig_2_072122.tiff",
       plot = ms_fiber_corr,
       width = 8,
       height = 5,
       dpi = 1000)
