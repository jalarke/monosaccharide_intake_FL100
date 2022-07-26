## Alpha and beta diversity of microbiome data from Kraken
## Jules Larke
## 042122
## Purpose: create plots for alpha/beta diversity across monosaccharide intake metrics (supplemental figure 3)

library(tidyverse)
library(vegan)
set.seed(1)
#Import data
metagenome_285 <- read.csv('data/metagenome/OTU_kraken.csv')
colnames(metagenome_285)[8:297] <- substr(colnames(metagenome_285)[8:297], 2,5)
glycan_data <- read.csv("data/glycan_metadata_041122.csv")
colnames(glycan_data)[1] = 'subject_id'
metagenome_285 <- metagenome_285 %>% select(Genus, `5001`:`9066`)
meta_genus <- metagenome_285 %>% group_by(Genus) %>% summarise(across(everything(), .fns = sum))
meta_genus <- as.data.frame(meta_genus)

meta <- setNames(data.frame(t(meta_genus[,-1])), meta_genus[,1])
meta <- meta %>% rownames_to_column('subject_id')

physio <- read_csv('../../mapping_foods/data/final/meta_clinical_data/physiology_vitals.csv')
physio <- physio %>% select(subject_id, age, bmi, sex.factor)
colnames(glycan_data)[1] = 'subject_id'
fecal_data <- read.csv('data/FL100_stool_variables.txt', sep = '\t')

# Sum total monosaccharides per participant
glycan_total <- glycan_data %>% group_by(subject_id, RecallNo) %>%
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
    fiber_consumed_g = sum(fiber_consumed_g),
    carb_consumed_g = sum(carb_consumed_g)
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
    fiber_consumed_g = mean(fiber_consumed_g),
    carb_consumed_g = mean(carb_consumed_g))
glycan_total = glycan_total %>%
  mutate(nonglc_ms = select(., Galactose:Ribose) %>% rowSums(na.rm = TRUE)) %>% mutate(total_ms = select(., Glucose:Ribose) %>% rowSums(na.rm = TRUE))

meta <- meta %>% filter(subject_id %in% glycan_total$subject_id)
physio <- physio %>% filter(subject_id %in% glycan_total$subject_id)

glycan_total <- merge(glycan_total, physio, by='subject_id')
glycan_all <- merge(glycan_total, fecal_data, by = 'subject_id')
glycan_all <- glycan_all %>% filter(subject_id %in% meta$subject_id)

taxa_mono <- merge(meta, glycan_total, by='subject_id')
taxa_mono <- taxa_mono %>% mutate(value = 1) %>% spread(sex.factor, value,  fill = 0) 


#Calculate alpha diversity
taxa_mono$shannon_div <- diversity(taxa_mono[,2:2284], "shannon") 
taxa_mono$richness <- specnumber(taxa_mono[,2:2284])
taxa_mono$evenness <- taxa_mono$shannon_div/log(taxa_mono$richness)
taxa_mono$nonglc_shannon <- diversity(taxa_mono[,2286:2294])
taxa_mono$mono_shannon <- diversity(taxa_mono[,2285:2294])

#Alpha diversity correlations
obs_nonglc_shannon <- lm(taxa_mono$nonglc_shannon ~ taxa_mono$richness)

summary(obs_nonglc_shannon)
plot(obs_nonglc_shannon, which = 2)
shapiro.test(obs_nonglc_shannon$residuals)

bn_observed <- bestNormalize::bestNormalize(taxa_mono$richness)
bn_taxa_shannon <- bestNormalize::boxcox(taxa_mono$shannon_div)
bn_mono_shannon <- bestNormalize::bestNormalize(taxa_mono$mono_shannon)
bn_nonglc_shannon <- bestNormalize::boxcox(taxa_mono$nonglc_shannon)

obs_nonglc_shannon <- lm(bn_nonglc_shannon$x.t ~ bn_observed$x.t)
summary(obs_nonglc_shannon)
plot(obs_nonglc_shannon, which = 2)
shapiro.test(obs_nonglc_shannon$residuals)

shan_nonglc_shannon <- lm(taxa_mono$nonglc_shannon ~ taxa_mono$shannon_div)
summary(shan_nonglc_shannon)
plot(shan_nonglc_shannon, which = 2)
shapiro.test(shan_nonglc_shannon$residuals)

shan_nonglc_shannon <- lm(bn_taxa_shannon$x.t ~ bn_nonglc_shannon$x.t)
summary(shan_nonglc_shannon)
plot(shan_nonglc_shannon, which = 2)
shapiro.test(shan_nonglc_shannon$residuals)

shan_mono_shannon <- lm(taxa_mono$mono_shannon ~ taxa_mono$shannon_div)
summary(shan_mono_shannon)
plot(shan_mono_shannon, which = 2)
shapiro.test(shan_mono_shannon$residuals)

shan_mono_shannon <- lm(bn_mono_shannon$x.t ~ bn_taxa_shannon$x.t)
summary(shan_mono_shannon)
plot(shan_mono_shannon, which = 2)
shapiro.test(shan_mono_shannon$residuals)

ppcor::pcor.test(x = bn_observed$x.t, y = bn_mono_shannon$x.t, z = taxa_mono[c('age', 'bmi', 'Female')])
ppcor::pcor.test(x = bn_observed$x.t, y = bn_nonglc_shannon$x.t, z = taxa_mono[c('age', 'bmi', 'Female')])
ppcor::pcor.test(x = bn_taxa_shannon$x.t, y = bn_mono_shannon$x.t, z = taxa_mono[c('age', 'bmi', 'Female')])
ppcor::pcor.test(x = bn_taxa_shannon$x.t, y = bn_nonglc_shannon$x.t, z = taxa_mono[c('age', 'bmi', 'Female')])

# Plot result:
r_label = paste0("italic(r) ==", 0.137)
p_label = paste0("italic(p) ==", 0.096)
c1 <- ggplot(taxa_mono, aes(x = nonglc_shannon, y = shannon_div)) +
  geom_point(color = "black",
             size = 2,
             shape = 16) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Shannon diversity of total\nnon-glucose monosaccharide intake') +
  labs(y = "Shannon diversity of microbiotas") +
  annotate(
    "text",
    x = 1.1,
    y = 4.22,
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
c1

r_label = paste0("italic(r) ==", 0.205)
p_label = paste0("italic(p) ==", 0.012)
c2 <- ggplot(taxa_mono, aes(x = nonglc_shannon, y = richness)) +
  geom_point(color = "black",
             size = 2,
             shape = 16) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Shannon diversity of total\nnon-glucose monosaccharide intake') +
  labs(y = "Observed taxa") +
  scale_x_continuous(breaks = c(1.0, 1.4, 1.8)) +
  annotate(
    "text",
    x = 1.1,
    y = 800,
    label = r_label,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 770,
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
c2

p_label = paste0("italic(p) ==", 0.99)
c3 <- ggplot(taxa_mono, aes(x = mono_shannon, y = shannon_div)) +
  geom_point(color = "black",
             size = 2,
             shape = 16) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Shannon diversity of total\nmonosaccharide intake') +
  labs(y = "Shannon diversity of microbiotas") +
  annotate(
    "text",
    x = 1.1,
    y = 4.2,
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
c3

figure_5 <- cowplot::plot_grid(mono_hei_plot, c2, nrow = 1, align = 'hv', labels = "AUTO", label_size = 8) # mono_hei_plot from 05_hei_ms_diversity_correlation

ggsave("output/correlation/figure_5_072122.tiff",
       plot = figure_5,
       width = 5.7,
       height = 5.7/2,
       dpi = 1000)

supp_fig_3 = cowplot::plot_grid(c3,c1, nrow = 1, labels = 'AUTO', label_size = 8)
ggsave("output/correlation/supp_fig_3AB_072122.tiff",
       plot = supp_fig_3,
       width = 5.7,
       height = 5.7/2,
       dpi = 1000)

#Beta diversity results
library(phyloseq)
library(microbiome)

source('scripts/08_genus_output_for_deseq.R')
#Age: Use study bins "18 - 33", "34 - 49", "50 - 65"
meta.a <- mutate(physio, Age_Category = "34_49")
meta.a$Age_Category[as.numeric(meta.a$age) <34] <- "18_33"
meta.a$Age_Category[as.numeric(meta.a$age) >49] <- "50_65"

#BMI bins
meta.b <- mutate(meta.a, BMI_Category = "25_to_29.9")
meta.b$BMI_Category[as.numeric(meta.b$bmi) < 25 ] <- "less_than_25"
meta.b$BMI_Category[as.numeric(meta.b$bmi) > 30 ] <- "30_to_41"

#Sequence count data
taxt <- meta_genus[,1] %>%
  as.matrix()
colnames(taxt) <- 'Genus'
rownames(taxt) <- taxt[,1]

meta_genus <- meta_genus %>% column_to_rownames(var = 'Genus')

#Build a phyloseq object with the data
physeq<-phyloseq(
  otu_table(meta_genus, taxa_are_rows = T), 
  # phy_tree(tree$data), 
  tax_table(taxt),
  sample_data(meta.b %>% as.data.frame() %>% column_to_rownames("subject_id"))
)

dim(otu_table(physeq))

physeq <- subset_samples(physeq, physeq@sam_data@row.names %in% glycan_total$subject_id)

physeq@sam_data$ara_quartile <- factor(ntile(glycan_all$Arabinose, 4))
physeq@sam_data$xyl_quartile <- factor(ntile(glycan_all$Xylose, 4))
physeq@sam_data$gala_quartile <- factor(ntile(glycan_all$GalA, 4))
physeq@sam_data$mann_quartile <- factor(ntile(glycan_all$Mannose, 4))
physeq@sam_data$total_ms_quartile <- factor(ntile(glycan_all$total_ms, 4))
physeq@sam_data$nonglc_quartile <- factor(ntile(glycan_all$nonglc_ms, 4))

physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq) # remove taxa not present in any samples

physeq <- microbiome::transform(physeq, "compositional") # covert to relative abundance


## Calculate distance metrics in samples of top/bottom quartiles of intake
set.seed(1)
bray_dist <- phyloseq::distance(physeq, method = "bray")
nmds <- vegan::metaMDS(bray_dist, trymax = 100, k = 3)
nmds

library(vegan)
##Total MS
#Beta dispersion test
nmdsdisp <- vegdist(bray_dist)
disper.all<-betadisper(nmdsdisp, physeq@sam_data$total_ms_quartile)
TukeyHSD(disper.all) 

#Permanova
adtest <- adonis(bray_dist ~ physeq@sam_data$total_ms_quartile, permutations = 999
)
adtest 

##Total non-glucose MS
#Beta dispersion test
disper.all<-betadisper(nmdsdisp, physeq@sam_data$nonglc_quartile)
TukeyHSD(disper.all) 
plot(disper.all)

#Permanova
adtest <- adonis(bray_dist ~ physeq@sam_data$nonglc_quartile, permutations = 999
)
adtest 

(total_ms_plot <- plot_ordination(physeq, nmds, "samples", color="total_ms_quartile"))
p1 <- total_ms_plot +
  stat_ellipse(type = "t") +
  theme_classic()+
  annotate(
    "text",
    x = 0,
    y = 0.3,
    label = 'betadisper, p > 0.83; adonis, p = 0.766',
    size = 3
  ) +
  theme(legend.position = 'bottom', axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))
p1
(nonglc_plot <- plot_ordination(physeq, nmds, "samples", color="nonglc_quartile"))
p2 <- nonglc_plot +
  stat_ellipse(type = "t") +
  theme_classic() +
  annotate(
    "text",
    x = 0,
    y = 0.3,
    label = 'betadisper, p > 0.29; adonis, p = 0.304',
    size = 3
  ) +
  theme(legend.position = 'bottom',   axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))

p2
g1 <- cowplot::plot_grid(p1, p2, labels = c("C","D"), label_size = 8)
g1

#combine all panels for supp figure 3
g2 <- cowplot::plot_grid(supp_fig_3, g1, nrow = 2, align = 'hv')
g2

ggsave("output/ordination/supp_fig_3_072122.tiff",
       plot = g2,
       width = 7,
       height = 5.7,
       dpi = 1000)
