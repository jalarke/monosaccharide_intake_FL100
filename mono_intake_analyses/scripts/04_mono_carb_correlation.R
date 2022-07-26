## Title: 04_mono_carb_correlation
## Author: Jules Larke
## Date: 041122
## Purpose: Correlations for MS and CHO

library(tidyverse)

##Set file paths
getwd()
wd = list()
wd$data = "/Users/jules.larke/work/project/glycan_library/mapping_foods/monosaccharide_FL100/data/"
wd$output = "/Users/jules.larke/work/project/glycan_library/mapping_foods/monosaccharide_FL100/output/correlation/"

data <- read_csv('data/all_items_unadjusted_040722.csv') # Not energy adjusted
colnames(data)[1] = 'subject_id'

data <- data %>% select(-c(GlcA, GalNAc, GlcNAc, Allose))

glycan_total <- data %>% group_by(subject_id, RecallNo) %>%
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
    total_carb = sum(carb_consumed_g)
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
    total_carb = mean(total_carb)) %>% mutate(total_ms = select(., Glucose:Ribose) %>% rowSums(na.rm = TRUE))


r2 <- cor.test(glycan_total$total_ms, glycan_total$total_carb, method = "pearson")
r2$p.value
round(r2$estimate, 3)
round(r2$p.value, 6)

r_label = paste0("italic(r) ==", round(r2$estimate, 3))
p_label = paste0("italic(p) ==", 6.8e-95)

n3 <- ggplot(glycan_total, aes(x = total_ms, y = total_carb)) +
  geom_point(color = "black",
             size = 2,
             shape = 16) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = 'Average monosaccharides consumed (g/day per 1000 kcal)\n(Glycopedia, analytically derived)') +
  labs(y = "Average carbohydrate consumed (g/day per 1000 kcal)\n(FNDDS, Standard Reference value)") +
  annotate(
    "text",
    x = 30,
    y = 475,
    label = r_label,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = 30,
    y = 450,
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
n3

ggsave(file.path(wd$output,"figure_2_072122.tiff"),
       plot = n3,
       width = 5.7,
       height = 4,
       dpi = 1000)
