## Title: 03_food_group_mono_intake
## Author: Jules Larke
## Date: 040822
## Purpose: compare monosaccharide intake across food groups

library(tidyverse)
library(RColorBrewer)
library(cowplot)

##Set file paths
getwd()
wd = list()
wd$output = "/Users/jules.larke/work/project/glycan_library/mapping_foods/monosaccharide_FL100/output/mono_intake/food_group/"


data <- read_csv("data/glycan_metadata_041122.csv") # using adjusted values for visualizing data for cohort
colnames(data)[1] = 'subject_id'

data1 <-
  pivot_longer(
    data,
    c(Glucose,
      Galactose,
      Fructose,
      Xylose,
      Arabinose,
      Fucose,
      Rhamnose,
      GalA,
      Mannose,
      Ribose
    ),
    names_to = "Glycan",
    values_to = "proportion"
  )

data1$glycan_food_class <- factor(
  data1$glycan_food_class,
  levels = rev(c(
    "Vegetables",
    "Fruits",
    "Grain Products",
    "Milk and Milk Products" ,
    "Beans, Peas, Legumes, Nuts, Seeds",
    "Meat, Poultry, Fish, and Mixtures",
    "Sugars, Sweets, and Beverages",
    "Fats, Oils, and Salad Dressings",
    "Eggs"
  )),
  ordered = TRUE
)

data1$Glycan <-
  factor(data1$Glycan, levels = rev(
    c("Glucose",
      "Fructose",
      "Galactose",
      "Arabinose",
      "Xylose",
      'GalA',
      "Mannose",
      "Rhamnose",
      "Ribose",
      "Fucose"
    )
  ), ordered = TRUE)

data1 = data1 %>%
  mutate(proportion_kg = (proportion/1000))

glc = data1 %>%
  filter(Glycan == 'Glucose')

glc = glc %>%
  group_by(glycan_food_class) %>%
  summarise(sum = sum(proportion_kg))

glc$sum[7]

arrow = arrow(ends = "both", angle = 90, length = unit(.1,"cm"))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector[4] = 'cornflowerblue'
col_vector[6] = "darkolivegreen4"
col_vector[8] = 'firebrick'
col_vector = rev(col_vector[1:14])
p1 <-
  ggplot(data = data1, aes(x = glycan_food_class, y = proportion_kg, fill = Glycan)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  #geom_hline(yintercept = c(5,10,15,20),  linetype = "longdash", color = 'gray') +
  # annotate('text', x = 9, y = 30, label = "50.46 kg of glucose for grain products") +
  # annotate("segment", x = 8.8, xend = 7, y = 30, yend = 50.457, arrow = arrow()) +
  scale_fill_manual(values = rev(col_vector[1:10]), labels = rev(c("Glucose",
                                                        "Fructose",
                                                        "Galactose",
                                                        "Arabinose",
                                                        "Xylose",
                                                        'GalA',
                                                        "Mannose",
                                                        "Rhamnose",
                                                        "Ribose",
                                                        "Fucose"
  ))) +
  labs(title = "", x = "", y = "Average monosaccharides consumed across recalls (kg/day per 1000 kcal)") +
  # scale_y_continuous(
  #   expand = c(0, 0),
  #   labels = scales::percent,
  #   breaks = c(0, 10, 20, 30, 40),
  #   limits = c(0, 1)
  # ) +
  #scale_x_discrete(expand = c(0, 0), breaks = c(0, 10, 20, 30, 40)) +
  theme_classic() +
  theme(title = element_text(size = 8),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)
  )
p1
legend <- get_legend(p1)

p1 <-
  ggplot(data = data1, aes(x = glycan_food_class, y = proportion_kg, fill = Glycan)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  #geom_hline(yintercept = c(5,10,15,20),  linetype = "longdash", color = 'gray') +
  # annotate('text', x = 9, y = 30, label = "50.46 kg of glucose for grain products") +
  # annotate("segment", x = 8.8, xend = 7, y = 30, yend = 50.457, arrow = arrow()) +
  scale_fill_manual(values = rev(col_vector[1:10]), labels = rev(c("Glucose",
                                                                   "Fructose",
                                                                   "Galactose",
                                                                   "Arabinose",
                                                                   "Xylose",
                                                                   'GalA',
                                                                   "Mannose",
                                                                   "Rhamnose",
                                                                   "Ribose",
                                                                   "Fucose"
  ))) +
  labs(title = "", x = "Food Category", y = "Average monosaccharides consumed\nacross recalls (kg/day per 1000 kcal)") +
  # scale_y_continuous(
  #   expand = c(0, 0),
  #   labels = scales::percent,
  #   breaks = c(0, 10, 20, 30, 40),
  #   limits = c(0, 1)
  # ) +
  #scale_x_discrete(expand = c(0, 0), breaks = c(0, 10, 20, 30, 40)) +
  theme_classic() +
  theme(title = element_text(size = 8),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8),
        legend.position = 'none'
  )
p1

# Look at mono composition after removing glucose, gal, fru
no_glucose = data1 %>%
  filter(Glycan != 'Glucose')

# no_glucose = no_glucose %>%
#   filter(Glycan != 'Galactose')
# 
# no_glucose = no_glucose %>%
#   filter(Glycan != 'Fructose')

p2 <-
  ggplot(data = no_glucose, aes(x = glycan_food_class, y = proportion, fill = Glycan)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  #geom_hline(yintercept = c(500,1000,1500,2000),  linetype = "longdash", color = 'gray') +
  scale_fill_manual(values = rev(col_vector[1:10]), labels = rev(c("Glucose",
                                                                   "Fructose",
                                                                   "Galactose",
                                                                   "Arabinose",
                                                                   "Xylose",
                                                                   'GalA',
                                                                   "Mannose",
                                                                   "Rhamnose",
                                                                   "Ribose",
                                                                   "Fucose"
  ))) +
    labs(title = "", x = "Food Category", y = "Average non-glucose monosaccharides\nconsumed across recalls (g/day per 1000 kcal)") +
  # scale_y_continuous(
  #   expand = c(0, 0),
  #   labels = scales::percent_format(accuracy = 1),
  #   breaks = c(0.2, 0.4, 0.6),
  #   limits = c(0, .6)
  # ) +
  theme_classic() +
  theme(title = element_text(size = 8),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none'

  )
p2

fig_3b <- cowplot::plot_grid(p1, p2, ncol = 1, align = 'hv')
fig_3b
ggsave(
  file = file.path(wd$output,"figure_3B_072122.eps"),
  plot = fig_3b,
  dpi = 1000,
  width = 4.7,
  height = 4,
  units = "in")
ggsave(
  file = file.path(wd$output,"figure_3B_legend_072122.tiff"),
  plot = legend,
  dpi = 1000,
  width = 1,
  height = 4,
  units = "in")
