## Title: 03_food_group_mono_intake_proportion
## Author: Jules Larke
## Date: 053122
## Purpose: visualize top food/ingredient contributors by proportion of monosaccharide intake (figure 4)

library(tidyverse)
library(RColorBrewer)
library(cowplot)

data <- read_csv("data/all_items_unadjusted_040722.csv") # using unadjusted values for visualizing data
colnames(data)[1] = 'subject_id'

######### Get sum total of each monosaccharide for unique foods/ingredients. Simple name is the glycopedia food name
data2 <- data %>%
  group_by(`Simple name`) %>%
  summarise(
    Glucose_g = sum(Glucose),
    Galactose_g = sum(Galactose),
    Fructose_g = sum(Fructose),
    Xylose_g = sum(Xylose),
    Arabinose_g = sum(Arabinose),
    Fucose_g = sum(Fucose),
    Rhamnose_g = sum(Rhamnose),
    GalA_g = sum(GalA),
    Mannose_g = sum(Mannose),
    Ribose_g = sum(Ribose)
  )

#Calculate total monosacchaide intake
total_group_mono <- as.data.frame(t(colSums(data2[2:11])))

#Divide each food by total monosaccharide to get proportion
data3 <- data2 %>% group_by(`Simple name`) %>% summarise(
  Glucose_prop = Glucose_g / total_group_mono$Glucose_g,
  Galactose_prop = Galactose_g / total_group_mono$Galactose_g,
  Fructose_prop = Fructose_g / total_group_mono$Fructose_g,
  Xylose_prop = Xylose_g / total_group_mono$Xylose_g,
  Arabinose_prop = Arabinose_g / total_group_mono$Arabinose_g,
  Fucose_prop = Fucose_g / total_group_mono$Fucose_g,
  Rhamnose_prop = Rhamnose_g / total_group_mono$Rhamnose_g,
  GalA_prop = GalA_g / total_group_mono$GalA_g,
  Mannose_prop = Mannose_g / total_group_mono$Mannose_g,
  Ribose_prop = Ribose_g /total_group_mono$Ribose_g
)

data3[2:11] <- data3[2:11]*100 #convert to percent

#Check math
colSums(data3[2:11]) #correct

# fix glycopedia food descriptions
data3$`Simple name` <- str_replace(data3$`Simple name`, " \\s*\\([^\\)]+\\)", "")
data3$`Simple name` <- str_replace(data3$`Simple name`, 'Enriched long white rice', 'Long grain white rice')
data3$`Simple name` <- str_replace(data3$`Simple name`, 'Whole Golden Del apple w/o seed', 'Golden Delicious Apple')
data3$`Simple name` <- str_replace(data3$`Simple name`, 'Yellow banana flesh only', 'Yellow banana')
data3$`Simple name` <- str_replace(data3$`Simple name`, 'Hone nut oat cereal', 'Honey nut oat cereal')

# create dataframe of top 3 food/ingredients for each monosaccharide 
glc <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(Glucose_prop))
glc <- glc[1:3,]

gal <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(Galactose_prop))
gal <- gal[1:3,]

gala <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(GalA_prop))
gala <- gala[1:3,]

fru <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(Fructose_prop))
fru <- fru[1:3,]

ara <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(Arabinose_prop))
ara<- ara[1:3,]

xyl <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(Xylose_prop))
xyl <- xyl[1:3,]

fuc <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(Fucose_prop))
fuc <- fuc[1:3,]

rha <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(Rhamnose_prop))
rha <- rha[1:3,]

gala <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(GalA_prop))
gala <- gala[1:3,]

man <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(Mannose_prop))
man <- man[1:3,]

rib <- data3 %>%
  group_by(`Simple name`) %>%
  arrange(desc(Ribose_prop))
rib <- rib[1:3,]

# create plots for each monosaccharide
p4 <-
  ggplot(data = glc, aes(x = fct_reorder(`Simple name`, Glucose_prop, .desc = TRUE), y = Glucose_prop, fill = `Simple name`)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_manual(values = rep('darkgray', 3)) +
  labs(title = "", x = "", y = "% of glucose intake") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), position = "right") +
  theme(title = element_text(size = 8),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.y.right = element_text(angle = 90),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        legend.title = element_blank(),
        legend.position = 'none'
  )
p4

p5 <-
  ggplot(data = gal, aes(x = fct_reorder(`Simple name`, Galactose_prop, .desc = TRUE), y = Galactose_prop, fill = `Simple name`)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_manual(values = rep('darkgray', 3)) +
  labs(title = "", x = "", y = "% of galactose intake") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), position = "right") +
  theme(title = element_text(size = 8),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.y.right = element_text(angle = 90),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        legend.title = element_blank(),
        legend.position = 'none'
  )
p5

p6 <-
  ggplot(data = fru, aes(x = fct_reorder(`Simple name`, Fructose_prop, .desc = TRUE), y = Fructose_prop, fill = `Simple name`)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_manual(values = rep('darkgray', 3)) +
  labs(title = "", x = "", y = "% of fructose intake") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), position = "right", breaks = c(0,2,4,6,8)) +
  theme(title = element_text(size = 8),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.y.right = element_text(angle = 90),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        legend.title = element_blank(),
        legend.position = 'none'
  )
p6

p7<-
  ggplot(data = man, aes(x = fct_reorder(`Simple name`, Mannose_prop, .desc = TRUE), y = Mannose_prop, fill = `Simple name`)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_manual(values = rep('darkgray', 3)) +
  labs(title = "", x = "", y = "% of mannose intake") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), position = "right") +
  theme(title = element_text(size = 8),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.y.right = element_text(angle = 90),
        # axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        legend.title = element_blank(),
        legend.position = 'none'
  )
p7

p8 <-
  ggplot(data = ara, aes(x = fct_reorder(`Simple name`, Arabinose_prop, .desc = TRUE), y = Arabinose_prop, fill = `Simple name`)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_manual(values = rep('darkgray', 3)) +
  labs(title = "", x = "", y = "% of arabinose intake") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), position = "right") +
  theme(title = element_text(size = 8),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.y.right = element_text(angle = 90),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        legend.title = element_blank(),
        legend.position = 'none'
  )
p8

p9 <-
  ggplot(data = xyl, aes(x = fct_reorder(`Simple name`, Xylose_prop, .desc = TRUE), y = Xylose_prop, fill = `Simple name`)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_manual(values = rep('darkgray', 3)) +
  labs(title = "", x = "", y = "% of xylose intake") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), position = "right") +
  theme(title = element_text(size = 8),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.y.right = element_text(angle = 90),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        legend.title = element_blank(),
        legend.position = 'none'
  )
p9

p10 <-
  ggplot(data = gala, aes(x = fct_reorder(`Simple name`, GalA_prop, .desc = TRUE), y = GalA_prop, fill = `Simple name`)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_manual(values = rep('darkgray', 3)) +
  labs(title = "", x = "", y = "% of GalA intake") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), position = "right", breaks = c(0,4,8,12)) +
  theme(title = element_text(size = 8),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.y.right = element_text(angle = 90),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        legend.title = element_blank(),
        legend.position = 'none'
  )
p10

p11 <-
  ggplot(data = rha, aes(x = fct_reorder(`Simple name`, Rhamnose_prop, .desc = TRUE), y = Rhamnose_prop, fill = `Simple name`)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_manual(values = rep('darkgray', 3)) +
  labs(title = "", x = "", y = "% of rhamnose intake") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), position = "right") +
  theme(title = element_text(size = 8),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.y.right = element_text(angle = 90),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        legend.title = element_blank(),
        legend.position = 'none'
  )
p11

p12 <-
  ggplot(data = rib, aes(x = fct_reorder(`Simple name`, Ribose_prop, .desc = TRUE), y = Ribose_prop, fill = `Simple name`)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_manual(values = rep('darkgray', 3)) +
  labs(title = "", x = "", y = "% of ribose intake") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), position = "right") +
  theme(title = element_text(size = 8),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.y.right = element_text(angle = 90),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        legend.title = element_blank(),
        legend.position = 'none'
  )
p12

p13 <-
  ggplot(data = fuc, aes(x = fct_reorder(`Simple name`, Fucose_prop, .desc = TRUE), y = Fucose_prop, fill = `Simple name`)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_manual(values = rep('darkgray', 3)) +
  labs(title = "", x = "", y = "% of fucose intake") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), position = "right") +
  theme(title = element_text(size = 8),
        # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.y.right = element_text(angle = 90),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.9, vjust = 0.3),
        legend.title = element_blank(),
        legend.position = 'none'
  )
p13


# combine plots
g1 = plot_grid(p10, p12, p13,p4, p5, p6, p7, p8, p9, p11, nrow = 2, align = 'hv')
g1
ggsave(
  file = file.path("output/mono_intake/food_group/figure_4_072122.tiff"),
  plot = g1,
  dpi = 1000,
  height = 5.7,
  width = 6,
  units = "in")
