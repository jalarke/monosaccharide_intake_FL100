library(tidyverse)
library(RColorBrewer)
library(cowplot)

# Load data
glycan_data <- read.csv("data/glycan_metadata_041122.csv")
colnames(carb_cov)[1] = 'subject_id'
glycan_data %>% group_by(subject_id, RecallNo)

glycan_data <- glycan_data %>% select(!c(GlcA, GalNAc, GlcNAc, Allose))


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
    Ribose = sum(Ribose)
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
    Ribose = mean(Ribose))

glycan_proportion = glycan_total %>%
  mutate(total_ms = select(., Glucose:Ribose) %>% rowSums(na.rm = TRUE)) %>%
  mutate(Glucose = Glucose/total_ms) %>%
  mutate(Galactose = Galactose/total_ms) %>%
  mutate(Fructose = Fructose/total_ms) %>%
  mutate(Xylose = Xylose/total_ms) %>%
  mutate(Arabinose = Arabinose/total_ms) %>%
  mutate(Fucose= Fucose/total_ms) %>%
  mutate(Rhamnose = Rhamnose/total_ms) %>%
  mutate(GalA = GalA/total_ms) %>%
  mutate(Mannose = Mannose/total_ms) %>%
  mutate(Ribose = Ribose/total_ms)

glycan_proportion <- glycan_proportion * 100

avg_mono_proportion <- glycan_proportion %>% summarise(across(everything(), mean))
sd_mono_proportion <- glycan_proportion %>% summarise(across(everything(), sd))
mean_sd <- rbind(avg_mono_proportion, sd_mono_proportion)
mean_sd$subject_id <- NULL
mean_sd$total_ms <- NULL
#write_csv(mean_sd, 'ms_proportion_mean_sd_020622.csv')

mean <- mean_sd[1,]
pie <-  pivot_longer(mean,
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
                     values_to = "grams"
)

pie <- arrange(pie, desc(grams))

pie$Glycan <-
  factor(pie$Glycan, levels = 
    c("Glucose",
      "Fructose",
      "Galactose",
      "Arabinose",
      "Xylose",
      'GalA',
      "Mannose",
      "Ribose",
      "Fucose",
      "Rhamnose"

  ), ordered = TRUE)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector[4] = 'cornflowerblue'
col_vector[6] = "darkolivegreen4"
col_vector[8] = 'firebrick'
col_vector = col_vector[1:14]
#col_vector = rev(col_vector)
library(ggrepel)
# Basic piechart
pie_chart <- ggplot(pie, aes(x = "", y = grams, fill = Glycan)) +
  geom_bar(stat = "identity",
           width = 1,
           color = "white") +
  scale_fill_manual(values = rev(col_vector), labels = 
                      rev(c("Glucose",
                        "Fructose",
                        "Galactose",
                        "Arabinose",
                        "Xylose",
                        'GalA',
                        "Mannose",
                        "Ribose",
                        "Fucose",
                        "Rhamnose"
                        ))) +
  geom_label(data = pie[1:7,],
                    aes(label = Glycan), hjust = 'top', vjust = 'top') +
  coord_polar("y") +
  theme_void() +
  theme(legend.position = "none")
pie_chart
round(mean_sd$GalA[2], 1)
