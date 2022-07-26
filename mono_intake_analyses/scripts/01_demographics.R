## Title: Participant count across recruitment strata for FL100
## Author: Jules Larke
## Date: 030422
## Purpose: Create table and figure for participant demographic data. Data for table 1, table 2 and supplemental figure 1

library(tidyverse)

# Load data
glycan_data <- read.csv("data/glycan_metadata_041122.csv")
glycan_total <- glycan_data %>% group_by(subject_id, RecallNo) %>%
  summarise(
    Glucose = sum(Glucose),
    Galactose = sum(Galactose),
    Fructose = sum(Fructose),
    Arabinose = sum(Arabinose),
    Xylose = sum(Xylose),
    Fucose = sum(Fucose),
    Rhamnose = sum(Rhamnose),
    GlcA = sum(GlcA),
    GalA = sum(GalA),
    GlcNAc = sum(GlcNAc),
    GalNAc = sum(GalNAc),
    Mannose = sum(Mannose),
    Allose = sum(Allose),
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
    GlcA = mean(GlcA),
    GalA = mean(GalA),
    GlcNAc = mean(GlcNAc),
    GalNAc = mean(GalNAc),
    Mannose = mean(Mannose),
    Allose = mean(Allose),
    Ribose = mean(Ribose))

data1 <-
  pivot_longer(
    glycan_total,
    c(Glucose,
      Galactose,
      Fructose,
      Xylose,
      Arabinose,
      Fucose,
      Rhamnose,
      GlcA,
      GalA,
      GlcNAc,
      GalNAc,
      Mannose,
      Allose,
      Ribose
    ),
    names_to = "Glycan",
    values_to = "grams"
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
      "Fucose",
      "GlcNAc",
      "GlcA",
      "GalNAc",
      "Allose"
    )
  ), ordered = TRUE)
levels(data1$subject_id)

data1$subject_id = factor(data1$subject_id)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector[4] = 'cornflowerblue'
col_vector[6] = "darkolivegreen4"
col_vector[8] = 'firebrick'

pheno <- distinct(glycan_data, subject_id, .keep_all = TRUE)
pheno_var <- pheno %>% select(subject_id, age, sex.factor, bmi)
#Age: Use study bins "18 - 33", "34 - 49", "50 - 65"
pheno_var <- mutate(pheno_var, Age_Category = "34_49")
pheno_var$Age_Category[as.numeric(pheno_var$age) < 34] <- "18_33"
pheno_var$Age_Category[as.numeric(pheno_var$age) > 49] <- "50_65"

#BMI bins
pheno_var <- mutate(pheno_var, BMI_Category = "Overweight")
pheno_var$BMI_Category[as.numeric(pheno_var$bmi) < 25 ] <- "Normal"
pheno_var$BMI_Category[as.numeric(pheno_var$bmi) > 30 ] <- "Obese"

pheno_var$Ethnicity <- pheno$Ethnicity 

pheno_var %>% summarise(mean_bmi = mean(bmi), sd_bmi = sd(bmi))
pheno_var %>% group_by(sex.factor) %>% summarise(mean_age = mean(age), sd_age = sd(age))

strata <- as.data.frame(xtabs(~ Age_Category + BMI_Category + sex.factor, data = pheno_var))

strata$BMI_Category <- factor(strata$BMI_Category, levels = c('Normal', 'Overweight', 'Obese'), ordered = TRUE)

p1 <- ggplot(data = strata) +
  aes(x = Age_Category, y = Freq, fill = sex.factor) +
  geom_bar(aes(fill = sex.factor), stat = 'identity', position = "dodge") +
  labs(x = 'Age range', y = 'Count', title = '') +
  theme_classic(base_size = 12) + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c('firebrick', 'royalblue')) +
  facet_wrap(.~BMI_Category, scales = "free")
p1
ggsave(
  file = file.path("output/participants/participant_strata_n180.tiff"),
  plot = p1,
  dpi = "retina",
  width = 8,
  height = 6,
  units = "in")

pheno_var$Ethnicity <- pheno$Ethnicity 
eth <- as.data.frame(xtabs(~ sex.factor + Ethnicity, data = pheno_var))
eth %>% filter(sex.factor == 'Male') %>% mutate(percent = Freq/sum(Freq)*100)
eth %>% filter(sex.factor == 'Female') %>% mutate(percent = Freq/sum(Freq)*100)
eth %>% group_by(Ethnicity) %>% mutate(sex.factor = NULL, percent = Freq/sum(Freq)*100)

eth <- as.data.frame(xtabs(~ Ethnicity, data = pheno_var))
eth %>% mutate(percent = Freq/sum(Freq)*100)
97/180

#Table 2: linear models to examine mono intake across covariates
# model each monosaccharide by age, BMI and sex. Age and BMI are scaled.

mono_pheno <- merge(glycan_total, pheno_var)

#Glucose
glucose_lm <- lm(Glucose ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
summary(glucose_lm)$coefficients
plot(glucose_lm)
shapiro.test(glucose_lm$residuals)

glc <- summary(glucose_lm)$coefficients

#Fructose
fructose_lm <- lm(Fructose ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(fructose_lm$residuals) # not normally distributed

fructose_lm <- lm(log10(Fructose) ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
summary(fructose_lm)$coefficient 
shapiro.test(fructose_lm$residuals)
plot(fructose_lm, which = 3)

fru <- summary(fructose_lm)$coefficients
fru

#Galactose
galactose_lm <- lm(Galactose ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
plot(galactose_lm)
shapiro.test(galactose_lm$residuals)

galactose_lm <- lm(log10(Galactose) ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(galactose_lm$residuals)
plot(galactose_lm, which = 3)

gal <- summary(galactose_lm)$coefficient


#Arabinose
arabinose_lm <- lm(Arabinose ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(arabinose_lm$residuals)

arabinose_lm <- lm(log10(Arabinose) ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(arabinose_lm$residuals)
plot(arabinose_lm, which = 3)

ara <- summary(arabinose_lm)$coefficient

#Xylose
xylose_lm <- lm(Xylose ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(xylose_lm$residuals)

xylose_lm <- lm(log10(Xylose) ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(xylose_lm$residuals)
plot(xylose_lm, which = 3)

xyl <- summary(xylose_lm)$coefficient

#Fucose
fucose_lm <- lm(Fucose ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(fucose_lm$residuals)

fucose_lm <- lm(log10(Fucose) ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(fucose_lm$residuals)
plot(fucose_lm, which = 3)

fuc<- summary(fucose_lm)$coefficient

#Rhamnose
rhamnose_lm <- lm(Rhamnose ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(rhamnose_lm$residuals)

rhamnose_lm <- lm(log10(Rhamnose) ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(rhamnose_lm$residuals)
plot(rhamnose_lm, which = 3)

rham <- summary(rhamnose_lm)$coefficient

#GalA
gala_lm <- lm(GalA ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(gala_lm$residuals)

gala_lm <- lm(sqrt(GalA) ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(gala_lm$residuals)
plot(gala_lm, which = 3)

gala <- summary(gala_lm)$coefficient

#Mannose
mann_lm <- lm(Mannose ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(mann_lm$residuals)

mann_lm <- lm(log10(Mannose) ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(mann_lm$residuals) 
plot(mann_lm, which = 3)

mann <- summary(mann_lm)$coefficient

#Ribose
rib_lm <- lm(Ribose ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(rib_lm$residuals)

rib_lm <- lm(sqrt(Ribose) ~ scale(age, scale = F) + scale(bmi, scale = F) + sex.factor, data = mono_pheno)
shapiro.test(rib_lm$residuals) 
plot(allo_lm, which = 3)

rib <- summary(rib_lm)$coefficient
rib <- str_c(round(rib[1,]),3)

lm_mono <- cbind(glc, fru, gal, ara, xyl, fuc, rham, gala, mann, rib)
lm_mono$
write.csv(lm_mono, 'lm_mono_covariates.csv')
