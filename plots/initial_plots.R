#################
#03/09/2026
#Initial plots to see if there's any trends
################

#load libraries
library(tidyverse)

lipid <- read.csv("./data/data/R_Master - Sheet1.csv")

#make grouping variables factors
lipid$Condition <- as.factor(lipid$Condition)
lipid$Site <- as.factor(lipid$Site)
lipid$Species <- as.factor(lipid$Species)

#check it out
str(lipid)
summary(lipid)

#lets look at just lipid data first (waxes)

lipid_clean <- lipid %>%
  drop_na(avg_lipid) #lost 5 ; same ones missing avg-Lipid are missing std dev obvs

#plot 
ggplot(lipid_clean, aes(x = Condition, y = avg_lipid)) +
  geom_boxplot() +
  geom_jitter(width = 0.1)

#break it down
ggplot(lipid_clean, aes(x = Condition, y = avg_lipids_per_AFDW, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  facet_grid(Species ~ Site) +
  theme_classic() +
  labs(
    x = "Bleaching Condition",
    y = "Total Lipids per AFDW",
    title = "Lipid Content by Bleaching Condition"
  )

#Tissue biomass vs lipid content
ggplot(lipid_clean, aes(x = AFDW_mg_cm2, y = avg_lipid, color = Condition)) +
  geom_point(size = 3) +
  facet_wrap(~Species) +
  theme_classic() +
  labs(
    x = "AFDW (mg/cm²)",
    y = "Average Lipid Content",
    title = "Relationship Between Tissue Biomass and Lipid Content"
  )

########################################
#now let's look at the rest of the lipids
########################################

#make long 
lipid_long <- lipid_clean %>%
  pivot_longer(
    cols = c(avg_tag, avg_st, avg_ampl, avg_pspi, avg_pc, avg_lpc, avg_phospholipids),
    names_to = "lipid_class",
    values_to = "value"
  )

#plot
ggplot(lipid_long, aes(x = Condition, y = value, fill = lipid_class)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Species) +
  theme_classic() +
  labs(
    y = "Proportion of Lipid Classes",
    title = "Lipid Class Composition by Bleaching Condition"
  )

#another check for species differences
ggplot(lipid_clean, aes(x = Species, y = avg_lipids_per_AFDW, fill = Condition)) +
  geom_boxplot() +
  theme_classic()


############################################################################
#quick stats just to see if there's a signal even though plots don't look like it 
############################################################################

model1 <- lm(avg_lipids_per_AFDW ~ Condition + Species + Site, data = lipid_clean)

anova(model1)
# Analysis of Variance Table
# 
# Response: avg_lipids_per_AFDW
# Df   Sum Sq Mean Sq F value  Pr(>F)  
# Condition  1    42858   42858  0.0320 0.85932  
# Species    1  8376078 8376078  6.2624 0.01895 *
#   Site       1  4130937 4130937  3.0885 0.09062 .
# Residuals 26 34775217 1337508                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(model1)
# Call:
#   lm(formula = avg_lipids_per_AFDW ~ Condition + Species + Site, 
#      data = lipid_clean)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1990.62  -576.42   -43.88   555.32  2480.31 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  6383.90     426.94  14.953 2.78e-14 ***
#   ConditionUB    10.72     429.73   0.025   0.9803    
# SpeciesPAST -1206.18     439.06  -2.747   0.0108 *  
#   SiteGMK      -761.04     433.04  -1.757   0.0906 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1157 on 26 degrees of freedom
# Multiple R-squared:  0.2652,	Adjusted R-squared:  0.1804 
# F-statistic: 3.128 on 3 and 26 DF,  p-value: 0.04283
###################################
#Interpretation
#Bleached and unbleached colonies have essentially the same lipid levels in this dataset (p = 0.859)
#Lipid content differs between species (p = 0.01)
#Site is not significant at 0.05 but suggestive (0.09) --> could be small sample size or slight envs effect
#bleached OFRA colonies at FTC have an average lipid level ≈ 6384
#Unbleached colonies only differ by ~11 units (+10.72)
#PAST colonies have about 1200 fewer lipids per AFDW than OFRA -->  a large biological difference?
#Colonies at GMK have about 760 fewer lipids than FTC --> kind of large but moderately sig (p=0.09)
#F-statistic p = 0.0428, At least one predictor significantly explains lipid variation --> gotta be species
###################################


#now go into species more and see if there's an effect
model2 <- lm(avg_lipids_per_AFDW ~ Condition * Species + Site, data = lipid_clean)
anova(model2)
# Analysis of Variance Table
# 
# Response: avg_lipids_per_AFDW
# Df   Sum Sq Mean Sq F value  Pr(>F)  
# Condition          1    42858   42858  0.0312 0.86118  
# Species            1  8376078 8376078  6.1012 0.02068 *
#   Site               1  4130937 4130937  3.0090 0.09512 .
# Condition:Species  1   453800  453800  0.3306 0.57047  
# Residuals         25 34321418 1372857                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(model2)
# Call:
#   lm(formula = avg_lipids_per_AFDW ~ Condition * Species + Site, 
#      data = lipid_clean)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2065.13  -694.80   -67.86   678.06  2391.12 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               6262.6      481.2  13.013 1.24e-12 ***
#   ConditionUB                221.2      568.8   0.389   0.7007    
# SpeciesPAST               -931.2      653.2  -1.426   0.1664    
# SiteGMK                   -775.7      439.5  -1.765   0.0898 .  
# ConditionUB:SpeciesPAST   -508.2      883.8  -0.575   0.5705    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1172 on 25 degrees of freedom
# Multiple R-squared:  0.2748,	Adjusted R-squared:  0.1587 
# F-statistic: 2.368 on 4 and 25 DF,  p-value: 0.07994


##################################################
#Interpretation
#Bleaching condition does not explain lipid differences (p = 0.861)
#The two coral species have different lipid levels (p = 0.02068)
#may be site-level environmental differences (p = 0.09512)
#Bleaching does not affect the species differently.
#Unbleached colonies have ~221 more lipids, but that difference is tiny relative to the variation.
#adding the interaction increases uncertainty and reduces statistical power -> but the direction of the effect stayed the same
#Colonies at GMK appear to have lower lipid reserves than FTC.
##################################################

ggplot(lipid_clean, aes(Species, avg_lipids_per_AFDW, fill = Species)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  theme_classic()
#OFRA slightly higher than PAST


#remove bleaching to see if that makes other drivers stand out : So it estimates lipid levels based on species and site only.
lm(avg_lipids_per_AFDW ~ Species + Site, data = lipid_clean)
# Call:
#   lm(formula = avg_lipids_per_AFDW ~ Species + Site, data = lipid_clean)
# 
# Coefficients:
#   (Intercept)  SpeciesPAST      SiteGMK  
# 6390.1      -1207.2       -760.3  

##################################################
#Interpretation
#OFRA colonies at FTC have ~ 6390 lipids per AFDW
#PAST colonies have ~1207 fewer lipids per AFDW than OFRA.
#Colonies at GMK have ~760 fewer lipids per AFDW than FTC.
##################################################
ggplot(lipid_clean, aes(Species, avg_lipids_per_AFDW, fill = Site)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  theme_classic() +
  labs(
    y = "Lipids per AFDW",
    title = "Lipid reserves by species and site"
  )



############################
#Lipid reserves are more species-dependent than bleaching-dependent, suggesting intrinsic physiological differences between coral species.
#condition has no detectable effect on lipid reserves
#OFRA stores substantially more lipids than PAST.
#GMK colonies may have slightly lower lipid reserves than FTC