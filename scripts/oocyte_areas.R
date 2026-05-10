##########################
#05/1/2026
#Thesis Oocyte Area
##########################

#load libraries 
library(tidyverse)
library(lme4)
library(lmerTest)   # gives p-values for lmer
library(readr)
library(patchwork)

#load data

oocytes <- read_csv("../Thesis_2026/data/data/oocyte_areas.csv")

#let's filter out the PAST for now since there are only 2 
oocytes_main <- oocytes %>%
  filter(Species == "OFRA")

#check distribution
ggplot(oocytes_main, aes(x = Oocyte_Area)) +
  geom_histogram(bins = 30)

#log transform
oocytes_main <- oocytes_main %>%
  mutate(logArea = log(Oocyte_Area))

ggplot(oocytes_main, aes(x = logArea)) +
  geom_histogram(bins = 30)
#still a little skewed 

##trying a mixed effects model first
model <- lmer(logArea ~ Status * Site + (1 | ID), data = oocytes_main)

summary(model)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: logArea ~ Status * Site + (1 | ID)
#    Data: oocytes_main
# 
# REML criterion at convergence: -209
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.6043 -0.3803  0.0912  0.5246  2.4088 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  ID       (Intercept) 0.008621 0.09285 
#  Residual             0.007997 0.08942 
# Number of obs: 123, groups:  ID, 13
# 
# Fixed effects:
#                   Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)       4.928135   0.049695  5.256688  99.167 8.33e-10 ***
# StatusUB         -0.007464   0.077280  5.245978  -0.097    0.927    
# SiteGMK           0.054891   0.091803  6.829324   0.598    0.569    
# StatusUB:SiteGMK -0.083043   0.123092  6.591829  -0.675    0.523    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) SttsUB SitGMK
# StatusUB    -0.643              
# SiteGMK     -0.541  0.348       
# SttsUB:SGMK  0.404 -0.628 -0.746

anova(model)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq   Mean Sq NumDF  DenDF F value Pr(>F)
# Status      0.0050659 0.0050659     1 6.5918  0.6335 0.4538
# Site        0.0003773 0.0003773     1 6.5918  0.0472 0.8346
# Status:Site 0.0036396 0.0036396     1 6.5918  0.4551 0.5229

########################################################################
#Interpretation
#overall: testing whether oocyte size differs by bleaching, site, and whether bleaching effects vary by site.
#Random effects: Colony-to-colony differences are about the same magnitude as within-colony (oocyte-level) variation
#Fixed Effects: Unbleached colonies have slightly smaller oocytes than bleached ones (-0.07) BUT the effect is basically zero and completely non-significant (p = 0.927) .....There was no significant effect of bleaching status on oocyte area.
#Site (GMK): GMK colonies trend slightly larger than FTC (+0.055) BUT also not significant (p= 0.569)....Oocyte area did not differ significantly between sites.
#Interaction: -0.083, p= 0.523 ....effect of bleaching does not differ between sites
########################################################################


#now view it 

theme_set(theme_minimal(base_size = 14))

#first oocyte area by condition and size 
ggplot(oocytes_main, aes(x = Status, y = logArea, fill = Status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_wrap(~Site) +
  labs(
    x = "Bleaching Status",
    y = "Log Oocyte Area",
    title = "Oocyte Area by Bleaching Status Across Sites"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold") 
  ) + scale_fill_manual(values = c("UB" = "lightpink", "B" = "magenta"))

#colony level means 
colony_means <- oocytes_main %>%
  group_by(ID, Site, Status) %>%
  summarize(mean_logArea = mean(logArea), .groups = "drop")

ggplot(colony_means, aes(x = Status, y = mean_logArea, fill = Status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  facet_wrap(~Site) +
  labs(
    x = "Bleaching Status",
    y = "Mean Log Oocyte Area (per colony)",
    title = "Colony-Level Oocyte Investment"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) + scale_fill_manual(values = c("UB" = "lightpink", "B" = "magenta"))


#raw scale
ggplot(oocytes_main, aes(x = Status, y = OocyteArea, fill = Status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_wrap(~Site) +
  scale_y_continuous(trans = "log10") +
  labs(
    x = "Bleaching Status",
    y = "Oocyte Area (log scale)",
    title = "Oocyte Area (Raw Values, Log Scaled Axis)"
  ) +
  theme(
    legend.position = "none"
  )


#make multipanel fig 
colors <- c("UB" = "lightpink", "B" = "magenta")

#panel A
p1 <- ggplot(oocytes_main, aes(x = Status, y = logArea, fill = Status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  facet_wrap(~Site) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Bleaching Status",
    y = "Log Oocyte Area",
    title = "A. Oocyte-Level Variation"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

#panel B
colony_means <- oocytes_main %>%
  group_by(ID, Site, Status) %>%
  summarize(mean_logArea = mean(logArea), .groups = "drop")

p2 <- ggplot(colony_means, aes(x = Status, y = mean_logArea, fill = Status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  facet_wrap(~Site) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Bleaching Status",
    y = "Mean Log Oocyte Area",
    title = "B. Colony-Level Means"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

#combine
combined_plot <- p1 / p2  # stacked vertically

combined_plot

#try side by side 
combined_plot2 <- p1 + p2

combined_plot2

#save 
ggsave("../Thesis_2026/plots/oocyte_area_combined2.png", combined_plot2, width = 10, height = 12, dpi = 300)
