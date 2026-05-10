##########################
#04/10/2026
#Thesis Fecundity
##########################

#load libraries 
library(tidyverse)
library(lme4)       # for mixed models
library(lmerTest)   # p-values for lmer
library(ggpubr)     # nicer plots
library(patchwork)

#########################################
# let's look at fecundity 
#########################################

#load data
fec<- read.csv("./data/data/Thesis_fec.csv")

fec <- fec %>%
  mutate(
    Species = as.factor(Species),
    Site = as.factor(Site),
    Status = as.factor(Status)
  )


#break it down into 3 sets since we had a lot of nonrepros (zeros) and ones that didnt have 2 fecund polyps

fec <- fec %>%
  mutate(
    # Binary reproduction (colony level) #where they repro yes or no
    repro = ifelse(F1 > 0 | F2 > 0, 1, 0), 
    
    # Mean fecundity including zeros #averaging all even with zeros 
    fecundity_all = (F1 + F2) / 2,
    
    # Conditional fecundity (only nonzero polyps) 
    fecundity_nonzero = case_when(
      F1 > 0 & F2 > 0 ~ (F1 + F2) / 2,
      F1 > 0 & F2 == 0 ~ F1,
      F1 == 0 & F2 > 0 ~ F2,
      TRUE ~ NA_real_
    )
  )

summary(fec)
table(fec$repro)
# 0  1 
# 25 10 #25 nonrepro, 10 repro 


#zfor the first model we'll look at the probability of a sample being reproductive

mod_repro <- glm(repro ~ Species * Site * Status,
                 data = fec,
                 family = "binomial")

summary(mod_repro) #are bleached corals less likely to reproduce?

# Call:
#   glm(formula = repro ~ Species * Site * Status, family = "binomial", 
#       data = fec)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                     1.386      1.118   1.240    0.215
# SpeciesPAST                   -19.952   2917.013  -0.007    0.995
# SiteGMK                        -2.485      1.607  -1.546    0.122
# StatusUB                       -2.079      1.414  -1.470    0.141
# SpeciesPAST:SiteGMK             2.485   7145.193   0.000    1.000
# SpeciesPAST:StatusUB           19.036   2917.013   0.007    0.995
# SiteGMK:StatusUB                2.485      2.021   1.230    0.219
# SpeciesPAST:SiteGMK:StatusUB  -19.442   8504.481  -0.002    0.998
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 41.879  on 34  degrees of freedom
# Residual deviance: 30.186  on 27  degrees of freedom
# AIC: 46.186
# 
# Number of Fisher Scoring iterations: 17

##########################################################
#Interpretation
#“What affects whether a coral reproduces at all?”
#model is struggling because PAST groups basically never reproduce -> why there's huge number (-19.95)
#the intercept is the log-odds of reproduction for the reference groups which are: Species = OFRA, Site = FTC, Status = B -> let's convert this to probability
exp(1.386) / (1 + exp(1.386)) #~80% -> so OFRA, FTC, Bleached corals have ~80% chance of reproducing
#don't really have enough data for this much less the interactions (big std errs and high p values)
##########################################################

# we need to simplify this model
glm(repro ~ Species + Site + Status,
    family = "binomial",
    data = fec)

# Call:  glm(formula = repro ~ Species + Site + Status, family = "binomial", 
#            data = fec)
# 
# Coefficients:
#   (Intercept)  SpeciesPAST      SiteGMK     StatusUB  
# 0.4592      -2.6360      -1.0894      -0.4457  
# 
# Degrees of Freedom: 34 Total (i.e. Null);  31 Residual
# Null Deviance:	    41.88 
# Residual Deviance: 33.99 	AIC: 41.99


#let's also interpret species seperately
fec_OFRA <- fec %>% filter(Species == "OFRA")

glm(repro ~ Site + Status,
    family = "binomial",
    data = fec_OFRA)


# Call:  glm(formula = repro ~ Site + Status, family = "binomial", data = fec_OFRA)
# 
# Coefficients:
#   (Intercept)      SiteGMK     StatusUB  
# 0.6931      -1.0296      -0.9163  
# 
# Degrees of Freedom: 20 Total (i.e. Null);  18 Residual
# Null Deviance:	    28.68 
# Residual Deviance: 26.4 	AIC: 32.4

##################################################################################
#Interpretation
#Intercept is the baseline group (OFRA+FTC+B)-> convert that to probability exp(0.6931) / (1 + exp(0.6931)) = 0.67 -> 67% of B OFRA colonies at FTC are repro
#SiteGMK is the effect of being at GMK (vs FTC) -> convert to odds exp(-1.0296) ≈ 0.36 -> 1-0.36 -> odds of reproducing at GMK are ~64% lower than FTC -> now convert that to probability -> logit = 0.6931 - 1.0296 = -0.3365  prob ≈ 0.42 -> OFRA at GMK have ~42% probability of reproducing vs 67% at FTC
#StatusUB -> exp(-0.9163) ≈ 0.40 -> odds of reproducing are ~60% lower in UB colonies -> now probability -> logit = 0.6931 - 0.9163 = -0.2232 -> prob ≈ 0.44 -> Unbleached OFRA have ~44% probability vs 67% in bleached
#Null deviance:     28.68  , Residual deviance: 26.4 , The model only explains a small amount of variation
##################################################################################

#check proportions
fec_OFRA %>%
  group_by(Site, Status) %>%
  summarise(
    prop_repro = mean(repro),
    n = n()
  )

# Site  Status prop_repro     n
# <fct> <fct>       <dbl> <int>
#   1 FTC   B           0.8       5
# 2 FTC   UB          0.333     6
# 3 GMK   B           0.25      4
# 4 GMK   UB          0.333     6

#maybe also try a Fisher's Exact Test
table(fec$Status, fec$repro)

fisher.test(table(fec$Status, fec$repro))

# Fisher's Exact Test for Count Data
# 
# data:  table(fec$Status, fec$repro)
# p-value = 0.7118
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.1188241 3.7872230
# sample estimates:
# odds ratio 
#  0.6745758 



#for the second model we'll look to see if among reproductive corals, does bleaching reduce egg count?

fec_repro <- fec %>% filter(repro == 1) #just get repro ones 

mod_fec <- lm(fecundity_nonzero ~ Species, data = fec_repro)

summary(mod_fec)

# Call:
#   lm(formula = fecundity_nonzero ~ Species, data = fec_repro)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -11.056  -4.931   2.222   4.444   9.444 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   12.556      2.538   4.947  0.00112 **
#   SpeciesPAST   -7.556      8.025  -0.941  0.37402   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 7.613 on 8 degrees of freedom
# Multiple R-squared:  0.09974,	Adjusted R-squared:  -0.01279 
# F-statistic: 0.8864 on 1 and 8 DF,  p-value: 0.374


#for this third model we'll do a mixed model at polyp level and include the avg ove 2 polyps 

#make data long first 
fec_long <- fec %>%
  pivot_longer(cols = c(F1, F2),
               names_to = "Polyp",
               values_to = "Eggs")

#run mixed model 
mod_mixed <- lmer(Eggs ~ Species * Site * Status + (1 | ID),
                  data = fec_long)

summary(mod_mixed)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Eggs ~ Species * Site * Status + (1 | ID)
#    Data: fec_long
# 
# REML criterion at convergence: 414.5
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.0431 -0.1827 -0.0425  0.0000  3.8996 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  ID       (Intercept) 33.04    5.748   
#  Residual             18.67    4.321   
# Number of obs: 70, groups:  ID, 35
# 
# Fixed effects:
#                              Estimate Std. Error     df t value Pr(>|t|)   
# (Intercept)                     8.600      2.911 27.000   2.954  0.00643 **
# SpeciesPAST                    -8.600      4.117 27.000  -2.089  0.04629 * 
# SiteGMK                        -7.600      4.367 27.000  -1.740  0.09319 . 
# StatusUB                       -2.100      3.942 27.000  -0.533  0.59858   
# SpeciesPAST:SiteGMK             7.600      8.362 27.000   0.909  0.37147   
# SpeciesPAST:StatusUB            2.933      5.575 27.000   0.526  0.60306   
# SiteGMK:StatusUB                4.683      5.762 27.000   0.813  0.42342   
# SpeciesPAST:SiteGMK:StatusUB   -5.517     10.597 27.000  -0.521  0.60691   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) SpPAST SitGMK SttsUB SpPAST:SGMK SPAST:SU SGMK:S
# SpeciesPAST -0.707                                                 
# SiteGMK     -0.667  0.471                                          
# StatusUB    -0.739  0.522  0.492                                   
# SpPAST:SGMK  0.348 -0.492 -0.522 -0.257                            
# SpcPAST:SUB  0.522 -0.739 -0.348 -0.707  0.364                     
# StGMK:SttUB  0.505 -0.357 -0.758 -0.684  0.396       0.484         
# SPAST:SGMK: -0.275  0.389  0.412  0.372 -0.789      -0.526   -0.544


#I just want to see proportions so I can add that to RF abstract
fec %>%
  group_by(Status) %>%
  summarise(
    prop_repro = mean(repro),
    n = n()
  )

#  Status prop_repro     n
# <fct>       <dbl> <int>
#   1 B           0.333    15
# 2 UB          0.25     20
#######################
#now try and visualize 
#######################

#first let's look at the prop repro
plot_repro <- fec %>%
  group_by(Species, Site, Status) %>%
  summarise(prop_repro = mean(repro),
            n = n()) %>%
  ggplot(aes(x = Status, y = prop_repro, fill = Status)) +
  geom_col() +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.5) +
  facet_grid(Species ~ Site) +
  scale_fill_manual(values = c("B" = "magenta", "UB" = "lightpink")) +
  ylim(0, 1) +
  labs(y = "Proportion Reproductive", x = "Bleaching Status") +
  theme_bw() +
  theme(legend.position = "none")

plot_repro


#now let's plot fec of only repro samples

fec_counts <- fec_repro %>%
  group_by(Species, Site, Status) %>%
  summarise(n = n(), .groups = "drop")

plot_fec <- fec_repro %>%
  ggplot(aes(x = Status, y = fecundity_nonzero, fill = Status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  geom_jitter(aes(color = Status), width = 0.08, size = 2, alpha = 0.7) +
  geom_text(data = fec_counts,
            aes(x = Status, y = Inf, label = paste0("n=", n)),
            vjust = 1.5, size = 3, inherit.aes = FALSE) +
  facet_grid(Species ~ Site) +
  scale_fill_manual(values = c("B" = "magenta", "UB" = "lightpink")) +
  scale_color_manual(values = c("B" = "magenta", "UB" = "lightpink")) +
  coord_cartesian(ylim = c(0, max(fec_repro$fecundity_nonzero) * 1.15)) +
  labs(y = "Eggs per Polyp", x = "Bleaching Status") +
  theme_bw() +
  theme(legend.position = "none")

plot_fec


#now let's look at raw fecundity that includes the samples with zeros
# get sample sizes for ALL colonies
fec_all_counts <- fec %>%
  group_by(Species, Site, Status) %>%
  summarise(n = n(), .groups = "drop")

plot_fec_all <- fec %>%
  ggplot(aes(x = Status, y = fecundity_all, fill = Status)) +
  
  # boxplot
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  
  # jitter points
  geom_jitter(aes(color = Status), width = 0.08, size = 2, alpha = 0.6) +
  
  # sample sizes
  geom_text(data = fec_all_counts,
            aes(x = Status, y = Inf, label = paste0("n=", n)),
            vjust = 1.5, size = 3, inherit.aes = FALSE) +
  
  facet_grid(Species ~ Site) +
  
  # consistent colors across ALL plots
  scale_fill_manual(values = c(
    "B" = "magenta",
    "UB" = "lightpink"
  )) +
  scale_color_manual(values = c(
    "B" = "magenta",
    "UB" = "lightpink"
  )) +
  
  # prevent label cutoff
  coord_cartesian(ylim = c(0, max(fec$fecundity_all) * 1.15)) +
  
  labs(y = "Eggs per Polyp (incl. zeros)", x = "Bleaching Status") +
  theme_bw() +
  theme(legend.position = "none")

plot_fec_all



#need to add sample sizes to plots I think

fec_summary <- fec %>%
  group_by(Species, Site, Status) %>%
  summarise(
    prop_repro = mean(repro),
    n = n()
  )



#combine these plots to see if it looks better 


figure_combined <- plot_repro / plot_fec / plot_fec_all +
  plot_annotation(tag_levels = "A")

figure_combined

figure_combined <- (plot_repro | plot_fec | plot_fec_all) +
  plot_layout(widths = c(1, 1.5, 1.5)) +
  plot_annotation(tag_levels = "A")

ggsave("../Thesis_2026/plots/fecundity_figure.png",
       figure_combined,
       width = 12,
       height = 10,
       dpi = 300)






