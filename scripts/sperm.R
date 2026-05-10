###############################
#4/11/2026
#Thesis sperm data
###############################

#load libraries
library(tidyverse)

#load data
sperm<- read.csv("./data/data/Thesis_sperm.csv")

#calculate proportion of colonies that had sperm
sperm$sperm_present <- as.numeric(sperm$stotal > 0)

sperm_mod <- glm(sperm_present ~ Species + Site + State, 
    data = sperm, 
    family = binomial)

#Call:  glm(formula = sperm_present ~ Species + Site + State, family = binomial, 
# data = sperm)
# 
# Coefficients:
#   (Intercept)  SpeciesPAST      SiteGMK      StateUB  
# -2.0218      -0.7558      -1.6516       1.7297  
# 
# Degrees of Freedom: 35 Total (i.e. Null);  32 Residual
# Null Deviance:	    32.44 
# Residual Deviance: 27.55 	AIC: 35.55

################################################################################################
#Interpretation
#Intercept: FTC+OFRA+B -> convert these log odds to probability -> odds = exp(-2.02) ≈ 0.13 -> ~13% of these have sperm 
#SpeciesPAST -> -0.76 -> decreases odds (exp(-0.76) ≈ 0.47) -> PAST colonies are about 50% less likely than OFRA to have sperm
#SiteGMK -> strong negative effect (exp(-1.65) ≈ 0.19) -> GMK colonies are about 80% less likely to have sperm than FTC
#StateUB -> strong positive effect (exp(1.73) ≈ 5.65) -> UB colonies are ~5–6× more likely to have sperm than bleached ones
#really small sample size so it's probably underpowered but we can at least see directions
################################################################################################

#get p-values
summary(sperm_mod)
#none are significant 

#prediction
sperm$pred <- predict(sperm_mod, type = "response")

#build prediction grid to plot
newdat <- expand.grid(
  Species = unique(sperm$Species),
  Site = unique(sperm$Site),
  State = unique(sperm$State)
)

newdat$pred <- predict(sperm_mod, newdat = newdat, type = "response")

#plot
predicted_sperm <- ggplot(newdat, aes(x = State, y = pred, fill = State)) +
  geom_col(position = position_dodge(), color = "black") +
  facet_grid(Species ~ Site) +
  scale_fill_manual(values = c(
    "UB" = "lightpink",
    "B" = "magenta"
  )) +
  labs(
    y = "Predicted probability of sperm presence",
    x = "Bleaching state",
    fill = "Condition"
  ) +
  theme_classic()

ggsave("../Thesis_2026/plots/predicted_sperm.png", predicted_sperm, width = 10, height = 12, dpi = 300)
