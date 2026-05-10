##########################
#04/9/2026
#Symbiont Data Comparisons
##########################

#load libraries 
library(tidyverse)
library (vegan)

#load data 
symb <- read.csv("./data/data/gillian_compiled_qpcr_data.csv")

#######################
#PERMANOVA
#######################
#make proportions matrix
symb_mat <- symb[, c("Prop_A","Prop_B","Prop_C","Prop_D")]

#get rid of 0s 
rowSums(symb_mat) == 0

symb <- symb[rowSums(symb[, c("Prop_A","Prop_B","Prop_C","Prop_D")]) > 0, ]

symb_mat <- symb[, c("Prop_A","Prop_B","Prop_C","Prop_D")]

#PERMANOVA
adonis2(symb_mat ~ Species * Site * Condition, data = symb, method = "bray")

# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = symb_mat ~ Species * Site * Condition, data = symb, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# Model     7   9.9011 0.94719 64.057  0.001 ***
#   Residual 25   0.5520 0.05281                  
# Total    32  10.4531 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


########################################################################
#INTERPRET
#R2 = 0.94719 -> means model explains ~94% of variation (suspiciously strong?-> could be driven by species differences)
#F = 64.06/ p = 0.001 -> symbiont community composition differs across predictors
#this is just model vs residual so can't tell what variable is driving the differences
########################################################################

#try and separate variables
adonis2(symb_mat ~ Species * Site * Condition, data = symb, method = "bray", by = "terms")

#Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = symb_mat ~ Species * Site * Condition, data = symb, method = "bray", by = "terms")
# Df SumOfSqs      R2        F Pr(>F)    
# Species                 1   3.7796 0.36157 171.1681  0.001 ***
#   Site                    1   4.3913 0.42010 198.8725  0.001 ***
#   Condition               1   0.0800 0.00765   3.6223  0.039 *  
#   Species:Site            1   1.0190 0.09749  46.1495  0.001 ***
#   Species:Condition       1   0.0237 0.00227   1.0739  0.338    
# Site:Condition          1   0.1616 0.01546   7.3201  0.003 ** 
#   Species:Site:Condition  1   0.4458 0.04265  20.1912  0.001 ***
#   Residual               25   0.5520 0.05281                    
# Total                  32  10.4531 1.00000                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

########################################################################
#INTERPRET
#Species & site are dominant drivers(R2 = 0.36; R2 = 0.42) -> site matters most
#SpeciesxSite interaction is next biggest (R2 = 0.10) -> The effect of site depends on species (e.g., OFRA responds differently across sites than PAST)
#Condition (R2 = 0.007) so less than 1% explains but it is significant -> suggests subtle shifts, not complete community restructuring? - ask Dusty
#Species × Site × Condition (3-way interaction) (R2 = 4.3, p = 0.001) -> suggest bleaching effects depend on both species AND site -> Some species at some sites may shift symbionts under bleaching? -> could also be that FTC colonies already have clade D and may not shift?
########################################################################

#Symbiont community composition is primarily structured by environmental differences between sites and host species identity, while bleaching condition exerts a smaller but significant influence that varies depending on site and species.


#check dispersion effects since PERMANOVA was used
dist <- vegdist(symb_mat, method = "bray")
bd <- betadisper(dist, symb$Condition)
anova(bd)

# Analysis of Variance Table
# 
# Response: Distances
# Df  Sum Sq  Mean Sq F value Pr(>F)
# Groups     1 0.01508 0.015081  0.5746 0.4542
# Residuals 31 0.81369 0.026248  

#no evidence of differences in dispersion (variance) between groups

#Differences in symbiont community composition were not driven by heterogeneity of dispersion (PERMDISP, p = 0.45), supporting the validity of PERMANOVA results.


#######################
#NMDS plot
#######################

#run NMDS
nmds <- metaMDS(symb_mat, distance = "bray", k = 2, trymax = 100)

#got warning so check stress
nmds$stress #0.000472301 -> we are on firreeeee 
#could also be that species/site had such a strong effect, only 4 variables, and they're mostly made up of clade D 

#NMDS ordination showed an excellent fit to the data (stress = 0.0005), indicating that the two-dimensional representation accurately reflects dissimilarities among samples.


#now plot it 
plot(nmds)

#but better

# Extract NMDS scores
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))

# Combine with metadata
nmds_scores <- cbind(nmds_scores, symb)

# Plot

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Species, shape = Condition), size = 3) +
  facet_wrap(~ Site) +
  theme_classic() +
  labs(title = "NMDS of Symbiont Community Composition",
       subtitle = paste("Stress =", round(nmds$stress, 4)))


#let's add the clades to it 
# Fit vectors
clade_fit <- envfit(nmds, symb_mat, permutations = 999)

# Check which vectors are significant
clade_fit

# Extract coordinates for arrows
vectors <- as.data.frame(scores(clade_fit, display = "vectors"))
vectors$clade <- rownames(vectors)

#also use jitter to seperate points a little 


library(ggplot2)
library(vegan)
library(grid)  # for arrow()

# Make sure Condition is a factor with proper labels
nmds_scores$Condition <- factor(nmds_scores$Condition, 
                                levels = c("UB", "B"),
                                labels = c("Unbleached", "Bleached"))

# Add clade vectors
clade_fit <- envfit(nmds, symb_mat, permutations = 999)
vectors <- as.data.frame(scores(clade_fit, display = "vectors"))
vectors$clade <- rownames(vectors)

# Make sure Condition is a factor with proper labels
nmds_scores$Condition <- factor(nmds_scores$Condition, 
                                levels = c("UB", "B"),
                                labels = c("Unbleached", "Bleached"))

# Add clade vectors
clade_fit <- envfit(nmds, symb_mat, permutations = 999)
vectors <- as.data.frame(scores(clade_fit, display = "vectors"))
vectors$clade <- rownames(vectors)

# Plot NMDS with arrows and jittered points
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Condition, shape = Species), 
             size = 3, 
             position = position_jitter(width = 0.02, height = 0.02)) +
  # Add clade arrows
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3,"cm")), color = "black") +
  geom_text(data = vectors, aes(x = NMDS1*1.1, y = NMDS2*1.1, label = clade),
            color = "black", fontface = "bold") +
  facet_wrap(~ Site) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  theme_classic() +
  labs(title = "NMDS of Symbiont Community Composition",
       subtitle = paste("Stress =", round(nmds$stress, 4)),
       color = "Bleaching Condition",
       shape = "Species") +
  theme(legend.position = "right")

  
  library(ggplot2)
library(vegan)
library(grid)  # for arrow()

# Make sure Condition is a factor with proper labels
nmds_scores$Condition <- factor(nmds_scores$Condition, 
                                levels = c("UB", "B"),
                                labels = c("Unbleached", "Bleached"))

# Add clade vectors
clade_fit <- envfit(nmds, symb_mat, permutations = 999)
vectors <- as.data.frame(scores(clade_fit, display = "vectors"))
vectors$clade <- rownames(vectors)

library(ggplot2)
library(vegan)
library(grid)  # for arrow()

# Make sure Condition is a factor with proper labels
nmds_scores$Condition <- factor(nmds_scores$Condition, 
                                levels = c("UB", "B"),
                                labels = c("Unbleached", "Bleached"))

# Add clade vectors
clade_fit <- envfit(nmds, symb_mat, permutations = 999)
vectors <- as.data.frame(scores(clade_fit, display = "vectors"))
vectors$clade <- rownames(vectors)

# Plot NMDS with arrows and jittered points
# Make sure Condition is a factor with proper labels
nmds_scores$Condition <- factor(nmds_scores$Condition, 
                                levels = c("UB", "B"),
                                labels = c("Unbleached", "Bleached"))

# Add clade vectors
clade_fit <- envfit(nmds, symb_mat, permutations = 999)
vectors <- as.data.frame(scores(clade_fit, display = "vectors"))
vectors$clade <- rownames(vectors)

nmds_scores <- as.data.frame(scores(nmds))

nmds_scores$Condition <- symb$Condition
nmds_scores$Species <- symb$Species
nmds_scores$Site <- symb$Site

# Plot NMDS with arrows and jittered points
nmds_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Condition, shape = Species), 
             size = 3, 
             position = position_jitter(width = 0.08, height = 0.08)) +
  # Add clade arrows
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3,"cm")), color = "black") +
  geom_text(data = vectors, aes(x = NMDS1*1.1, y = NMDS2*1.1, label = clade),
            color = "black", fontface = "bold") +
  facet_wrap(~ Site) +
  scale_color_manual(values = c("steelblue", "lightpink")) +  # Blue = Unbleached, Pink = Bleached
  theme_classic() +
  labs(title = "NMDS of Symbiont Community Composition",
       subtitle = paste("Stress =", round(nmds$stress, 4)),
       color = "Bleaching Condition",
       shape = "Species") +
  theme(legend.position = "right")

nmds_plot

# Save to plots folder

ggsave(filename = "plots/NMDS_symbionts.png", 
       plot = nmds_plot, 
       width = 14, height = 7, dpi = 300)


##############################################################
#See which clades are dominant in each site, species, condition
##############################################################


#get summary data
symb_summary <- symb %>%
  group_by(Site, Species, Condition) %>%
  summarise(
    n = n(),
    mean_A = mean(Prop_A, na.rm = TRUE),
    sd_A = sd(Prop_A, na.rm = TRUE),
    mean_B = mean(Prop_B, na.rm = TRUE),
    sd_B = sd(Prop_B, na.rm = TRUE),
    mean_C = mean(Prop_C, na.rm = TRUE),
    sd_C = sd(Prop_C, na.rm = TRUE),
    mean_D = mean(Prop_D, na.rm = TRUE),
    sd_D = sd(Prop_D, na.rm = TRUE)
  )

#####################################
#Interpret
#####################################
#which clades are dominant
symb_summary_dom <- symb_summary %>%
  rowwise() %>%
  mutate(
    dominant_clade = c("A","B","C","D")[which.max(c(mean_A, mean_B, mean_C, mean_D))]
  )


#make a barplot


# Convert to long format
symb_long <- symb %>%
  pivot_longer(cols = starts_with("Prop_"), 
               names_to = "Clade", 
               values_to = "Prop") %>%
  mutate(Clade = gsub("Prop_", "", Clade),
         Condition = factor(Condition, levels = c("UB","B"), labels = c("Unbleached","Bleached")))

# Summarize mean proportion per group
symb_summary <- symb_long %>%
  group_by(Site, Species, Condition, Clade) %>%
  summarise(mean_prop = mean(Prop, na.rm = TRUE), .groups = "drop")

# Stacked bar plot with bleaching side-by-side
symb_prop <- ggplot(symb_summary, aes(x = Species, y = mean_prop, fill = Clade)) +
  geom_bar(stat = "identity", position = position_stack(reverse = FALSE)) +
  facet_grid(Site ~ Condition) +
  scale_fill_manual(values = c("A"="gold","B"="magenta","C"="pink","D"="pink4")) +
  theme_classic() +
  labs(y = "Average Proportion", x = "Species", fill = "Clade") +
  theme(strip.text = element_text(face = "bold"))


symb_prop

ggsave("../Thesis_2026/plots/symbionts_prop.png", symb_prop, width = 10, height = 12, dpi = 300)
ggsave("../Thesis_2026/plots/symbionts_prop.png", symb_prop, width = 10, height = 12, dpi = 300)
