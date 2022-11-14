#!/usr/bin/Rscript

###############################################################
#       Programm producing boxplots of normalized DTWm        #
# Author: Alice Laigle - M2BB - Nantes University (2021/2022) #
###############################################################

# Load libraries
library(readr)
library(ggplot2) 
theme_set(theme_bw()) # Define theme for plots

# Load data
distri_0_25_0_05 <- read_csv("./analysis_anonym_meth1_prop-level_0.25_perturb-level_0.05/distri_dissim_norm_meth1_prop-level_0.25_perturb-level_0.05.csv")
distri_0_25_0_25 <- read_csv("./analysis_anonym_meth1_prop-level_0.25_perturb-level_0.25/distri_dissim_norm_meth1_prop-level_0.25_perturb-level_0.25.csv")
distri_0_50_0_05 <- read_csv("./analysis_anonym_meth1_prop-level_0.50_perturb-level_0.05/distri_dissim_norm_meth1_prop-level_0.50_perturb-level_0.05.csv")
distri_0_50_0_25 <- read_csv("./analysis_anonym_meth1_prop-level_0.50_perturb-level_0.25/distri_dissim_norm_meth1_prop-level_0.50_perturb-level_0.25.csv")

# Make a list of the 4 normalized dissimilarity lists 
lists_distri <- list(x1 = distri_0_25_0_05, x2 = distri_0_25_0_25, x3 = distri_0_50_0_05, x4 = distri_0_50_0_25)
lists_distri

# Make a data.frame of these 4 normalized dissimilarity lists 
df_distri <- data.frame(lists_distri = unlist(lists_distri), Parameters = rep(c("0.25 - 0.005","0.25 - 0.25", "0.50 - 0.005","0.50 - 0.25")))

## Boxplot ##
boxplot <- ggplot(df_distri, aes(x = Parameters, y = lists_distri, color = Parameters)) + geom_boxplot() 

# Change labels 
black.text <- element_text(color = "black", size = "10")
black.bold.text <- element_text(face="bold", color = "black")
boxplot2 <- boxplot + theme(title=black.bold.text, axis.title=black.text) + labs(title="Normalized dissimilarities between anonymized and real patients",
                 x="Parameters (Proportion level - Perturbation level)", y = "Dissimilarity values")

#Â Save the file 
ggsave(boxplot2, file="./tools/boxplot_meth1.png", width=9, height=7)
