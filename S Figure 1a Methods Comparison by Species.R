# ---- Script header ----


# Method Comparisons by Species
# Date of last edit: 16/01/2026
# Description: This code is meant to be ran altogether at once. For plotting differences in Rubisco dark 
  # dark inhibition across species that have data from all three methods used in literature for 
  # identifying dark inhibition values. 

#### Load packages 
library(ggplot2)
library(multcompView)
library(dplyr)
library(readxl)
library(agricolae)
library(tibble)
library(forcats)
library(car)
library(FSA)         # For dunnTest()
library(rcompanion)  # For cldList()
library(viridisLite)


# ---- data Loading and prep -----


#### Read the Excel file into a data frame
df <- read_excel("./NehlsRamos_etal_JxB_BriefComms_Code_and_Dataset/NehlsRamos_etal_JxB_BriefComms_Data.xlsx", 
                 sheet = "Full Dataset")

#filter for species that had multiple methods for dark inhibition measurement
df <- df %>%
  group_by(Species) %>%
  filter(n_distinct(Method) >= 2) %>%
  ungroup()

# Ensure order of species remains by phylogenetic proximity
df$Order <- factor(df$Order, levels = unique(df$Order))
df$Species <- factor(df$Species, levels = unique(df$Species))

#for clean x axis species labels, removing "_"
Species_Labels <- c("Alocasia macrorrhiza", "Triticum aestivum", "Megathyrsus maximus", "Zea mays", 
                    "Helianthus annuus", "Petunia hvbrida", "Solanum lycopersicum", "Beta vulgaris",  
                    "Spinacia oleracea", "Phaseolus vulgaris", "Arabidopsis thaliana")


# ---- Final Plot ----


ggplot(data = subset(df, !is.na(Species)), aes(x = Species, y = Dark_Inhibition)) +
  geom_hline(yintercept = c(0.18, 0.44, 0.77), 
             linetype = "dashed", color = "grey70") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6, size = 3, aes(x = Species, y = Dark_Inhibition, color = Method)) +  # Show individual points
  scale_color_viridis_d(option = "viridis") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, face = "italic"), 
        axis.title = element_text(size = 18, face = "bold"), legend.position = "bottom", legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"), plot.title = element_text(size = 18, face = "bold")) +
  scale_x_discrete(labels = Species_Labels) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Species", y = "Dark Inhibition") +
  guides(colour = guide_legend(override.aes = list(size=4), title = "Method: "))


# ---- End ----

