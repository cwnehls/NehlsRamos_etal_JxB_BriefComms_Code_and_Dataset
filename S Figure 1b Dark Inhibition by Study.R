# ---- Script header ----


# Dark Inhibition by Study
# Date of last edit: 12/11/2025
# Description: To plot and compare the Rubisco dark inhibition values of leaf / algae samples across different Orders and studies

#### Load packages 
library(ggplot2)
library(dplyr)
library(readxl)
library(tibble)
library(viridisLite)


# ---- data Loading and prep -----

 
#### Read the Excel file into a data frame
df <- read_excel("./NehlsRamos_etal_JxB_BriefComms_Code_and_Dataset/NehlsRamos_etal_JxB_BriefComms_Data.xlsx", 
                 sheet = "Full Dataset")

# Filter data to only include Radiolabelled CO2 assay samples
df_filtered <- df %>%
  filter(Method == "Radiolabelled CO2")

df_filtered$Order <- factor(df_filtered$Order, levels = unique(df_filtered$Order))

df_filtered$Species <- factor(df_filtered$Species, levels = unique(df_filtered$Species))

# Order studies by year of publication
df_filtered <- df_filtered %>%
  mutate(Lit_Source = reorder(Lit_Source, Year, FUN = mean))


# ---- Final Plot ----


ggplot() +
  geom_hline(yintercept = c(0.18, 0.44, 0.77), 
             linetype = "dashed", color = "grey70") +
  geom_jitter(data = df_filtered, aes(x = Lit_Source, y = Dark_Inhibition, color = Order),
              width = 0.2, alpha = 0.6) +
  scale_color_viridis_d(option = "turbo") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), 
        axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"), plot.title = element_text(size = 16.5, face = "bold")) +
  labs(title = "",
       x = "Study", y = "Dark Inhibition")


# ---- End ----
