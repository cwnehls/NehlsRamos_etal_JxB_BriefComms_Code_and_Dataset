# ---- Script header ----


# Kruskal Wallis Analyses for Tribes with 5 or More Data Points
# Date of last edit: 16/01/2026
# Description: This code is meant to be ran altogether at once. This code is made to analyse 
  # the variation of Rubisco dark inhibition across flowering plant tribes. Optional 
  # tests for normality and variation are included. This code uses the Kruskal-Wallis paired  
  # with Dunn test and compact letter display (cld) to produce a plot of the variation of
  # Rubisco dark inhibition, colored by order. only cld results and box plots are shown
  # for orders with 5 or more data points. For all others only the median is shown. tribes 
  # are ordered by phylogenetic proximity, with the more basal on the left

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

# Filter out for only Rabiolabbeled CO2 method and flowering plants
df_filtered <- df %>%
  filter(Method == "Radiolabelled CO2", !Tribe == "NA")

# Filter for tribes with 5 or more data points
df_5plus <- df_filtered %>%
  group_by(Tribe) %>%
  filter(n() >= 5, ) %>%
  ungroup()


# ---- Build model and test for normality ----


# The commented code below is optional and is used for determining the normality and variation of the data

#model <- aov(Dark_Inhibition ~ Tribe, data = df_filtered)

# Histogram
#hist(residuals(model))

# Q-Q plot
#qqnorm(residuals(model)); qqline(residuals(model))

# Shapiro-Wilk test
#shapiro.test(residuals(model))

# Levene’s Test 
#leveneTest(Dark_Inhibition ~ Tribe, data = df_filtered)

# Perform Kruskal wallis test to indentify significant differences between Tribe median dark inhibition levels
kruskal.test(Dark_Inhibition ~ Tribe, data = df_5plus)

# Perform Dunn test on Kruskal-wallis results to compare significances between Tribe
dunn_results <- dunnTest(Dark_Inhibition ~ Tribe, data = df_5plus, method = "bonferroni")

# Extract the p-values
dunn_df <- dunn_results$res

# Create the CLD to represent Tribe differences
cld <- cldList(P.adj ~ Comparison,
               data = dunn_df,
               threshold = 0.05)


# ---- Prepare data for plotting ----


# Make summary dataframe
df_summary <- df_5plus %>%
  group_by(Tribe) %>%
  summarise(
    mean = mean(Dark_Inhibition),
    sd = sd(Dark_Inhibition),
    n = n(),
    se = sd / sqrt(n)) %>%
  left_join(cld, by = c("Tribe" = "Group"))

# Prepare CLD label positions
df_labels <- df_5plus %>%
  group_by(Tribe) %>%
  summarise(Q4 = quantile(Dark_Inhibition, 1),
            SD = sd(Dark_Inhibition), .groups = "drop") %>%
  left_join(cld, by = c("Tribe" = "Group")) %>%
  mutate(
    cld_y = Q4 + 0.025,                             
    sd_y  = Q4 + 0.075                        
  )

df_5plus$Order <- factor(df_5plus$Order, levels = unique(df_5plus$Order))
df_5plus$Tribe <- factor(df_5plus$Tribe, levels = unique(df_5plus$Tribe))


# ---- Final Plot ----


ggplot(data = subset(df_5plus, !is.na(Tribe)), aes(x = Tribe, y = Dark_Inhibition)) +
  geom_hline(yintercept = c(0.18, 0.44, 0.77), 
             linetype = "dashed", color = "grey70") +
  geom_boxplot(outlier.shape = NA, fill = "white") +  # Don't plot outliers twice
  geom_jitter(width = 0.15, height = 0, alpha = 0.6, aes(x = Tribe, y = Dark_Inhibition, color = Order)) +  # Show individual points
  geom_text(data = df_labels, aes(x = Tribe, y = cld_y, label = Letter),
            size = 5, vjust = 0) +
  scale_color_viridis_d(option = "turbo") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, face = "italic"), 
        axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"), plot.title = element_text(size = 16.5, face = "bold")) +
  labs(title = "",
       x = "Tribe", y = "Dark Inhibition") +
  guides(colour = guide_legend(override.aes = list(size=4), title = "Order: "))


# ---- End ----

