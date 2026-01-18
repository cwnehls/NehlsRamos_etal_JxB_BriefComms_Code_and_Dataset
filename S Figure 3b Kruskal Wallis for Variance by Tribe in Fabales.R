# ---- Script header ----


# Kruskal Wallis for Variance by Tribe in Fabales
# Date of last edit: 18/01/2026
# Description: This code is meant to be ran altogether at once. This code is made to analyse 
  # the variation of Rubisco dark inhibition within data for the order Fabales. Optional 
  # tests for normality and variation are included. This code uses the Kruskal-Wallis paired  
  # with Dunn test and compact letter display (cld) to produce a plot of the variation of
  # Rubisco dark inhibition across Fabales genera, colored by tribe. only cld results and box 
  # plots are shown for genera with 6 or more data points. For all others only the median is 
  # shown. Genera are ordered by phylogenetic proximity, with the more basal genera on the left

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

# Filter data to only include Radiolabelled CO2 in Fabales

df_filtered <- df %>%
  filter(Method == "Radiolabelled CO2", Order == "Fabales")

df_filtered$Genus <- factor(df_filtered$Genus, levels = unique(df_filtered$Genus))

df_6plus <- df_filtered %>%
  group_by(Genus) %>%
  filter(n() >= 6) %>%
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

# Perform Kruskal wallis test to indentify significant differences between genera median dark inhibition levels
kruskal.test(Dark_Inhibition ~ Genus, data = df_6plus)

# Perform Dunn test on Kruskal-wallis results to compare significances between genera
dunn_results <- dunnTest(Dark_Inhibition ~ Genus, data = df_6plus, method = "bonferroni")

# Show results of Kruskal Wallis paired with Dunn test
print(dunn_results$res)

# Extract the p-values from dunn test
dunn_df <- dunn_results$res

# Create the CLD to represent genera differences
cld <- cldList(P.adj ~ Comparison,
               data = dunn_df,
               threshold = 0.05)


# ---- Prepare data for plotting ----


# Make summary dataframe
df_summary <- df_filtered %>%
  group_by(Genus) %>%
  summarise(
    mean = mean(Dark_Inhibition),
    sd = sd(Dark_Inhibition),
    n = n(),
    se = sd / sqrt(n)) %>%
  left_join(cld, by = c("Genus" = "Group"))


# Make dataframe of labels for graph. 
df_labels <- df_filtered %>%
  group_by(Genus, Tribe) %>%
  summarise(Q4 = quantile(Dark_Inhibition, 1),
            SD = sd(Dark_Inhibition), .groups = "drop") %>%
  left_join(cld, by = c("Genus" = "Group")) %>%
  mutate(
    cld_y = Q4 + 0.025,                             
    sd_y  = Q4 + 0.075                        
  )

# filter for significant groups
significant_groups <- unique(c(dunn_results$group1[dunn_results$p.adj < 0.05],
                               dunn_results$group2[dunn_results$p.adj < 0.05]))

# Add significance to labels dataframe
df_labels_filtered <- df_labels %>%
  filter(Genus %in% significant_groups)

# Filter for groups with 6 or more datapoints that can have boxplots and further analyses
df_box <- df_filtered %>%
  group_by(Genus) %>%
  filter(n() >= 6) %>%
  ungroup()

# Filter for all other groups so as to only show medians in plot aside from data points
df_small <- df_filtered %>%
  group_by(Genus) %>%
  filter(n() < 6) %>%
  summarise(mean_value = mean(Dark_Inhibition), .groups = "drop")

# Ensure the order of Genera are by phylogeny proximity, as in spreadsheet
original_order <- unique(df_filtered$Genus)

df_filtered$Genus <- factor(df_filtered$Genus, levels = original_order)
df_box$Genus     <- factor(df_box$Genus, levels = original_order)
df_small$Genus   <- factor(df_small$Genus, levels = original_order)

# Ensure the order of Tribes are by phylogeny proximity, as in spreadsheet
Tribe_order <- unique(df_filtered$Tribe) 

# Generate turbo colors by tribe
n_groups <- length(Tribe_order)
colors <- viridisLite::turbo(n_groups)

# Stagger first half and second half of colors
staggered_colors <- colors[c(seq(1, n_groups, by = 2), seq(2, n_groups, by = 2))]

# Name the colors by the original group order
names(staggered_colors) <- Tribe_order

df_filtered$Tribe <- factor(df_filtered$Tribe, levels = Tribe_order)


# ---- Final Plot ----


ggplot() +
  geom_hline(yintercept = c(0.18, 0.44, 0.77), 
             linetype = "dashed", color = "grey70") +
  geom_jitter(data = df_filtered, size = 2, aes(x = Genus, y = Dark_Inhibition, color = Tribe), width = 0.15, alpha = 0.6) +
  geom_boxplot(data = df_box, aes(x = Genus, y = Dark_Inhibition), 
               outlier.shape = NA, fill = "transparent") +  # Don't plot outliers twice
  geom_crossbar(data = df_small,
                aes(x = Genus, y = mean_value, ymin = mean_value, ymax = mean_value),
                width = 0.4, color = "black", fatten = 2) + # Show individual points
  geom_text(data = df_labels, aes(x = Genus, y = cld_y, label = Letter),
            size = 6, vjust = 0, na.rm = TRUE) +
  scale_color_viridis_d(option = "turbo") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, face = "italic"), 
        axis.title = element_text(size = 22, face = "bold"), legend.position = "bottom", legend.text = element_text(size = 13, face = "italic"),
        legend.title = element_text(size = 16, face = "bold"), plot.title = element_text(size = 18, face = "bold")) +
  labs(title = "",
       x = "Genus", y = "Dark Inhibition")


# ---- End ----


