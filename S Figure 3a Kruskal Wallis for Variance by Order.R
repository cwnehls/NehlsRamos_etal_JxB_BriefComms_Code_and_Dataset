# ---- Script header ----


# Kruskal Wallis for Variance by Order
# Date of last edit: 18/01/2026
# Description: This code is meant to be ran altogther at once. This code is made to analyse 
  # the variation of Rubisco dark inhibition across flowering plant orders. Optional 
  # tests for normality and variation are included. This code uses the Kruskal-Wallis paired  
  # with Dunn test and compact letter display (cld) to produce a plot of the variation of
  # Rubisco dark inhibition, colored by dark inhibition level. only cld results and box 
  # plots are shown for orders with 6 or more data points. For all others only the median is 
  # shown. Orders are ordered by phylogenetic proximity, with the more basal genera on the left

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


# ---- data Loading and prep -----


#### Read the Excel file into a data frame
df <- read_excel("./NehlsRamos_etal_JxB_BriefComms_Code_and_Dataset/NehlsRamos_etal_JxB_BriefComms_Data.xlsx", 
                 sheet = "Full Dataset")

# Filter out for only Rabiolabbeled CO2 method and flower plants
df_filtered <- df %>%
  filter(Method == "Radiolabelled CO2", !Order %in% c("Sphaeropleales", "Chlamydomonadales"))

df_filtered$Order <- factor(df_filtered$Order, levels = unique(df_filtered$Order))

# Select for only averaged inhibition values per species
Averages <- df_filtered %>%
  group_by(Order, Species) %>%
  summarise(avg_value = mean(Dark_Inhibition), .groups = "drop")

# Filter for Orders with 6 or more data points for boxplots and further analyses
df_6plus <- Averages %>%
  group_by(Order) %>%
  filter(n() >= 6) %>%
  ungroup()


# ---- Build model and test for normality ----


# The commented code below is optional and is used for determining the normality and variation of the data

# model <- aov(avg_value ~ Order, data = df_6plus)

# Histogram
#hist(residuals(model))

# Q-Q plot
#qqnorm(residuals(model)); qqline(residuals(model))

# Shapiro-Wilk test
#shapiro.test(residuals(model))

# Levene’s Test 
#leveneTest(avg_value ~ Order, data = df_6plus)

# Perform Kruskal wallis test to indentify significant differences between order median dark inhibition levels
kruskal.test(avg_value ~ Order, data = df_6plus)

# Perform Dunn test on Kruskal-wallis results to compare significances between orders
dunn_results <- dunnTest(avg_value ~ Order, data = df_6plus, method = "bonferroni")

# Show results of Kruskal Wallis paired with Dunn test
print(dunn_results$res)

# Extract the p-values
dunn_df <- dunn_results$res

# Create the CLD to represent order differences
cld <- cldList(P.adj ~ Comparison,
               data = dunn_df,
               threshold = 0.05)


# ---- Prepare data for plotting ----

# Make summary dataframe
df_summary <- df_6plus %>%
  group_by(Order) %>%
  summarise(
    mean = mean(avg_value),
    sd = sd(avg_value),
    n = n(),
    se = sd / sqrt(n)) %>%
  left_join(cld, by = c("Order" = "Group"))

#identify placement for letters
df_q4 <- Averages %>%
  group_by(Order) %>%
  summarise(Q4 = quantile(avg_value, 1))

df_letters <- df_q4 %>%
  left_join(cld, by = c("Order" = "Group")) %>%
  mutate(y_pos = Q4 + 0.1)  # Adjust vertical position slightly

# Separate data by levels of dark inhibition
df_colored <- Averages %>%
  mutate(threshold_group = case_when(
    avg_value < 0.18 ~ "Low (<18%)",
    avg_value >= 0.18 & avg_value < 0.44 ~ "Moderate (18% - 44%)",
    avg_value >= 0.44 & avg_value < 0.77 ~ "High (44% - 77%)",
    avg_value >= 0.77 ~ "Very High (>77%)"
  ))

df_colored$threshold_group <- factor(
  df_colored$threshold_group,
  levels = c( "Very High (>77%)", "High (44% - 77%)", "Moderate (18% - 44%)", "Low (<18%)")
)

# Filter for groups with 6 or more datapoints that can have boxplots and further analyses
df_box <- df_colored %>%
  group_by(Order) %>%
  filter(n() >= 6) %>%
  ungroup()

# Filter for all other groups so as to only show medians in plot aside from data points
df_small <- df_colored %>%
  group_by(Order) %>%
  filter(n() < 6) %>%
  summarise(mean_value = mean(avg_value), .groups = "drop")

# Ensure data retains original order by phylogeny proximity
original_order <- unique(df_colored$Order)

df_colored$Order <- factor(df_colored$Order, levels = original_order)
df_box$Order     <- factor(df_box$Order, levels = original_order)
df_small$Order   <- factor(df_small$Order, levels = original_order)


# ---- Final Plot ----


ggplot() +
  geom_hline(yintercept = c(0.18, 0.44, 0.77), 
             linetype = "dashed", color = "grey70") +
  geom_jitter(data = df_colored, width = 0.175, alpha = 0.6, size = 2.5, aes(x = Order, y = avg_value, color = threshold_group)) +  # Show individual points
  geom_boxplot(data = df_box, aes(x = Order, y = avg_value), 
               outlier.shape = NA, fill = "transparent") +  # Don't plot outliers twice
  geom_crossbar(data = df_small,
                aes(x = Order, y = mean_value, ymin = mean_value, ymax = mean_value),
                width = 0.4, color = "black", fatten = 2) + # Show individual points
  scale_color_manual(name = "Dark Inhibition Levels", values = c(
    "Low (<18%)" = "yellow2",
    "Moderate (18% - 44%)" = "#96d831",
    "High (44% - 77%)" = "#32a354",
    "Very High (>77%)" = "#3b4bcd")) +  
  geom_text(data = df_letters, aes(x = Order, y = y_pos, label = Letter),
            size = 6, vjust = 0, na.rm = TRUE) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), 
        axis.title = element_text(size = 22, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text  = element_text(size = 13), 
        panel.grid.minor.y = element_blank(), legend.position ="bottom") +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5), reverse = TRUE)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position="bottom") +
  labs(title = "",
       x = "Order", y = "Dark Inhibition")


# ---- End ----
