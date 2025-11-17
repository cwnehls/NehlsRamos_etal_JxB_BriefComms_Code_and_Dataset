# ---- Script header ----


# Comparing Radiolabbled CO2 Methods
# Date of last edit: 12/11/2025
# Description: Code for analysing each of the components of the Rubisco Activity Assays and sample collection which
  # could impact Rubisco Dark Inhibition readings across plant species. NOTE: This code is meant to be ran at once,
  # and the desired condition to analyse must be de-commented prior to running. Additionally, as certain species only
  # have data for particular conditions due to differences in methods, the mean of the condition is used as a placeholder
  # to allow for the analyses to be performed.

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


# ---- data Loading, condition selection, and data prep -----


#### Read the Excel file into a data frame
df <- read_excel("./NehlsRamos_etal_JxB_BriefComms_Code_and_Dataset/NehlsRamos_etal_JxB_BriefComms_Data.xlsx", 
                 sheet = "Full Dataset")

### Filter out for only Rabiolabbeled CO2 method
df <- df %>%
  filter(Method %in% c("Radiolabelled CO2"))

####Filter for desired condition -- NOTE: remove # for the condition being analysed
#df <- df %>% rename_at('Rubisco_Activation_Time_min', ~'condition')
#df <- df %>% rename_at('Rubisco_Activity_Time_min', ~'condition')
#df <- df %>% rename_at('Buffer_Type', ~'condition')
#df <- df %>% rename_at('Buffer_mM', ~'condition')
#df <- df %>% rename_at('MgCl2_mM', ~'condition')
#df <- df %>% rename_at('DTT_mM', ~'condition')
#df <- df %>% rename_at('B_mercaptoethanol_mM', ~'condition')
#df <- df %>% rename_at('EDTA_mM', ~'condition')
#df <- df %>% rename_at('RuBP_mM', ~'condition')
#df <- df %>% rename_at('NaH14CO3_20mM', ~'condition')
#df <- df %>% rename_at('KHCO3_mM', ~'condition')
#df <- df %>% rename_at('ATP_mM', ~'condition')
#df <- df %>% rename_at('Phosphocreatine_mM', ~'condition')
#df <- df %>% rename_at('NADH_mM', ~'condition')
#df <- df %>% rename_at('pH', ~'condition')
#df <- df %>% rename_at('Temperature_C', ~'condition')

####Fill empty cells with average value of the column-- this is necessary for further analyses
# Replace NA with column mean
df <- df %>%
  mutate(condition = as.numeric(as.character(condition)),
         condition = ifelse(is.na(condition), mean(condition, na.rm = TRUE), condition))

#convert condition column to factor
df$condition <- as.factor(df$condition)


# ---- (optional) Test for normality ----


# The commented code below is optional and is used for determining the normality and variation of the data

####test assumptions to show Kruskal Wallis is the right test -- remove # from code as needed
#model <- aov(Dark_Inhibition ~ condition, data = df)

### Test for normality
## Histogram
#hist(residuals(model))
## Q-Q plot
#qqnorm(residuals(model)); qqline(residuals(model))
## Shapiro-Wilk test
#shapiro.test(residuals(model))

###Test for Equal Variance
## Levene’s Test 
#leveneTest(Dark_Inhibition ~ condition, data = df)


# ---- Kruskal-Walllis and Dunn Tests ----


####Kruskal-Wallis test 
kruskal.test(Dark_Inhibition ~ condition, data = df)

dunn_results <- dunnTest(Dark_Inhibition ~ condition, data = df, method = "bonferroni")

## Extract the p-values
dunn_df <- dunn_results$res

# Filter significant results (p-value < 0.05)
significant_results <- dunn_df[dunn_df$P.adj < 0.05, ]


# ---- Print results for significant treatments within a condition ----


#Check for significantly different condition treatments-- if mean treatment is shown, do not consider significant
print(significant_results)

# ---- End ----
