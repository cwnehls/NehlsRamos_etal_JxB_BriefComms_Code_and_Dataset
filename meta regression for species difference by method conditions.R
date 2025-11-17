# ---- Script header ----


# Meta-Regression for Species Differences within the Radiolabbeled CO2 Assay Data
# Date of last edit: 12/11/2025
# Description: This code is meant to be ran altogether at once. This code is made to analyse 
  # for significant differences between assay and sampling conditions within the dark inhibition 
  # data for the Radiolabbeled CO2 assay method. The Analyses factor in the natural variation 
  # found with species differences as well as by study source. Results are reported within the console.
  # Included as well is a robust variance estimation in addition, factor for species differences further.

#### Load packages 
library(ggplot2)
library(dplyr)
library(readxl)
library(metafor)
library(clubSandwich)
library(robumeta)


# ---- data Loading and prep -----


#### Read the Excel file into a data frame
df <- read_excel("./NehlsRamos_etal_JxB_BriefComms_Code_and_Dataset/NehlsRamos_etal_JxB_BriefComms_Data.xlsx", 
                 sheet = "Full Dataset")

### Filter out for only Rabiolabbeled CO2 method
df <- df %>%
  filter(Method %in% c("Radiolabelled CO2"))

# assume (covariance) CV = 0.2 for enzyme activity measures
assumed_cv <- 0.2

# approximate variance using CV
df <- df %>%
  mutate(
    vi = (assumed_cv^2)
  )

# Convert numeric data to numeric
df <- df %>%
  mutate(across(c( Light_Level_umol_m2_s, Dark_Level_umol_m2_s, Rubisco_Activation_Time_min,
                    Rubisco_Activity_Time_min, Buffer_mM, MgCl2_mM, DTT_mM,
                    B_mercaptoethanol_mM, EDTA_mM, RuBP_mM, NaH14CO3_20mM, KHCO3_mM, 
                    ATP_mM, Phosphocreatine_mM, NADH_mM, pH, Temperature_C), as.numeric))

# Convert non-numeric data as factors
df <- df %>%
  mutate(across(c(Activity_Measurement, Method, Buffer_Type), as.factor))

yi = df$Dark_Inhibition


# ---- Build model and test for normality ----


# Fit regression model
res.ml <- rma.mv(yi, vi, mods = ~ Light_Level_umol_m2_s + Dark_Level_umol_m2_s + Rubisco_Activation_Time_min +
                    Rubisco_Activity_Time_min + Buffer_Type + Buffer_mM + MgCl2_mM + DTT_mM +
                    B_mercaptoethanol_mM + EDTA_mM + RuBP_mM + NaH14CO3_20mM + pH + Temperature_C, random = ~1 | Species/Lit_Source, data = df, method="REML")

summary(res.ml)

# Robust variance estimation factoring for species heterogeneity
res_cv_robust <- coef_test(res.ml, vcov = "CR2", cluster = df$Species)

# print results of robust model
res_cv_robust


# ---- End ----