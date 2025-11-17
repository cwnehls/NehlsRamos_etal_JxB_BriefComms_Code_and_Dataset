# ---- Script header ----


# KDE Modality Plot for Dark Inhibition Levels
# Date of last edit: 12/11/2025
# Description: To analyse the clustering of  Rubisco dark inhibition Flowering Plant data points across all 
  # Radiolabbeled CO2 assay Flowering Plant samples,within Fabales, and within samples not including Fabales.
  # NOTE: You must select your dataset by de-commenting the desired filter method below. Optionally you can
  # you can decomment the mclust modelling to help identify the predicted modes within the data for creating
  # the final kernel density estimate (KDE) plot.

#### Load packages 
library(readxl)
library(dplyr)
library(pracma)
library(diptest)
library(mclust)
library(multimode)


# ---- data Loading, selection and prep -----


#### Read the Excel file into a data frame
df <- read_excel("./NehlsRamos_etal_JxB_BriefComms_Code_and_Dataset/NehlsRamos_etal_JxB_BriefComms_Data.xlsx", 
                 sheet = "Full Dataset")

#### For the dataset with all Radiolabelled CO2 data, remove the #'s from below
df_filtered <- df %>%
  filter(Method == "Radiolabelled CO2", !Order %in% c("Chlamydomonadales", "Sphaeropleales"))

#### For only Fabales data, remove the #'s from below
#df_filtered <- df %>%
  #filter(Method == "Radiolabelled CO2", Order == "Fabales")

#### For the dataset without Fabales data remove the #'s from below
#df_filtered <- df %>%
  filter(Method == "Radiolabelled CO2", !Order %in% c("Chlamydomonadales", "Sphaeropleales", "Fabales"))

# Use only averaged species values
Averages <- df_filtered %>%
  group_by(Species) %>%
  summarise(avg_value = mean(Dark_Inhibition, na.rm = TRUE))


# ----- (Optional) McLust Modelling for Modes and Peaks ----


#Remove the comments below to perform McLust Modeling. NOTE: this section should be ran by itself, after loading data

#McLust modeling for number of modes; remove #'s below to use
#model <- Mclust(df_filtered$Dark_Inhibition)
#summary(model)
#plot(model)


# ---- Prepare  plot ----


# Kernel Density Estimate (KDE)
density_est <- density(Averages$avg_value, bw = 0.05669) #NOTE: set bw to the critical bandwidth 

# Plot KDE
plot(density_est, main = "Kernel Density Estimation with Peaks")
rug(Averages$avg_value)

# find peaks
peaks <- findpeaks(density_est$y)
peak_locations <- density_est$x[peaks[,2]]

# Add peaks to plot
points(peak_locations, peaks[,1], col = "orange", pch = 19)

#Hartigan's Dip Test:
  # If the p-value is small (< 0.05), the null hypothesis of unimodality is rejected, it is likely multimodal.
dip.test(df_filtered$Dark_Inhibition)

# locmode
result <- locmodes(Averages$avg_value, mod0 = 6, display = TRUE) #NOTE: set mod0 as needed, comparing with McLust results for expected number of peaks / modes

# bandwidth
bw_val <- as.numeric(result$cbw[3])

# show density bands
d <- density(Averages$avg_value, bw = bw_val)
plot(d, main = "", xlab = "Dark Inhibition", ylab = "Density", xlim = c(-1,1), ylim = c(0,2.5))
rug(Averages$avg_value)


# Find local minima (antimodes)
den <- d$y
x_vals <- d$x
antimode_indices <- which(diff(sign(diff(den))) == 2) + 1
antimodes_found <- x_vals[antimode_indices]
antimodes_y     <- den[antimode_indices]  # y-values of minima

# Add points at antimodes
points(antimodes_found, antimodes_y, col = "orange", pch = 19)

# Label antimodes directly on dips
text(antimodes_found, antimodes_y, 
     labels = round(antimodes_found, 2),
     pos = 3, cex = 0.8, col = "black")

# prepare peak data and labels
peaks <- findpeaks(d$y)
peak_locations <- d$x[peaks[,2]]

text(peak_locations, peaks[,1], 
     labels = round(peak_locations, 2),
     pos = 1, cex = 0.8, col = "darkgreen")

# Add peaks to plot
points(peak_locations, peaks[,1], col = "green", pch = 19)

# Vertical dotted lines from x-axis to curve height
segments(x0 = antimodes_found, y0 = 0, 
         x1 = antimodes_found, y1 = antimodes_y, 
         col = "black", lty = 2)


# ---- End ----
