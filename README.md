extr### Alienz Insurance
## 2026 SOA Case Study: Actuaries in Space - The Final Frontier
By: Amelia Chung, Asrith Devarapalli, Ho Yin Lam, Daniel Song, Nathan Tan

# Objective Overview
As part of the collaboration between Galaxy General Insurance Company and Cosmic Quarry Mining Company, our team has been tasked with developing an insurance pricing strategy to best cover the risks of four areas (Business Interruptions, Cargo Loss, Workers' Compensation and Equipment Failure) in three different solar systems: Helionis Cluster, Bayesia System and Oryn Delta.

# Libraries, Data Cleaning, Data Limitations
For our analysis, we utilised the following libraries:
```r
library(readxl)
library(stringr)
library(dplyr)
library(ggplot2)
library(fitdistrplus)
library(MASS)
library(caret)
library(regclass)
library(glmmTMB)
library(DHARMa)
library(gap)
```

The raw historical data used as the basis to price our product contained various data integrity issues such as inconsistent values and wrong data classes, thus cleaning was conducted for each dataset. Following the provided Data Dictionary, each column was checked to ascertain it aligned with the required format and otherwise corrected with the following code snippet:
```r
### Functions to help optimise the code

# Clean issues in IDs. For example, "BI-000093_???5643" -> "BI-000093"
clean_id <- function(x) str_replace(x, "_[?][?][?].*", "") |> str_trim()

#For numeric columns: values that are clearly out of the defined range, are set to NA
#The specific low, high thresholds are specified later when the function is used.
clean_numeric <- function(x, lo = -Inf, hi = Inf) {x[!is.na(x) & (x < lo | x > hi)] <- NA
  x}

#Round integers, and see if they are in the valid levels, or else set as NA
clean_levels <- function(x, valid = 1:5) {
  x_int <- round(as.numeric(x))
  x_int[!x_int %in% valid] <- NA
  x_int}

# Clean the "???" issues in character columns, e.g. "Epsilon_???1063" -> "Epsilon"
clean_charac <- function(x) str_replace(x, "_[?][?][?].*", "") |> str_trim()
```

To deal with NA factors affecting frequency model training, if the proportion of NA values in the dataset was high, median imputation would be applied to the factor and character class data to best maintain central tendency and ensure useful data would not be lost from omission. However, for low proportion of NAs where bias risk would be minimal and imputation would cause reduced variance and distort the distribution, NA value omission was applied instead: 
```r
clean_data <- function(df) {
  char_cols <- sapply(df, is.character)
  df[char_cols] <- lapply(df[char_cols], factor)
  
  factor_cols <- sapply(df, is.factor)
  for (col in names(df)[factor_cols]) {
    mode_val <- names(sort(table(df[[col]]), decreasing = TRUE))[1]
    df[[col]][is.na(df[[col]])] <- mode_val
  } # impute NA factors and characters (numeric handled by caret)
  
  df_clean <- df[!is.na(df$claim_count), ] # remove NA values from target variable
  return(df_clean)
}
```

Moving onto the data limitations, a primary limitation was the historical dataset mismatch of solar systems, with claims data being available for Helionis, Epsilon and Zeta instead of Helionis, Bayesia and Oryn Delta where future operations would occur. To overcome this, a proxy calibration approach was used along with the provided information on each system in the Encyclopedia. Despite the structural similarity however, this approach introduces parameter uncertainty. Therefore, sensitivity testing was done to assess the impact of proxy misspecification on aggregate loss tails and capital requirements. Other limitations included:
- Limited data on extremely correlated solar storm events
- Potential inflation misalignment over long term transit windows in Oryn Delta
- Limited granularity in human capital risk variables

# Product Design

# Aggregate Loss Modelling

# Risk Assessment

Use these below in order to link and display different formats:
[data](player_data_salaries_2020.csv), [code](sample-data-clean.ipynb), and [images](ACC.png) [link](www.google.com)

This file is written using Markdown.
