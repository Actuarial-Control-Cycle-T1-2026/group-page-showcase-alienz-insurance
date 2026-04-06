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
To effectively managage the vastly different risk profiles across Cosmic Quarry's operations, the portfolio is bifurcated into two layers: Tier 1 (Attritional) and Tier 2 (Catastrophic)
- **Tier 1 (Attritional Risk Layer)**: Covers Workers' Compensation (WC) and Equipment Failure (EF), tailored for the high-frequency, low-severity claims characteristic of the Helionis Cluster's high-traffic and high-debris environment.
- **Tier 2 (Catastrophic Risk Layer)**: Covers Business Interruption (BI) and Cargo Loss (CL). This layer focuses on systemic, correlated and totl-loss threats in the Bayesia System and Oryn Delta, where extreme isolation demands unique parametric and agreed-value structures.

Premiums are driven by Generalised Linear Models (GLMs) that isolate predictive operational variables, allowing the rating engine to automatically adjust to the relaties of each system. 
| Product | Primary Rating Variables | System Weighting & Actuarial Justification |
| :--- | :--- | :--- |
| **Workers' Compensation (WC)** | `gravity_level`, `safety_training_index` | **Highest in Helionis (Frequency):** High gravity environments directly inflate the frequency and severity of musculoskeletal claims, demanding heavy premium loading. |
| **Equipment Failure (EF)** | `usage_intensity`, `solar_radiation` | **Highest in Bayesia (Tail Risk):** Radiation spikes cause instantaneous electronic degradation, shifting EF from a predictable wear-and-tear risk to a volatile, heavy-tailed exposure. |
| **Business Interruption (BI)** | `energy_backup_score`, `supply_chain_index` | **Highest in Bayesia/Oryn (Correlation):** Poor backup scores heavily penalise premiums, as a single solar event can trigger simultaneous, system-wide communication and production stoppages. |
| **Cargo Loss (CL)** | `transit_duration`, `vessel_age`, `pilot_exp` | **Highest in Oryn Delta (Severity):** The 240 AU distance and 60-month duration maximise the probability of a total loss, requiring peak risk margins to cover extreme uncertainty. |

To implement this operational design in our model, we established a central sys_params matrix to map environmental constraints directly into the model: 
```r
# Parameters retained from Online Encyclopedia:
# route_risk : Helionis=3, Bayesia=2, Oryn=4
# debris_density : consistent with route_risk ordering
# solar_radiation: Helionis=0.50; Bayesia=0.70; Oryn=0.30

sys_params <- data.frame( cq_system = cq_systems, proxy_sys = unname(sys_proxy),
# WC: gravity
gravity      = c(1.125, 1.300, 1.113),
    
 # WC: psych_stress: training-data median per proxy system
 psych_stress = as.numeric(get_wc_param("psych_stress_med")),
    
# Cargo: route_risk 
route_risk_n = c(3L, 2L, 4L),   
route_risk_c = c("3","2","4"),   
    
 # Cargo: debris_density and solar_radiation
debris_density  = c(0.60, 0.30, 0.70),
solar_radiation = c(0.50, 0.70, 0.30),
    
# ... [Additional BI and Cargo Parameters] ...
 stringsAsFactors = FALSE
  )

# Pricing & Commerical Strategy 
The pricing framework is designed to capture the baseline expected 10-year present value cost while aggressively funding the liquid capital reserves required to survive extreme tail-risk volatility. 

The fundamental pricing equation applied across the portfolio is:
         **Gross Premium = Expected Loss + Cost of Capital + Expenses + Profit Margin**

Because CL represents an extreme capital burden (generating a standard deviation of 38.02 Billion Ð), a uniform pricing approach would critically undercapitalise Galaxy General. Therefore, we propose a **Modular Pricing Structure:**
1. **Core Coverage Package:**Includes BI, WC, and EF. Because these hazards exhibit stable finanical variance, they are priced with a predictable, moderate cost-of-capital margin.
2. **Cargo Loss Coverage (Modular Addition):**Priced separately as an optional or standalone capital module due to its massive tail risk. To optimise costs, this module utilises a layered risk-retention approach
   * Deductible: Policyholders retain an initial 34,000K Ð (5% of max cargo value) to absorb high-frequency attritional damage.
   * Primary Retention: Galaxy General covers severity up to the 680,000K Ð maximum fleet exposure limit.
   * Risk Transfer: Extreme total-loss events are partially transferred via CAT XL and Facultative reinsurance.

**10-year Comprehensive Pricing Structure**
Applying this framework to the fully underwritten comrpehensive portfolio (Core + Cargo) over a 10-year projection yields the following structure:
| Pricing Component | Value (10-Year PV) | Actuarial Rationale |
| :--- | :--- | :--- |
| **Expected Loss** | 90.21 Billion Ð | Baseline 10-year present value cost derived from aggregate Monte Carlo simulations. |
| **Cost of Capital (20%)** | 18.04 Billion Ð | Risk margin dedicated to funding reserves against the 197.97 Billion Ð $VaR_{0.99}$ threshold. |
| **Capital-Adjusted Pure Premium** | 108.25 Billion Ð | The true cost of absorbing Cosmic Quarry's operational risk. |
| **Target Gross Premium** | 166.54 Billion Ð | Factors in a 30% expense ratio (administration/reinsurance) and a 5% corporate profit margin (65% permissible loss ratio). |

**Long-Term Projection Implementation**
To simulate these commerical returns, we projected the portfolio out 10 years, accounting for compounding exposure growth, interplantary inflation and investment float yields. 
# ======
  # PRICING, EXPENSE, AND ECONOMIC PARAMETERS
  # ======
  expense_ratio <- 0.30   # 30% of premium for admin, reinsurance, commissions
  profit_margin <- 0.05   # 5% target profit margin
  load_factor   <- 1 / (1 - expense_ratio - profit_margin)  # Yields 1.538 (65% permissible loss ratio)
  
  # Economic Assumptions (5-yr trailing avg 2170–2174)
  inf_fwd   <- 0.0423     # 4.23% Claims inflation rate
  r1yr_fwd  <- 0.0371     # 3.71% Short-term reserve investment rate
  r10yr_fwd <- 0.0376     # 3.76% Long-term discount rate

  # ======
  # 10-YEAR FINANCIAL PROJECTION LOOP
  # ======
  for (t in seq_len(N_YEARS)) {
    
    exp_growth_t  <- (1 + wtd_growth)^(t - 1)    # Exposure multiplier
    inf_factor_t  <- (1 + inf_fwd)^(t - 1)       # Severity inflation multiplier
    disc_factor_t <- 1 / (1 + r10yr_fwd)^t       # PV discount factor
    
    # 1. Expected total loss for year t 
    E_S_t <- total_el_y1 * exp_growth_t * inf_factor_t
    
    # 2. Gross Premium (priced at start of year based on expected loss and load factor)
    P_t   <- E_S_t * load_factor
    
    # 3. Expenses (Deterministic % of premium)
    Exp_t <- P_t * expense_ratio
    
    # 4. Investment Income: 0.5 × E[loss] held as float, earned at short-term rate
    Inv_t <- 0.5 * E_S_t * r1yr_fwd
    
    # --- Simulate aggregate losses (S_t) for each line ---
    # [Simulation execution functions omitted for brevity]
    S_t  <- sim_t_bi + sim_t_wc + sim_t_ef + sim_t_cargo
    
    # 5. Net Revenue Calculation (per simulation)
    NR_t <- P_t + Inv_t - S_t - Exp_t
    
    # Store Present Values
    proj_pv_cost <- proj_pv_cost + S_t  * disc_factor_t
    proj_pv_rev  <- proj_pv_rev  + NR_t * disc_factor_t
  }
# Aggregate Loss Modelling

# Risk Assessment

Use these below in order to link and display different formats:
[data](player_data_salaries_2020.csv), [code](sample-data-clean.ipynb), and [images](ACC.png) [link](www.google.com)

This file is written using Markdown.
