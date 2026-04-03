  
  # Frequency + Severity + Aggregate Loss Modelling
  
  library(dplyr)
  library(ggplot2)
  library(fitdistrplus)
  library(MASS)
  library(caret)
  library(regclass)
  library(glmmTMB)
  library(DHARMa)
  library(gap)
  
  
  # =======
  # FREQUENCY MODELLING
  
  ### Datasets
  bi    <- read.csv("outputs/cleaned/bi_freq_clean.csv")
  cargo <- read.csv("outputs/cleaned/cargo_freq_clean.csv")
  ef    <- read.csv("outputs/cleaned/ef_freq_clean.csv")
  wc    <- read.csv("outputs/cleaned/wc_freq_clean.csv")
  
  ### Functions
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
  
  ### Models
  set.seed(5100)
  train_control <- trainControl(method = "cv", number = 10)
  
  ## Business Interruptions
  bi_clean <- na.omit(bi)
  bi_vars      <- setdiff(names(bi_clean), c("policy_id", "station_id", "claim_count", "X", "exposure"))
  joined_names <- paste(bi_vars, collapse = " + ")
  formula      <- as.formula(paste("claim_count", "~", joined_names, "+ offset(log(exposure))"))
  bi_poisson   <- train(formula,
                        data       = bi_clean,
                        method     = "glm",
                        family     = poisson,
                        trControl  = train_control)
  summary(bi_poisson) # only avg_crew_exp insignificant
  
  ## Cargo
  cargo_clean <- na.omit(cargo)
  cargo_vars   <- setdiff(names(cargo), c("policy_id", "shipment_id", "claim_count", "X", "exposure", "cargo_value"))
  joined_names <- paste(cargo_vars, collapse = " + ")
  formula      <- as.formula(paste("claim_count", "~", joined_names, "+ offset(log(exposure))"))
  cargo_nb     <- train(formula,
                        data      = cargo_clean,
                        method    = "glm.nb",
                        trControl = train_control,
                        tuneGrid  = expand.grid(link = "log"))
  summary(cargo_nb)
  
  ## Equipment Failure
  ef_clean <- clean_data(ef)
  ef_vars      <- setdiff(names(ef), c("policy_id", "equipment_id", "claim_count", "X", "exposure"))
  joined_names <- paste(ef_vars, collapse = " + ")
  formula      <- as.formula(paste("claim_count", "~", joined_names, "+ offset(log(exposure))"))
  ef_poisson   <- train(formula,
                        data       = ef_clean,
                        method     = "glm",
                        family     = poisson,
                        preProcess = c("medianImpute"),
                        na.action  = na.pass,
                        trControl  = train_control)
  summary(ef_poisson) # all significant
  
  ## Workers' Compensation
  wc_clean <- na.omit(wc)
  wc_vars      <- setdiff(names(wc), c("policy_id", "worker_id", "claim_count", "X", "exposure", "station_id"))
  joined_names <- paste(wc_vars, collapse = " + ")
  formula      <- as.formula(paste("claim_count", "~", joined_names, "+ offset(log(exposure))"))
  wc_poisson   <- train(formula,
                        data      = wc_clean,
                        method    = "glm",
                        family    = poisson,
                        trControl = train_control)
  summary(wc_poisson) # all significant except hours_per_week
  
  
  # =======
  # SEVERITY MODELLING
  
  bi_sev    <- read.csv("outputs/cleaned/bi_sev_clean.csv",    TRUE)
  cargo_sev <- read.csv("outputs/cleaned/cargo_sev_clean.csv", TRUE)
  ef_sev    <- read.csv("outputs/cleaned/ef_sev_clean.csv",    TRUE)
  wc_sev    <- read.csv("outputs/cleaned/wc_sev_clean.csv",    TRUE)
  
  # Cleaning
  bi_sev <- bi_sev %>%
    filter(solar_system != "NA", !is.na(solar_system),
           !is.na(claim_amount), claim_amount > 0) %>%
    mutate(claim_m    = claim_amount / 1e6,
           log_claim  = log(claim_amount))
  
  cargo_sev <- cargo_sev %>%
    filter(claim_amount > 0, !is.na(claim_amount)) %>%
    mutate(claim_amt_scaled = claim_amount / 1000000,
           container_type   = as.factor(container_type),
           cargo_type        = as.factor(cargo_type),
           route_risk        = as.factor(route_risk))
  
  ef_sev <- ef_sev %>%
    filter(claim_amount > 0, !is.na(claim_amount)) %>%
    mutate(claim_amt_scaled = claim_amount / 1000,
           equipment_type   = as.factor(equipment_type),
           solar_system     = as.factor(solar_system))
  
  wc_sev <- wc_sev %>%
    filter(claim_amount > 0, !is.na(claim_amount)) %>%
    mutate(claim_amt_scaled      = claim_amount / 1000,
           solar_system          = as.factor(solar_system),
           occupation            = as.factor(occupation),
           employment_type       = as.factor(employment_type),
           accident_history_flag = as.factor(accident_history_flag),
           psych_stress_index    = as.factor(psych_stress_index),
           safety_training_index = as.factor(safety_training_index),
           protective_gear_quality = as.factor(protective_gear_quality),
           injury_type           = as.factor(injury_type),
           injury_cause          = as.factor(injury_cause))
  
  ### Business Interruption: marginal lognormal (fitted on $millions)
  fit_bi_ln <- fitdist(bi_sev$claim_m, "lnorm")
  summary(fit_bi_ln)
  
  ### Cargo Loss: Gamma GLM
  cols_for_model  <- c("claim_amount", "cargo_type", "cargo_value", "weight",
                       "route_risk", "distance", "transit_duration", "vessel_age",
                       "debris_density", "pilot_experience", "container_type", "solar_radiation")
  cargo_sev_complete <- na.omit(cargo_sev[, cols_for_model])
  cargo_glm_full     <- glm(claim_amount ~ ., data = cargo_sev_complete, family = Gamma(link = "log"))
  print(summary(cargo_glm_full))
  
  cargo_glm <- step(cargo_glm_full, direction = "backward", trace = 0)
  print(summary(cargo_glm))
  # Significant covariates: cargo_type + cargo_value + weight + route_risk + debris_density + solar_radiation
  
  ### Equipment Failure: Gamma GLM
  # Note: equipment_age intentionally excluded to prevent massive row deletion
  cols_for_ef    <- c("claim_amount", "equipment_type", "solar_system", "maintenance_int", "usage_int")
  ef_sev_complete <- na.omit(ef_sev[, cols_for_ef])
  ef_glm_full     <- glm(claim_amount ~ ., data = ef_sev_complete, family = Gamma(link = "log"))
  
  ef_glm <- step(ef_glm_full, direction = "backward")
  print(summary(ef_glm))
  # Significant covariates: equipment_type + solar_system + usage_int
  
  ### Workers' Compensation: marginal lognormal (fitted on $thousands)
  fit_wc_ln <- fitdist(wc_sev$claim_amt_scaled, "lnorm")
  
  
  # =======
  # AGGREGATE LOSS MODELLING
  
  # Requires the objects fitted above. Solar system proxy mapping:
  #   Helionis Cluster --> Helionis Cluster (direct)
  #   Bayesia System   --> Epsilon          (binary star radiation environment)
  #   Oryn Delta       --> Zeta             (frontier dwarf-star system)
  #   Cargo has no solar_system covariate, so no proxy needed
  
  N_SIM <- 10000
  
  # ======
  # UTILITY FUNCTIONS
  
  # simulate_aggregate(): Collective risk model: S = sum_{i=1..N} X_i
  #
  # mu_pool   : severity mean parameter pool
  #               lnorm --> meanlog (length 1, no covariates for BI/WC)
  #               gamma --> E[X] from predict(glm, type="response") per row
  # lambda    : expected annual claim count (rate × exposure volume)
  # freq_dist : "poisson", "negbin"
  # theta     : NegBin size parameter (required for negbin)
  # prob_w    : sampling weights for mu_pool (NULL = uniform)
  #             EF uses claim-count weights per equipment type
  # sev_dist  : "lnorm", "gamma"
  # sigma     : sdlog (required for lnorm)
  # phi       : Gamma dispersion = summary(model)$dispersion (required for gamma)
  #             Parameterisation: X ~ Gamma(shape=1/phi, scale=phi*mu_i)
  #             so E[X_i]=mu_i and Var[X_i]=phi*mu_i^2
  # scale     : 1e6 for BI (fitted on $millions), 1000 for WC (fitted on $thousands),
  #             1 for EF and Cargo (fitted on raw $)
  
  simulate_aggregate <- function(mu_pool, lambda, N_SIM,
                                 freq_dist = "poisson", theta = NULL,
                                 prob_w    = NULL,
                                 sev_dist  = "lnorm",
                                 sigma     = NULL,
                                 phi       = NULL,
                                 scale     = 1) {
    if (freq_dist == "negbin" && is.null(theta))
      stop("theta is required for negbin simulation")
    if (sev_dist == "lnorm"  && is.null(sigma))
      stop("sigma (sdlog) is required for lognormal simulation")
    if (sev_dist == "gamma"  && (is.null(phi) || phi <= 0))
      stop("phi > 0 is required for gamma simulation")
    
    pool_n   <- length(mu_pool)
    out      <- numeric(N_SIM)
    marginal <- (pool_n == 1L)
    
    for (iter in seq_len(N_SIM)) {
      n <- if (freq_dist == "poisson")
        rpois(1L, lambda)
      else
        rnbinom(1L, mu = lambda, size = theta)
      
      if (n == 0L) { out[iter] <- 0.0; next }
      
      if (sev_dist == "lnorm") {
        ml <- if (marginal) mu_pool[1L] else
          mu_pool[sample.int(pool_n, n, replace = TRUE, prob = prob_w)]
        x  <- rlnorm(n, meanlog = ml, sdlog = sigma)
      } else {
        mu_i <- if (marginal) mu_pool[1L] else
          mu_pool[sample.int(pool_n, n, replace = TRUE, prob = prob_w)]
        x    <- rgamma(n, shape = 1.0 / phi, scale = phi * mu_i)
      }
      
      out[iter] <- sum(x) * scale
    }
    out
  }
  
  summarise_sim <- function(x, label) {
    q99 <- as.numeric(quantile(x, 0.99))
    data.frame(
      Line   = label,
      Mean   = mean(x),
      SD     = sd(x),
      CV     = round(sd(x) / mean(x), 3),
      P50    = as.numeric(quantile(x, 0.50)),
      P75    = as.numeric(quantile(x, 0.75)),
      P95    = as.numeric(quantile(x, 0.95)),
      P99    = q99,
      TVaR99 = mean(x[x >= q99])
    )
  }
  
  fmt_tbl <- function(df, title) {
    cat(sprintf("\n=== %s ===\n", title))
    df2     <- df
    big_col <- sapply(df2, is.numeric) & (names(df2) != "CV")
    df2[big_col] <- lapply(df2[big_col], function(x)
      formatC(round(x), format = "f", digits = 0, big.mark = ","))
    print(df2, row.names = FALSE)
  }
  
  # ======
  # VALIDATE REQUIRED OBJECTS
  required_objs <- c("bi_poisson", "cargo_nb", "ef_poisson", "wc_poisson",
                     "fit_bi_ln",  "cargo_glm", "ef_glm",    "fit_wc_ln",
                     "bi_clean",   "cargo_clean", "ef_clean", "wc_clean",
                     "cargo_sev_complete")
  
  missing <- required_objs[!vapply(required_objs, exists, logical(1))]
  if (length(missing) > 0)
    stop(paste("Missing required objects:", paste(missing, collapse = ", ")))
  
  cat("cargo_glm formula (after backward selection):\n"); print(formula(cargo_glm))
  cat("\nef_glm formula (after backward selection):\n");  print(formula(ef_glm))
  
  # =======
  # 3. PARAMETERS
  # BI/WC/Cargo parameters
  # Parameters that are manually adjusted are directly supported by the Online Encyclopedia (route_risk, debris_density, solar_radiation, supply_chain, gravity for Bayesia).
  #Other vars are derived from training-data medians/modes per proxy solar system.
  
  cq_systems <- c("Helionis Cluster", "Bayesia System", "Oryn Delta")
  
  sys_proxy <- c("Helionis Cluster" = "Helionis Cluster",
                 "Bayesia System"   = "Epsilon",
                 "Oryn Delta"       = "Zeta")
  
  # Exposure counts
  # n_workers:  36,000 employees allocated proportional to mine counts (30:15:10)
  #             Source: CQ encyclopedia
  # n_vessels:  sum of all vessel types per system from CQ inventory sheet
  # n_stations: directly stated in CQ encyclopedia (30 / 15 / 10 mines)
  n_workers  <- c("Helionis Cluster" = 19636L, "Bayesia System" = 9818L, "Oryn Delta" = 6545L)
  n_vessels  <- c("Helionis Cluster" = 1160L,  "Bayesia System" = 1128L, "Oryn Delta" =  774L)
  n_stations <- c("Helionis Cluster" = 30L,    "Bayesia System" = 15L,   "Oryn Delta" =   10L)
  
  # ======
  # STEP 1: Derive training-data medians/modes per proxy solar system
  # Assumption: CQ's operational profile in each system matches the historical distribution of the corresponding proxy system in the training data.
  
  proxy_systems <- c("Helionis Cluster", "Epsilon", "Zeta")
  
  # BI parameters by proxy system
  bi_data_params <- do.call(rbind, lapply(proxy_systems, function(ps) {
    rows <- bi_clean[bi_clean$solar_system == ps & !is.na(bi_clean$solar_system), ]
    data.frame(
      proxy_sys      = ps,
      prod_load_med  = median(rows$production_load,     na.rm = TRUE),
      maint_freq_med = median(rows$maintenance_freq,    na.rm = TRUE),
      energy_bk_mode = as.numeric(names(sort(table(rows$energy_backup_score), dec = TRUE))[1]),
      stringsAsFactors = FALSE
    )
  }))
  
  # WC psych_stress median by proxy system
  wc_psi_params <- do.call(rbind, lapply(proxy_systems, function(ps) {
    rows <- wc_clean[wc_clean$solar_system == ps & !is.na(wc_clean$solar_system), ]
    data.frame(
      proxy_sys        = ps,
      psych_stress_med = median(rows$psych_stress_index, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  
  # Cargo: pilot_exp and vessel_age: no solar_system covariate in cargo model.
  # Use training-data portfolio-wide medians as baseline.
  # vessel_age is then scaled by relative system maturity, proxied by the
  # inventory-derived weighted mean equipment age per system:
  #   Helionis WM age = 11.30 yrs  -->  scale = 1.00  (matches training median)
  #   Bayesia  WM age =  8.15 yrs  -->  scale = 0.73
  #   Oryn     WM age =  4.13 yrs  -->  scale = 0.36
  cargo_pilot_med  <- median(cargo_clean$pilot_experience, na.rm = TRUE)  # = 15.0
  cargo_vessel_med <- median(cargo_clean$vessel_age,       na.rm = TRUE)  # = 20.0
  
  vessel_age_scale <- c("Helionis Cluster" = 1.00,
                        "Bayesia System"   = 0.73,
                        "Oryn Delta"       = 0.36)
  
  pilot_exp_vals  <- setNames(rep(cargo_pilot_med, 3), cq_systems)
  vessel_age_vals <- setNames(round(cargo_vessel_med * vessel_age_scale[cq_systems], 1), cq_systems)
  
  cat("=== Data-derived BI parameters per proxy system ===\n"); print(bi_data_params)
  cat("\n=== Data-derived WC psych_stress per proxy system ===\n"); print(wc_psi_params)
  cat(sprintf("\nCargo pilot_exp (training median): %.1f [all systems]\n", cargo_pilot_med))
  cat(sprintf("Cargo vessel_age (scaled): H=%.1f  B=%.1f  O=%.1f\n",
              vessel_age_vals["Helionis Cluster"],
              vessel_age_vals["Bayesia System"],
              vessel_age_vals["Oryn Delta"]))
  
  # ======
  # STEP 2: Helper functions: look up training-data params by CQ system --> proxy
  get_bi_param <- function(param_col) {
    setNames(sapply(cq_systems, function(cs) {
      ps <- sys_proxy[cs]
      bi_data_params[bi_data_params$proxy_sys == ps, param_col]
    }), cq_systems)
  }
  
  get_wc_param <- function(param_col) {
    setNames(sapply(cq_systems, function(cs) {
      ps <- sys_proxy[cs]
      wc_psi_params[wc_psi_params$proxy_sys == ps, param_col]
    }), cq_systems)
  }
  
  # ======
  # STEP 3: Build sys_params
  #
  # Parameters retained from Online Encyclopedia (justified below):
  #   route_risk     : Helionis=3 (erratic debris, micro-collisions),
  #                     Bayesia=2  (heavily mapped, stable orbits),
  #                     Oryn=4     (asymmetric ring, orbital shear)
  #   debris_density : consistent with route_risk ordering
  #   solar_radiation: Helionis G2V low flare=0.50; Bayesia binary EM spikes=0.70;
  #                     Oryn M3V dim/mild=0.30
  #   supply_chain   : Bayesia "well-established orbital stations"=0.75;
  #                     Oryn "infrastructure rapidly evolving"=0.60;
  #                     training medians all ~0.495, so CQ above-median is reasonable
  #   gravity (Bayesia=1.300): directly stated "high-gravity" in encyclopedia
  #   gravity (Helionis=1.125, Oryn=1.113): assumed; flagged below
  # ======
  sys_params <- data.frame(
    cq_system = cq_systems,
    proxy_sys = unname(sys_proxy),
    
    # WC: gravity
    # Bayesia=1.300 encyclopedia-stated ("high-gravity with thin magnetosphere")
    # Helionis=1.125. Oryn=1.113:assumed (no explicit encyclopedia figure); flagged
    gravity      = c(1.125, 1.300, 1.113),
    
    # WC: psych_stress: training-data median per proxy system (replaces assumed 3,4,4)
    psych_stress = as.numeric(get_wc_param("psych_stress_med")),
    
    # Cargo: route_risk: encyclopedia-supported (see note above)
    route_risk_n = c(3L, 2L, 4L),   # numeric for cargo_nb
    route_risk_c = c("3","2","4"),   # character (--> factor) for cargo_glm
    
    # Cargo: debris_density and solar_radiation: encyclopedia-supported
    debris_density  = c(0.60, 0.30, 0.70),
    solar_radiation = c(0.50, 0.70, 0.30),
    
    # Cargo: pilot_exp: training-data median; no system-level breakdown available
    pilot_exp     = as.numeric(pilot_exp_vals[cq_systems]),
    
    # Cargo: vessel_age: training-data median scaled by inventory-derived system maturity
    vessel_age_cq = as.numeric(vessel_age_vals[cq_systems]),
    
    # BI: supply_chain: encyclopedia-supported (see note above)
    supply_chain  = c(0.80, 0.75, 0.60),
    
    # BI: maint_freq: training-data median per proxy system (replaces assumed 4,3,2)
    maint_freq_bi = as.numeric(get_bi_param("maint_freq_med")),
    
    # BI: prod_load: training-data median per proxy system (replaces assumed 0.80,0.75,0.85)
    prod_load     = as.numeric(get_bi_param("prod_load_med")),
    
    # BI: energy_backup: training-data mode per proxy system (replaces assumed all=4)
    #     Helionis=5, Epsilon(Bayesia)=2, Zeta(Oryn)=2
    energy_backup = as.numeric(get_bi_param("energy_bk_mode")),
    
    stringsAsFactors = FALSE
  )
  rownames(sys_params) <- cq_systems
  
  cat("\n=== Final sys_params ===\n"); print(sys_params)
  
  # ======
  # STEP 4: Equipment type map: explicit named mapping, replaces positional c()
  #
  # Training-data equipment_type levels (exact strings from ef_clean):
  #   "Fexstram Carrier", "Flux Rider", "Graviton Extractor",
  #   "Ion Pulverizer", "Quantum Bore", "ReglAggregators"
  # These match the CQ inventory names exactly.
  #
  # n_units       : CQ inventory sheet (unit counts per system)
  # usage_int     : (% in operation from inventory) × 24 hrs/day
  # maintenance_int: maintenance schedule in hours from inventory
  # equipment_age : inventory-derived weighted mean service years per system
  #                 H=11.30, B=8.15, O=4.13  (midpoints: <5 --> 2.5, 5-9 --> 7, 10-14 --> 12, etc.)
  # ======
  ef_type_map <- data.frame(
    label = c("Fexstram Carrier", "Flux Rider", "Graviton Extractor",
              "Ion Pulverizer",   "Quantum Bore", "ReglAggregators"),
    # Unit counts: Helionis / Bayesia / Oryn
    n_H = c(150, 1500, 240,  90, 300, 300),
    n_B = c( 75,  750, 120,  45, 150, 150),
    n_O = c( 50,  500,  80,  30, 100, 100),
    # % in operation × 24 hrs
    u_H = c(0.90, 0.80, 0.95, 0.50, 0.95, 0.80),
    u_B = c(0.75, 0.80, 0.80, 0.60, 0.80, 0.75),
    u_O = c(0.70, 0.75, 0.75, 0.50, 0.75, 0.70),
    # Maintenance schedule (hrs)
    m_H = c( 375, 1500,  750, 1000,  750, 1500),
    m_B = c( 400, 1000,  600,  750,  600, 1000),
    m_O = c( 250,  300,  500,  500,  500,  300),
    stringsAsFactors = FALSE
  )
  
  # Runtime validation: confirm all map labels exist in the fitted model
  ef_model_levels <- levels(factor(ef_clean$equipment_type))
  ef_model_levels <- ef_model_levels[ef_model_levels != "NA"]
  missing_types   <- setdiff(ef_type_map$label, ef_model_levels)
  extra_types     <- setdiff(ef_model_levels, ef_type_map$label)
  if (length(missing_types) > 0)
    warning("Equipment types in map NOT in model: ", paste(missing_types, collapse = ", "))
  if (length(extra_types)   > 0)
    warning("Equipment types in model NOT in map (ignored): ", paste(extra_types, collapse = ", "))
  cat("\nEquipment type validation:\n")
  cat("  Model levels:", paste(sort(ef_model_levels), collapse = ", "), "\n")
  cat("  Map labels:  ", paste(sort(ef_type_map$label), collapse = ", "), "\n")
  cat("  All matched: ", length(missing_types) == 0 && length(extra_types) == 0, "\n")
  
  
  # =======
  # 4. BUSINESS INTERRUPTION
  
  # 4A. BI FREQUENCY
  # Covariates: solar_system, production_load, energy_backup_score,
  #   supply_chain_index, avg_crew_exp, maintenance_freq, safety_compliance
  # One prediction row per system; rate/station × n_stations = lambda
  
  bi_crew_med   <- median(bi_clean$avg_crew_exp,      na.rm = TRUE)
  bi_safety_med <- median(bi_clean$safety_compliance, na.rm = TRUE)
  
  bi_lambda <- setNames(numeric(3), cq_systems)
  
  cat("BI Frequency (Poisson):\n")
  for (sys in cq_systems) {
    p  <- sys_params[sys, ]
    nd <- data.frame(
      solar_system        = p$proxy_sys,
      production_load     = p$prod_load,
      energy_backup_score = p$energy_backup,
      supply_chain_index  = p$supply_chain,
      avg_crew_exp        = bi_crew_med,
      maintenance_freq    = p$maint_freq_bi,
      safety_compliance   = bi_safety_med,
      exposure            = 1.0
    )
    rate           <- as.numeric(predict(bi_poisson, newdata = nd, type = "raw"))
    bi_lambda[sys] <- rate * n_stations[sys]
    cat(sprintf("  %-20s  rate/station=%.5f  n=%d  lambda=%.3f\n",
                sys, rate, n_stations[sys], bi_lambda[sys]))
  }
  
  # 4B. BI SEVERITY: marginal lognormal (no covariates; fitted on $millions)
  bi_meanlog <- fit_bi_ln$estimate["meanlog"]
  bi_sdlog   <- fit_bi_ln$estimate["sdlog"]
  bi_mean_x  <- exp(bi_meanlog + bi_sdlog^2 / 2) * 1e6
  
  cat(sprintf("\nBI Severity (LogNormal): meanlog=%.4f  sdlog=%.4f  E[claim]=Đ%.0f\n",
              bi_meanlog, bi_sdlog, bi_mean_x))
  
  # 4C. BI SIMULATION
  cat("Simulating BI...\n")
  bi_sim <- setNames(vector("list", 3), cq_systems)
  for (sys in cq_systems) {
    bi_sim[[sys]] <- simulate_aggregate(
      mu_pool = bi_meanlog, lambda = bi_lambda[sys], N_SIM = N_SIM,
      freq_dist = "poisson", sev_dist = "lnorm", sigma = bi_sdlog, scale = 1e6
    )
  }
  
  bi_sim_all <- Reduce("+", bi_sim)
  bi_det     <- sum(bi_lambda * bi_mean_x)
  cat(sprintf("  BI sanity: det=Đ%.0f  sim_mean=Đ%.0f  ratio=%.3f\n",
              bi_det, mean(bi_sim_all), mean(bi_sim_all) / bi_det))
  
  fmt_tbl(rbind(
    summarise_sim(bi_sim[["Helionis Cluster"]], "BI - Helionis"),
    summarise_sim(bi_sim[["Bayesia System"]],   "BI - Bayesia"),
    summarise_sim(bi_sim[["Oryn Delta"]],       "BI - Oryn"),
    summarise_sim(bi_sim_all,                   "BI - All Systems")
  ), "BI Aggregate Loss (Đ)")
  
  
  # =======
  # 5. WORKERS' COMPENSATION
  
  # 5A. WC FREQUENCY
  # Strategy: build one row per occupation (11 rows) with all other covariates at representative values (median/mode from training data); override gravity and psych_stress per system. Predict rate/worker, then take exposure-weighted average across occupations. 11 rows × 3 systems = 33 calls.
  
  wc_occ_lvls <- sort(unique(wc_clean$occupation))
  
  occ_exp_wts <- tapply(wc_clean$exposure, wc_clean$occupation, sum)
  occ_exp_wts <- occ_exp_wts[wc_occ_lvls]
  occ_wts     <- as.numeric(occ_exp_wts / sum(occ_exp_wts))
  
  wc_emp_mode <- names(sort(table(wc_clean$employment_type),      dec = TRUE))[1]
  wc_acc_mode <- as.numeric(names(sort(table(wc_clean$accident_history_flag), dec = TRUE))[1])
  wc_hrs_mode <- as.numeric(names(sort(table(wc_clean$hours_per_week),        dec = TRUE))[1])
  wc_exp_med  <- median(wc_clean$experience_yrs,    na.rm = TRUE)
  wc_sup_med  <- median(wc_clean$supervision_level, na.rm = TRUE)
  wc_sal_med  <- median(wc_clean$base_salary,       na.rm = TRUE)
  wc_safety_bas <- 5  # mode = 5 (50% of records)
  wc_gear_bas   <- 5  # mode = 5 (44% of records)
  
  cat("WC Frequency (Poisson):\n")
  cat(sprintf("  Occupation weights (exposure-wtd): %s\n",
              paste(sprintf("%s=%.1f%%", wc_occ_lvls, occ_wts * 100), collapse = ", ")))
  
  wc_lambda <- setNames(numeric(3), cq_systems)
  
  for (sys in cq_systems) {
    p <- sys_params[sys, ]
    nd <- data.frame(
      solar_system            = p$proxy_sys,
      occupation              = wc_occ_lvls,
      employment_type         = wc_emp_mode,
      experience_yrs          = wc_exp_med,
      accident_history_flag   = wc_acc_mode,
      psych_stress_index      = p$psych_stress,
      hours_per_week          = wc_hrs_mode,
      supervision_level       = wc_sup_med,
      gravity_level           = p$gravity,
      safety_training_index   = wc_safety_bas,
      protective_gear_quality = wc_gear_bas,
      base_salary             = wc_sal_med,
      exposure                = 1.0
    )
    preds <- as.numeric(predict(wc_poisson, newdata = nd, type = "raw"))
    if (any(is.na(preds))) {
      warning(sprintf("WC %s: %d NA predictions", sys, sum(is.na(preds))))
      preds[is.na(preds)] <- 0
    }
    rate           <- sum(occ_wts * preds)
    wc_lambda[sys] <- rate * n_workers[sys]
    cat(sprintf("  %-20s  rate/worker=%.6f  n=%d  lambda=%.1f\n",
                sys, rate, n_workers[sys], wc_lambda[sys]))
  }
  
  # 5B. WC SEVERITY: marginal lognormal (no covariates; fitted on $thousands)
  wc_meanlog <- fit_wc_ln$estimate["meanlog"]
  wc_sdlog   <- fit_wc_ln$estimate["sdlog"]
  wc_mean_x  <- exp(wc_meanlog + wc_sdlog^2 / 2) * 1000
  
  cat(sprintf("\nWC Severity (LogNormal): meanlog=%.4f  sdlog=%.4f  E[claim]=Đ%.0f\n",
              wc_meanlog, wc_sdlog, wc_mean_x))
  
  # 5C. WC SIMULATION
  cat("Simulating WC...\n")
  wc_sim <- setNames(vector("list", 3), cq_systems)
  for (sys in cq_systems) {
    wc_sim[[sys]] <- simulate_aggregate(
      mu_pool = wc_meanlog, lambda = wc_lambda[sys], N_SIM = N_SIM,
      freq_dist = "poisson", sev_dist = "lnorm", sigma = wc_sdlog, scale = 1000
    )
  }
  
  wc_sim_all <- Reduce("+", wc_sim)
  wc_det     <- sum(wc_lambda * wc_mean_x)
  cat(sprintf("  WC sanity: det=Đ%.0f  sim_mean=Đ%.0f  ratio=%.3f\n",
              wc_det, mean(wc_sim_all), mean(wc_sim_all) / wc_det))
  
  fmt_tbl(rbind(
    summarise_sim(wc_sim[["Helionis Cluster"]], "WC - Helionis"),
    summarise_sim(wc_sim[["Bayesia System"]],   "WC - Bayesia"),
    summarise_sim(wc_sim[["Oryn Delta"]],       "WC - Oryn"),
    summarise_sim(wc_sim_all,                   "WC - All Systems")
  ), "WC Aggregate Loss (Đ)")
  
  
  # =======
  # 6. EQUIPMENT FAILURE
  
  # 6A. EF INVENTORY
  # ef_inventory built from ef_type_map (explicit named mapping defined in
  # Step 4 above). predict.glm() silently ignores columns not in ef_glm formula,
  # so including all inventory columns is safe and future-proof if step() changes.
  
  ef_solar_lvls <- sort(unique(na.omit(as.character(ef_clean$solar_system))))
  
  ef_inventory <- do.call(rbind, lapply(cq_systems, function(sys) {
    ps  <- sys_proxy[sys]
    sfx <- switch(sys,
                  "Helionis Cluster" = "_H",
                  "Bayesia System"   = "_B",
                  "Oryn Delta"       = "_O")
    data.frame(
      equipment_type  = factor(ef_type_map$label, levels = ef_model_levels),
      cq_system       = sys,
      solar_system    = factor(ps, levels = ef_solar_lvls),
      n_units         = ef_type_map[[paste0("n", sfx)]],
      usage_int       = ef_type_map[[paste0("u", sfx)]] * 24,   # % --> hrs/day
      maintenance_int = ef_type_map[[paste0("m", sfx)]],
      equipment_age   = switch(sys,
                               "Helionis Cluster" = 11.30,
                               "Bayesia System"   =  8.15,
                               "Oryn Delta"       =  4.13),
      exposure        = 1.0,
      stringsAsFactors = FALSE
    )
  }))
  
  # 6B. EF FREQUENCY
  ef_inventory$rate_per_unit <- as.numeric(
    predict(ef_poisson, newdata = ef_inventory, type = "raw"))
  ef_inventory$row_exp_n <- ef_inventory$rate_per_unit * ef_inventory$n_units
  
  # 6C. EF SEVERITY: Gamma GLM
  phi_ef <- summary(ef_glm)$dispersion
  cat(sprintf("EF Gamma phi=%.6f\n", phi_ef))
  
  ef_inventory$mu_sev       <- as.numeric(predict(ef_glm, newdata = ef_inventory, type = "response"))
  ef_inventory$row_exp_loss <- ef_inventory$row_exp_n * ef_inventory$mu_sev
  
  cat("EF inventory (by system × equipment type):\n")
  print(
    ef_inventory %>%
      dplyr::select(cq_system, equipment_type, n_units,
                    rate_per_unit, row_exp_n, mu_sev, row_exp_loss) %>%
      mutate(across(c(rate_per_unit, row_exp_n, mu_sev, row_exp_loss), ~round(., 2))),
    row.names = FALSE
  )
  
  # 6D. EF SIMULATION
  cat("\nSimulating EF...\n")
  ef_sim <- setNames(vector("list", 3), cq_systems)
  for (sys in cq_systems) {
    rows    <- ef_inventory[ef_inventory$cq_system == sys, ]
    lam_sys <- sum(rows$row_exp_n)
    prob_w  <- rows$row_exp_n / lam_sys
    cat(sprintf("  %-20s  lambda=%.1f  E[loss]=Đ%.0f\n", sys, lam_sys, sum(rows$row_exp_loss)))
    ef_sim[[sys]] <- simulate_aggregate(
      mu_pool = rows$mu_sev, lambda = lam_sys, N_SIM = N_SIM,
      freq_dist = "poisson", prob_w = prob_w, sev_dist = "gamma", phi = phi_ef
    )
  }
  
  ef_sim_all <- Reduce("+", ef_sim)
  ef_det     <- sum(ef_inventory$row_exp_loss)
  cat(sprintf("  EF sanity: det=Đ%.0f  sim_mean=Đ%.0f  ratio=%.3f\n",
              ef_det, mean(ef_sim_all), mean(ef_sim_all) / ef_det))
  
  fmt_tbl(rbind(
    summarise_sim(ef_sim[["Helionis Cluster"]], "EF - Helionis"),
    summarise_sim(ef_sim[["Bayesia System"]],   "EF - Bayesia"),
    summarise_sim(ef_sim[["Oryn Delta"]],       "EF - Oryn"),
    summarise_sim(ef_sim_all,                   "EF - All Systems")
  ), "EF Aggregate Loss (Đ)")
  
  
  # =======
  # 7. CARGO LOSS
  
  phi_cargo     <- summary(cargo_glm)$dispersion
  cargo_rr_lvls <- levels(cargo_sev_complete$route_risk)
  cat(sprintf("Cargo bootstrap pool: %d rows\n", nrow(cargo_sev_complete)))
  cat(sprintf("Cargo Gamma phi=%.6f  |  route_risk levels: %s\n",
              phi_cargo, paste(cargo_rr_lvls, collapse = ", ")))
  
  # 7A. CARGO SEVERITY: BOOTSTRAP MU PRECOMPUTATION
  # Override system-level parameters (route_risk, debris_density, solar_radiation)
  # in the bootstrap pool; cargo_type, cargo_value, and weight vary naturally,
  # preserving the cargo_value distribution and its correlation with weight.
  
  #Precomputing cargo bootstrap mu pools
  cargo_mu_pool <- setNames(vector("list", 3), cq_systems)
  for (sys in cq_systems) {
    p      <- sys_params[sys, ]
    pool   <- cargo_sev_complete
    rr_str <- as.character(p$route_risk_c)
    if (!rr_str %in% cargo_rr_lvls)
      stop(sprintf("route_risk='%s' not valid for %s", rr_str, sys))
    pool$route_risk      <- factor(rr_str, levels = cargo_rr_lvls)
    pool$debris_density  <- p$debris_density
    pool$solar_radiation <- p$solar_radiation
    mu_vec <- predict(cargo_glm, newdata = pool, type = "response")
    mu_vec <- mu_vec[!is.na(mu_vec)]
    cargo_mu_pool[[sys]] <- mu_vec
    cat(sprintf("  %-20s  pool_n=%d  mean_sev=Đ%.0f\n", sys, length(mu_vec), mean(mu_vec)))
  }
  
  # 7B. CARGO FREQUENCY
  ct_lvls  <- sort(unique(cargo_clean$cargo_type))
  ct_wts   <- as.numeric(prop.table(table(cargo_clean$cargo_type)[ct_lvls]))
  med_cont <- names(sort(table(cargo_clean$container_type), dec = TRUE))[1]
  med_wt   <- median(cargo_clean$weight,           na.rm = TRUE)
  med_dist <- median(cargo_clean$distance,         na.rm = TRUE)
  med_td   <- median(cargo_clean$transit_duration, na.rm = TRUE)
  
  cat("\nCargo Frequency (NegBin):\n")
  cargo_lambda <- setNames(numeric(3), cq_systems)
  for (sys in cq_systems) {
    p <- sys_params[sys, ]
    rates_by_type <- vapply(ct_lvls, function(ct) {
      nd <- data.frame(
        cargo_type       = ct,
        weight           = med_wt,
        route_risk       = as.numeric(p$route_risk_n),
        distance         = med_dist,
        transit_duration = med_td,
        pilot_experience = p$pilot_exp,
        vessel_age       = p$vessel_age_cq,
        container_type   = med_cont,
        solar_radiation  = p$solar_radiation,
        debris_density   = p$debris_density,
        exposure         = 1.0
      )
      as.numeric(predict(cargo_nb, newdata = nd, type = "raw"))
    }, numeric(1))
    rate_per_vessel    <- sum(ct_wts * rates_by_type)
    cargo_lambda[sys]  <- rate_per_vessel * n_vessels[sys]
    cat(sprintf("  %-20s  rate/vessel=%.5f  n=%d  lambda=%.1f\n",
                sys, rate_per_vessel, n_vessels[sys], cargo_lambda[sys]))
  }
  
  cargo_theta <- cargo_nb$finalModel$theta
  cat(sprintf("  NegBin theta=%.4f\n", cargo_theta))
  
  # 7C. CARGO SIMULATION
  cargo_sim <- setNames(vector("list", 3), cq_systems)
  for (sys in cq_systems) {
    mu_pool_sys <- cargo_mu_pool[[sys]]
    cat(sprintf("  %-20s  lambda=%.1f  E[loss]=Đ%.0f\n",
                sys, cargo_lambda[sys], cargo_lambda[sys] * mean(mu_pool_sys)))
    cargo_sim[[sys]] <- simulate_aggregate(
      mu_pool = mu_pool_sys, lambda = cargo_lambda[sys], N_SIM = N_SIM,
      freq_dist = "negbin", theta = cargo_theta, sev_dist = "gamma", phi = phi_cargo
    )
  }
  
  cargo_sim_all <- Reduce("+", cargo_sim)
  cargo_det     <- sum(vapply(cq_systems, function(s)
    cargo_lambda[s] * mean(cargo_mu_pool[[s]]), numeric(1)))
  cat(sprintf("  Cargo sanity: det=Đ%.0f  sim_mean=Đ%.0f  ratio=%.3f\n",
              cargo_det, mean(cargo_sim_all), mean(cargo_sim_all) / cargo_det))
  
  fmt_tbl(rbind(
    summarise_sim(cargo_sim[["Helionis Cluster"]], "Cargo - Helionis"),
    summarise_sim(cargo_sim[["Bayesia System"]],   "Cargo - Bayesia"),
    summarise_sim(cargo_sim[["Oryn Delta"]],       "Cargo - Oryn"),
    summarise_sim(cargo_sim_all,                   "Cargo - All Systems")
  ), "Cargo Aggregate Loss (Đ)")
  
  
  # =======
  # 8. PORTFOLIO SUMMARY
  
  full_portfolio <- bi_sim_all + wc_sim_all + ef_sim_all + cargo_sim_all
  
  exp_tbl <- data.frame(
    System = cq_systems,
    BI     = bi_lambda * bi_mean_x,
    WC     = wc_lambda * wc_mean_x,
    EF     = vapply(cq_systems, function(s)
      sum(ef_inventory$row_exp_loss[ef_inventory$cq_system == s]), numeric(1)),
    Cargo  = vapply(cq_systems, function(s)
      cargo_lambda[s] * mean(cargo_mu_pool[[s]]), numeric(1))
  ) %>% mutate(Total = BI + WC + EF + Cargo)
  
  cat("Expected Annual Losses by Line and System (Đ):\n")
  print(
    exp_tbl %>% mutate(across(where(is.numeric),
                              ~formatC(round(.), format = "f", digits = 0, big.mark = ","))),
    row.names = FALSE
  )
  cat(sprintf("Total expected annual loss: Đ%s\n",
              formatC(round(sum(exp_tbl$Total)), format = "f", digits = 0, big.mark = ",")))
  
  fmt_tbl(summarise_sim(full_portfolio, "FULL PORTFOLIO"),
          "Full Portfolio Aggregate Loss Distribution (Đ)")
  
  stress <- rbind(
    summarise_sim(bi_sim_all,     "BI"),
    summarise_sim(wc_sim_all,     "WC"),
    summarise_sim(ef_sim_all,     "EF"),
    summarise_sim(cargo_sim_all,  "Cargo"),
    summarise_sim(full_portfolio, "Portfolio")
  ) %>%
    dplyr::select(Line, Mean, P99, TVaR99) %>%
    mutate(
      P99_to_Mean  = round(P99    / Mean, 2),
      TVaR_to_Mean = round(TVaR99 / Mean, 2),
      Mean   = formatC(round(Mean),   format = "f", digits = 0, big.mark = ","),
      P99    = formatC(round(P99),    format = "f", digits = 0, big.mark = ","),
      TVaR99 = formatC(round(TVaR99), format = "f", digits = 0, big.mark = ",")
    )
  
  cat("\n=== STRESS TEST: 1-IN-100 YEAR (P99 and TVaR99) ===\n")
  print(stress, row.names = FALSE)
  
  
  # =======
  # 9. SHORT-TERM AND LONG-TERM PROJECTIONS
  #    Addresses Project Objective 1: short and long-term ranges for aggregate
  #      (a) costs, (b) returns (investment income), (c) net revenue
  #    Addresses Project Objective 2: expected values, variances, tail behaviours
  #
  # DATA SOURCES:
  #   Interest & Inflation: srcsc2026interestandinflation.xlsx (2160–2174)
  #   Exposure growth:      In the encyclopedia it mentions 25% over 10yr for Helionis/Bayesia;
  #                          15% for Oryn Delta
  #
  # FRAMEWORK:
  #   Short-term  = Year 1 (already simulated in Sections 4–8)
  #   Long-term   = 10-year cumulative present value
  #
  #   Each year t = 1..10 is resimulated with updated parameters:
  #     lambda_t  = lambda_1 × (1 + CAGR)^(t-1)         [exposure growth]
  #     mu_t      = mu_1 × (1 + inflation)^(t-1)          [severity growth]
  #     lnorm: shift meanlog by log(inflation_factor)
  #     gamma:  multiply mu_pool by inflation_factor
  #
  # LOADING ASSUMPTIONS (documented):
  #   expense_ratio  = 0.30  (30% of premium: admin, reinsurance, commissions:
  #                            typical for commercial space P&C)
  #   profit_margin  = 0.05  (5% target profit margin: conservative for new entrant)
  #   => load_factor = 1 / (1 - 0.30 - 0.05) = 1.538
  #
  # INVESTMENT INCOME:
  #   Reserve float about 0.5 × E[annual loss] (mid-year approximation: premiums
  #   collected upfront, claims paid throughout the year, so about half the year's
  #   expected loss is held on average). Invested at the short-term risk-free rate.
  # =======
  
  library(readxl)
  
  # ======
  # 9A. LOAD INTEREST & INFLATION DATA
  int_inf_path <- "srcsc-2026-interest-and-inflation.xlsx"  
  
  # skip=2: skip title row AND header row; assign names manually by position
  # This avoids messy long column names from the original headers
  int_inf_raw <- read_excel(int_inf_path, sheet = "Sheet1",
                            col_names = FALSE, skip = 2)
  int_inf <- int_inf_raw[, 1:5]
  colnames(int_inf) <- c("year", "inflation", "overnight", "r1yr", "r10yr")
  int_inf <- int_inf %>%
    mutate(across(everything(), as.numeric)) %>%
    filter(!is.na(year))                      # drop any blank rows
  
  # Forward projection: trailing 5-year averages (2170–2174)
  recent5   <- int_inf %>% filter(year >= 2170)
  inf_fwd   <- mean(recent5$inflation, na.rm = TRUE)  # 4.23%: claims inflation
  r1yr_fwd  <- mean(recent5$r1yr,     na.rm = TRUE)  # 3.71%: investment rate
  r10yr_fwd <- mean(recent5$r10yr,    na.rm = TRUE)  # 3.76%: discount rate
  
  #ECONOMIC ASSUMPTIONS (5-yr trailing avg 2170–2174)
  cat(sprintf("  Claims inflation rate:    %.2f%%\n", inf_fwd   * 100))
  cat(sprintf("  Short-term invest rate:  %.2f%%\n",  r1yr_fwd  * 100))
  cat(sprintf("  Long-term discount rate: %.2f%%\n",  r10yr_fwd * 100))
  
  # ======
  # 9B. PRICING AND EXPENSE PARAMETERS
  expense_ratio <- 0.30
  profit_margin <- 0.05
  load_factor   <- 1 / (1 - expense_ratio - profit_margin)  # = 1.538
  
  #PRICING PARAMETERS
  cat(sprintf("  Expense ratio:   %.0f%%\n",    expense_ratio * 100))
  cat(sprintf("  Profit margin:   %.0f%%\n",    profit_margin * 100))
  cat(sprintf("  Load factor:     %.4f\n",      load_factor))
  
  # ======
  # 9C. EXPOSURE GROWTH RATES (CQ encyclopedia)
  #   Helionis + Bayesia: 25% over 10 years  -->  CAGR = 1.25^(1/10) − 1 = 2.26%
  #   Oryn Delta:         15% over 10 years  -->  CAGR = 1.15^(1/10) − 1 = 1.41%
  growth_rate <- c(
    "Helionis Cluster" = 1.25^(1/10) - 1,
    "Bayesia System"   = 1.25^(1/10) - 1,
    "Oryn Delta"       = 1.15^(1/10) - 1
  )
  
  # Portfolio-weighted CAGR (weighted by year-1 expected loss contribution)
  el_by_sys <- setNames(exp_tbl$Total, cq_systems)
  wtd_growth <- sum(el_by_sys * growth_rate[cq_systems]) / sum(el_by_sys)
  
  #EXPOSURE GROWTH (from  encyclopedia)
  for (sys in cq_systems)
    cat(sprintf("  %-20s  CAGR = %.2f%%\n", sys, growth_rate[sys] * 100))
  cat(sprintf("  Portfolio-weighted CAGR: %.2f%%\n", wtd_growth * 100))
  
  # ======
  # funciton summarise_sim_rev()
  # For NET REVENUE, risk lies in the LEFT tail (bad = low/negative net revenue).
  # summarise_sim() reports the upper tail (P99 = best outcome), which is
  # misleading for revenue. This function instead reports lower-tail risk:
  #   P05, P01  = 1-in-20 and 1-in-100 year bad revenue outcomes
  #   TVaR01    = average net revenue conditional on being in worst 1%
  # ======
  summarise_sim_rev <- function(x, label) {
    q01 <- as.numeric(quantile(x, 0.01))
    data.frame(
      Line       = label,
      Mean       = mean(x),
      SD         = sd(x),
      Variance   = var(x),
      P50        = as.numeric(quantile(x, 0.50)),
      P25        = as.numeric(quantile(x, 0.25)),
      P05        = as.numeric(quantile(x, 0.05)),   # 1-in-20 bad year
      P01        = q01,                              # 1-in-100 bad year
      TVaR01     = mean(x[x <= q01]),               # avg of worst 1%
      Pr_Loss    = mean(x < 0)                      # probability of underwriting loss
    )
  }
  
  # Add Variance to cost summarise for completeness
  summarise_sim_cost <- function(x, label) {
    q99 <- as.numeric(quantile(x, 0.99))
    data.frame(
      Line     = label,
      Mean     = mean(x),
      SD       = sd(x),
      Variance = var(x),
      CV       = round(sd(x) / mean(x), 3),
      P50      = as.numeric(quantile(x, 0.50)),
      P75      = as.numeric(quantile(x, 0.75)),
      P95      = as.numeric(quantile(x, 0.95)),
      P99      = q99,
      TVaR99   = mean(x[x >= q99])
    )
  }
  
  # ======
  # precomputing the cargo mu pool weights 
  # Cargo severity pool sampled proportional to each system's lambda,
  # not uniformly. Uniform sampling would underweight Helionis (37.9% of claims)
  # and overweight Oryn (25.3%), since all pools have equal row count.
  # Hence, build a single weighted pool by resampling each system's mu_pool in
  # proportion to its lambda share.
  # ======
  cargo_lam_total <- sum(cargo_lambda)
  cargo_lam_wts   <- cargo_lambda / cargo_lam_total   # H=37.9%, B=36.8%, O=25.3%
  
  cargo_pool_target <- length(cargo_mu_pool[[cq_systems[1]]])
  
  set.seed(5100)
  cargo_mu_combined <- unlist(lapply(cq_systems, function(sys) {
    n_draw <- round(cargo_pool_target * cargo_lam_wts[sys])
    pool   <- cargo_mu_pool[[sys]]
    sample(pool, size = n_draw, replace = TRUE)
  }))
  
  cat(sprintf("\nCargo combined mu pool: %d rows (lambda-weighted: H=%.1f%% B=%.1f%% O=%.1f%%)\n",
              length(cargo_mu_combined),
              cargo_lam_wts["Helionis Cluster"] * 100,
              cargo_lam_wts["Bayesia System"]   * 100,
              cargo_lam_wts["Oryn Delta"]        * 100))
  
  # ======
  # 9F. 10-YEAR PROJECTION LOOP
  N_YEARS <- 10
  
  proj_costs   <- matrix(0, N_SIM, N_YEARS)   # simulated claim cost per year
  proj_premium <- numeric(N_YEARS)             # deterministic gross premium per year
  proj_net_rev <- matrix(0, N_SIM, N_YEARS)   # simulated net revenue per year
  proj_pv_cost <- numeric(N_SIM)              # PV of cumulative costs
  proj_pv_rev  <- numeric(N_SIM)              # PV of cumulative net revenue
  
  yr_exp_cost    <- numeric(N_YEARS)
  yr_exp_premium <- numeric(N_YEARS)
  yr_exp_net_rev <- numeric(N_YEARS)
  
  # Portfolio-level year-1 inputs
  total_lam_y1 <- sum(bi_lambda) + sum(wc_lambda) +
    sum(ef_inventory$row_exp_n) + sum(cargo_lambda)
  total_el_y1  <- sum(exp_tbl$Total)
  
  # EF: combined mu pool weighted by row_exp_n (done inside simulate_aggregate via prob_w)
  ef_prob_w_base <- ef_inventory$row_exp_n / sum(ef_inventory$row_exp_n)
  
  cat(sprintf("\nProjecting %d years (each year: %d simulations × 4 lines)...\n",
              N_YEARS, N_SIM))
  
  set.seed(5100)
  
  for (t in seq_len(N_YEARS)) {
    
    exp_growth_t  <- (1 + wtd_growth)^(t - 1)    # exposure multiplier
    inf_factor_t  <- (1 + inf_fwd)^(t - 1)        # severity inflation multiplier
    disc_factor_t <- 1 / (1 + r10yr_fwd)^t         # PV discount factor
    
    # Expected total loss for year t (used for premium pricing and investment float)
    E_S_t <- total_el_y1 * exp_growth_t * inf_factor_t
    
    # Gross premium (priced at start of year based on expected loss)
    P_t   <- E_S_t * load_factor
    
    # Expenses (deterministic % of premium)
    Exp_t <- P_t * expense_ratio
    
    # Investment income: 0.5 × E[loss] held on average as float, earned at r1yr
    # (premiums in at start, claims out throughout year --> ~half-year average float)
    Inv_t <- 0.5 * E_S_t * r1yr_fwd
    
    # --- Simulate each line for year t ---
    
    # BI: lognormal: shift meanlog by log(inflation)
    sim_t_bi <- simulate_aggregate(
      mu_pool   = bi_meanlog + log(inf_factor_t),
      lambda    = sum(bi_lambda) * exp_growth_t,
      N_SIM     = N_SIM, freq_dist = "poisson",
      sev_dist  = "lnorm", sigma = bi_sdlog, scale = 1e6
    )
    
    # WC: lognormal: same approach
    sim_t_wc <- simulate_aggregate(
      mu_pool   = wc_meanlog + log(inf_factor_t),
      lambda    = sum(wc_lambda) * exp_growth_t,
      N_SIM     = N_SIM, freq_dist = "poisson",
      sev_dist  = "lnorm", sigma = wc_sdlog, scale = 1000
    )
    
    # EF: Gamma: inflate mu_pool; weights unchanged (equipment mix stable)
    sim_t_ef <- simulate_aggregate(
      mu_pool   = ef_inventory$mu_sev * inf_factor_t,
      lambda    = sum(ef_inventory$row_exp_n) * exp_growth_t,
      N_SIM     = N_SIM, freq_dist = "poisson",
      prob_w    = ef_prob_w_base,
      sev_dist  = "gamma", phi = phi_ef
    )
    
    # Cargo: Gamma: use lambda-weighted combined pool, inflated
    sim_t_cargo <- simulate_aggregate(
      mu_pool   = cargo_mu_combined * inf_factor_t,
      lambda    = cargo_lam_total * exp_growth_t,
      N_SIM     = N_SIM, freq_dist = "negbin",
      theta     = cargo_theta,
      sev_dist  = "gamma", phi = phi_cargo
    )
    
    S_t  <- sim_t_bi + sim_t_wc + sim_t_ef + sim_t_cargo
    
    # Net revenue (per simulation): premium + investment income − claims − expenses
    NR_t <- P_t + Inv_t - S_t - Exp_t
    
    # Store
    proj_costs[, t]   <- S_t
    proj_premium[t]   <- P_t
    proj_net_rev[, t] <- NR_t
    
    proj_pv_cost <- proj_pv_cost + S_t  * disc_factor_t
    proj_pv_rev  <- proj_pv_rev  + NR_t * disc_factor_t
    
    yr_exp_cost[t]    <- mean(S_t)
    yr_exp_premium[t] <- P_t
    yr_exp_net_rev[t] <- mean(NR_t)
    
    cat(sprintf(
      "  Y%2d (%d): E[Cost]=Đ%s  Premium=Đ%s  InvInc=Đ%s  E[NetRev]=Đ%s\n",
      t, 2175 + t - 1,
      formatC(round(mean(S_t)), format="f", digits=0, big.mark=","),
      formatC(round(P_t),       format="f", digits=0, big.mark=","),
      formatC(round(Inv_t),     format="f", digits=0, big.mark=","),
      formatC(round(mean(NR_t)),format="f", digits=0, big.mark=",")))
  }
  
  # =======
  # 9G. SHORT-TERM RESULTS (Year 1)
  
  fmt_tbl_cost <- function(df, title) {
    cat(sprintf("\n=== %s ===\n", title))
    df2     <- df
    num_col <- sapply(df2, is.numeric) & !(names(df2) %in% c("CV","Pr_Loss"))
    df2[num_col] <- lapply(df2[num_col], function(x)
      formatC(round(x), format="f", digits=0, big.mark=","))
    if ("Pr_Loss" %in% names(df2))
      df2$Pr_Loss <- sprintf("%.1f%%", as.numeric(df$Pr_Loss) * 100)
    print(df2, row.names=FALSE)
  }
  
  fmt_tbl_cost(
    summarise_sim_cost(proj_costs[, 1], "Y1: Portfolio Costs"),
    "SHORT-TERM (Year 1): Aggregate Costs (Đ)"
  )
  
  fmt_tbl_cost(
    summarise_sim_rev(proj_net_rev[, 1], "Y1: Net Revenue"),
    "SHORT-TERM (Year 1): Net Revenue (Đ) [left-tail = bad outcomes]"
  )
  
  cat(sprintf("\nYear 1 (2175) Income Statement:\n"))
  cat(sprintf("  Gross Premium:             Đ%s\n",
              formatC(round(yr_exp_premium[1]), format="f", digits=0, big.mark=",")))
  cat(sprintf("  Expected Claims:           Đ%s\n",
              formatC(round(yr_exp_cost[1]),    format="f", digits=0, big.mark=",")))
  cat(sprintf("  Expenses (%.0f%% of prem):   Đ%s\n", expense_ratio*100,
              formatC(round(yr_exp_premium[1] * expense_ratio), format="f", digits=0, big.mark=",")))
  cat(sprintf("  Investment Income (float): Đ%s\n",
              formatC(round(0.5 * total_el_y1 * r1yr_fwd), format="f", digits=0, big.mark=",")))
  cat(sprintf("  Expected Net Revenue:      Đ%s\n",
              formatC(round(yr_exp_net_rev[1]), format="f", digits=0, big.mark=",")))
  cat(sprintf("  Combined Ratio:            %.1f%%\n",
              (yr_exp_cost[1] + yr_exp_premium[1] * expense_ratio) / yr_exp_premium[1] * 100))
  cat(sprintf("  Pr(Underwriting Loss Y1):  %.1f%%\n", mean(proj_net_rev[,1] < 0) * 100))
  
  # =======
  # 9H. LONG-TERM RESULTS (10-Year)
  
  cat("\n\n=== LONG-TERM ANNUAL TRAJECTORY (Nominal Đ, Expected Values) ===\n")
  traj <- data.frame(
    Year        = 2175:(2175 + N_YEARS - 1),
    E_Cost      = yr_exp_cost,
    Premium     = yr_exp_premium,
    E_Net_Rev   = yr_exp_net_rev,
    Cum_Cost    = cumsum(yr_exp_cost),
    Cum_Premium = cumsum(yr_exp_premium),
    Cum_Net_Rev = cumsum(yr_exp_net_rev)
  )
  traj_fmt <- traj
  traj_fmt[, -1] <- lapply(traj[, -1], function(x)
    formatC(round(x), format="f", digits=0, big.mark=","))
  print(traj_fmt, row.names=FALSE)
  
  # PV distributions
  fmt_tbl_cost(
    summarise_sim_cost(proj_pv_cost, "10-Yr PV Costs"),
    "LONG-TERM (10-Year PV): Aggregate Costs (Đ)"
  )
  
  fmt_tbl_cost(
    summarise_sim_rev(proj_pv_rev, "10-Yr PV Net Revenue"),
    "LONG-TERM (10-Year PV): Net Revenue (Đ) [left-tail = bad outcomes]"
  )
  
  cat(sprintf("\n10-Year Cumulative (Nominal):\n"))
  cat(sprintf("  Total Premium:            Đ%s\n",
              formatC(round(sum(yr_exp_premium)), format="f", digits=0, big.mark=",")))
  cat(sprintf("  Expected Total Cost:      Đ%s\n",
              formatC(round(sum(yr_exp_cost)),    format="f", digits=0, big.mark=",")))
  cat(sprintf("  Expected Net Revenue:     Đ%s\n",
              formatC(round(sum(yr_exp_net_rev)), format="f", digits=0, big.mark=",")))
  cat(sprintf("\n10-Year PV @ %.2f%% discount rate:\n", r10yr_fwd * 100))
  cat(sprintf("  PV(Costs)   Mean:          Đ%s\n",
              formatC(round(mean(proj_pv_cost)), format="f", digits=0, big.mark=",")))
  cat(sprintf("  PV(NetRev)  Mean:          Đ%s\n",
              formatC(round(mean(proj_pv_rev)),  format="f", digits=0, big.mark=",")))
  cat(sprintf("  PV(NetRev)  P05 (1-in-20): Đ%s\n",
              formatC(round(quantile(proj_pv_rev, 0.05)), format="f", digits=0, big.mark=",")))
  cat(sprintf("  PV(NetRev)  P01 (1-in-100):Đ%s\n",
              formatC(round(quantile(proj_pv_rev, 0.01)), format="f", digits=0, big.mark=",")))
  cat(sprintf("  Pr(PV Net Revenue < 0):    %.1f%%\n", mean(proj_pv_rev < 0) * 100))
  
  # =======
  # 9I. TAIL BEHAVIOUR SUMMARY: ALL TIME HORIZONS
  #     Addresses Objective 2c: tail behaviours for a variety of risks
  
  cat("\n\n=== TAIL BEHAVIOUR: COSTS (upper tail = bad) ===\n")
  cost_tail <- rbind(
    summarise_sim_cost(proj_costs[,  1], "Y1  Costs"),
    summarise_sim_cost(proj_costs[,  5], "Y5  Costs"),
    summarise_sim_cost(proj_costs[, 10], "Y10 Costs"),
    summarise_sim_cost(proj_pv_cost,     "10yr PV Costs")
  )
  cost_tail_fmt <- cost_tail %>%
    mutate(across(c(Mean, SD, Variance, P50, P75, P95, P99, TVaR99),
                  ~formatC(round(.), format="f", digits=0, big.mark=",")))
  print(cost_tail_fmt, row.names=FALSE)
  
  cat("\n=== TAIL BEHAVIOUR: NET REVENUE (lower tail = bad) ===\n")
  rev_tail <- rbind(
    summarise_sim_rev(proj_net_rev[,  1], "Y1  Net Revenue"),
    summarise_sim_rev(proj_net_rev[,  5], "Y5  Net Revenue"),
    summarise_sim_rev(proj_net_rev[, 10], "Y10 Net Revenue"),
    summarise_sim_rev(proj_pv_rev,        "10yr PV Net Revenue")
  )
  rev_tail_fmt <- rev_tail %>%
    mutate(across(c(Mean, SD, Variance, P50, P25, P05, P01, TVaR01),
                  ~formatC(round(.), format="f", digits=0, big.mark=","))) %>%
    mutate(Pr_Loss = sprintf("%.1f%%", Pr_Loss * 100))
  print(rev_tail_fmt, row.names=FALSE)
  
  # =======
  # 10. EXTREME SCENARIO STRESS TESTING (1-in-100 Year Event)
  # Scenario: "Carrington-Class Coronal Mass Ejection"
  # =======
  
  cat("\n--- Running 1-in-100 Year Stress Test ---\n")
  
  # 1. Define the Shock Multipliers based on qualitative risk assessment
  # Helionis is shielded but suffers supply chain shocks. Bayesia/Oryn get hit directly.
  freq_shocks <- c("Helionis Cluster" = 1.25, "Bayesia System" = 3.50, "Oryn Delta" = 3.00)
  sev_shocks  <- c("Helionis Cluster" = 1.15, "Bayesia System" = 2.50, "Oryn Delta" = 2.00)
  
  # Initialize lists to hold the stressed simulation results
  bi_sim_stress <- list()
  wc_sim_stress <- list()
  ef_sim_stress <- list()
  cargo_sim_stress <- list()
  
  for (sys in cq_systems) {
    # Extract specific system shocks
    f_shock <- freq_shocks[sys]
    s_shock <- sev_shocks[sys]
    
    # --- BI STRESS ---
    # Lognormal adjustment: add log(severity_shock) to meanlog
    bi_sim_stress[[sys]] <- simulate_aggregate(
      mu_pool = bi_meanlog + log(s_shock), 
      lambda = bi_lambda[sys] * f_shock, 
      N_SIM = N_SIM, freq_dist = "poisson", sev_dist = "lnorm", sigma = bi_sdlog, scale = 1e6
    )
    
    # --- WC STRESS ---
    wc_sim_stress[[sys]] <- simulate_aggregate(
      mu_pool = wc_meanlog + log(s_shock), 
      lambda = wc_lambda[sys] * f_shock, 
      N_SIM = N_SIM, freq_dist = "poisson", sev_dist = "lnorm", sigma = wc_sdlog, scale = 1000
    )
    
    # --- EF STRESS ---
    # Gamma adjustment: multiply mu_pool by severity_shock
    rows <- ef_inventory[ef_inventory$cq_system == sys, ]
    lam_sys_stressed <- sum(rows$row_exp_n) * f_shock
    prob_w <- rows$row_exp_n / sum(rows$row_exp_n) # Weights remain proportional
    
    ef_sim_stress[[sys]] <- simulate_aggregate(
      mu_pool = rows$mu_sev * s_shock, 
      lambda = lam_sys_stressed, 
      N_SIM = N_SIM, freq_dist = "poisson", prob_w = prob_w, sev_dist = "gamma", phi = phi_ef
    )
    
    # --- CARGO STRESS ---
    cargo_sim_stress[[sys]] <- simulate_aggregate(
      mu_pool = cargo_mu_pool[[sys]] * s_shock, 
      lambda = cargo_lambda[sys] * f_shock, 
      N_SIM = N_SIM, freq_dist = "negbin", theta = cargo_theta, sev_dist = "gamma", phi = phi_cargo
    )
  }
  
  # 2. Aggregate the Stressed Enterprise Loss
  enterprise_stress <- Reduce("+", bi_sim_stress) + Reduce("+", wc_sim_stress) + 
    Reduce("+", ef_sim_stress) + Reduce("+", cargo_sim_stress)
  
  # 3. Compare Baseline vs Stressed
  baseline_mean <- mean(full_portfolio)
  baseline_p99  <- quantile(full_portfolio, 0.99)
  
  stress_mean <- mean(enterprise_stress)
  stress_p99  <- quantile(enterprise_stress, 0.99)
  
  cat(sprintf("Baseline Expected Loss: %s\n", formatC(baseline_mean, format="f", big.mark=",", digits=0)))
  cat(sprintf("Baseline 1-in-100 (P99): %s\n", formatC(baseline_p99, format="f", big.mark=",", digits=0)))
  cat("--------------------------------------------------\n")
  cat(sprintf("STRESSED Expected Loss: %s\n", formatC(stress_mean, format="f", big.mark=",", digits=0)))
  cat(sprintf("STRESSED 1-in-100 (P99): %s\n", formatC(stress_p99, format="f", big.mark=",", digits=0)))
  
  # 4. Visualization (Plotting the Tail Risk)
  # We can use ggplot2 to create a highly professional overlapping density plot for the report
  plot_data <- data.frame(
    Loss = c(full_portfolio, enterprise_stress),
    Scenario = rep(c("Baseline", "Stressed (Solar Flare)"), each = N_SIM)
  )
  
  ggplot(plot_data, aes(x = Loss, fill = Scenario)) +
    geom_density(alpha = 0.5) +
    scale_x_log10(labels = scales::comma) +
    scale_fill_manual(values = c("Baseline" = "steelblue", "Stressed (Solar Flare)" = "firebrick")) +
    theme_minimal() +
    labs(title = "Aggregate Enterprise Loss: Baseline vs. 1-in-100 Year Stress Scenario",
         x = "Total Aggregate Loss",
         y = "Density",
         #subtitle = "Illustrating the catastrophic tail shift caused by a system-wide solar event."
    ) +
    theme(legend.position = "bottom")
  