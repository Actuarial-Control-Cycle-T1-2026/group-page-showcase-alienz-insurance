
# Cleans all 4 claims datasets (frequency and severity), so that every column satisfies its data dictionary definition

library(readxl)
library(dplyr)
library(stringr)

#========================================
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

# =================================
# 1. Business Interruption Dataset – Cleaning the frequency data 

bi_freq_raw <- read_excel(
  'srcsc-2026-claims-business-interruption.xlsx',
  sheet = 'freq'
)
 
bi_freq <- bi_freq_raw |>
  mutate(
    policy_id  = clean_id(coalesce(policy_id, "")),
    policy_id  = if_else(grepl("^BI-[0-9]{6}$", policy_id), policy_id, NA_character_),
    station_id = clean_charac(coalesce(station_id, "")),
    station_id = if_else(station_id == "", NA_character_, station_id),
    
    # cleaning and ensuring valid solar systems: {Helionis Cluster, Epsilon, Zeta}
    solar_system = clean_charac(coalesce(solar_system, "")),
    solar_system = if_else(
      solar_system %in% c("Helionis Cluster", "Epsilon", "Zeta"),
      solar_system, NA_character_),
    
    # cleaning production load (ratio 0-1), there are many value below 0 and above 1
    production_load = clean_numeric(production_load, lo = 0, hi = 1),
    
    # energy backup score (ordered category 1-5)
    # issue is that, there are float values and negatives
    # hence, will round to nearest int, and see if it is within {1,2,3,4,5}, or else rest -> NA
    energy_backup_score = clean_levels(energy_backup_score, 1:5),
    
    # Supply chain index (ratio 0 to 1), there are many values below 0 and above 1 
    supply_chain_index = clean_numeric(supply_chain_index, lo = 0, hi = 1),
    
    # Average crew experience (between 1-30 years),
    avg_crew_exp = clean_numeric(avg_crew_exp, lo = 1, hi = 30),
    
    # Maintenance frequency (between 0 to 6), 
    maintenance_freq = clean_numeric(maintenance_freq, lo = 0, hi = 6),
    
    # Safety compliance (ordered category 1-5), 
    safety_compliance = clean_levels(safety_compliance, 1:5),
    
    # exposure (ratio between 0 to 1), 
    exposure = clean_numeric(exposure, lo = 0, hi = 1),
    
    # claim count (integer 0 to 4)
    claim_count = clean_levels(claim_count, 0:4)
  ) |>
  # Drop rows where exposure is NA (because we can't compute rate without exposure)
  filter(!is.na(exposure))

summary(bi_freq)

#number of rows removed, due to NA exposure 
#nrow(bi_freq_raw) - nrow(bi_freq)

# ============================
# 2. Business Interruption - cleaning the severity data 

bi_sev_raw <- read_excel(
  'srcsc-2026-claims-business-interruption.xlsx',
  sheet = 'sev'
)

bi_sev <- bi_sev_raw |>
  mutate(
    # cleaning id's 
    claim_id  = clean_id(coalesce(claim_id, "")),
    claim_id  = if_else(grepl("^BI-C-[0-9]+$", claim_id), claim_id, NA_character_),
    policy_id = clean_id(coalesce(policy_id, "")),
    policy_id = if_else(grepl("^BI-[0-9]{6}$", policy_id), policy_id, NA_character_),
    station_id = clean_charac(coalesce(station_id, "")),
    station_id = if_else(station_id == "", NA_character_, station_id),
    
    # Claim seq should be positive integer 1,2,3,4, but there are negative values (-3,-2,-1) and floats (6.5, 41.6)
    claim_seq = as.integer(round(claim_seq)),
    claim_seq = if_else(claim_seq >= 1 & claim_seq <= 4, claim_seq, NA_integer_),
    
    # solar systems 
    solar_system = clean_charac(coalesce(solar_system, "")),
    solar_system = if_else(
      solar_system %in% c("Helionis Cluster", "Epsilon", "Zeta"),
      solar_system, NA_character_
    ),
    
    # production load  
    production_load = clean_numeric(production_load, lo = 0, hi = 1),
    
    # energy backup score  
    energy_backup_score = clean_levels(energy_backup_score, 1:5),
    # exposure (ratio between 0 to 1), 
    exposure = clean_numeric(exposure, lo = 0, hi = 1),
    # safety compliance  
    safety_compliance = clean_levels(safety_compliance, 1:5),
    
    # claim amount: dict states ~28K – 1,426K [FIX 1: replaced positivity check with range bounds]
    claim_amount = clean_numeric(claim_amount, lo = 28000, hi = 1426000)
  ) |>
  # for analysis, we need to have valid positive claim amounts 
  filter(!is.na(claim_amount))

#nrow(bi_sev_raw) - nrow(bi_sev)

# ====================================
# Cargo loss data - cleaning the frequency data 

cargo_freq_raw <- read_excel(
  'srcsc-2026-claims-cargo.xlsx',
  sheet = 'freq'
)

cargo_freq <- cargo_freq_raw |>
  mutate(
    policy_id   = clean_id(coalesce(policy_id, "")),
    policy_id   = if_else(grepl("^CL-[0-9]{6}$", policy_id), policy_id, NA_character_),
    shipment_id = clean_id(coalesce(shipment_id, "")),
    shipment_id = if_else(grepl("^S-[0-9]{6}$", shipment_id), shipment_id, NA_character_),
    
    # cargo type -- CARGO TYPE (categories: cobalt, gold, lithium, platinum, rare earths,
    #                            supplies, titanium) --
    cargo_type = clean_charac(coalesce(cargo_type, "")),
    cargo_type = if_else(
      cargo_type %in% c("cobalt","gold","lithium","platinum",
                        "rare earths","supplies","titanium"),
      cargo_type, NA_character_
    ),
    
    # cargo value: dict states ~50K – 680,000K [FIX 2: replaced positivity check with range bounds]
    cargo_value = clean_numeric(cargo_value, lo = 50000, hi = 680000000),
    
    # weight, cleaning the negative weights
    #also according to the data dictionary, valid weights are 1.5K – 250K
    weight = clean_numeric(weight, lo = 1500, hi = 250000),
    
    # route risk (ordered category 1 to 5) 
    route_risk = clean_levels(route_risk, 1:5),
    
    # distance (according to data dictionary, the range should only go up to 100) 
    distance = clean_numeric(distance, lo = 0.01, hi = 100),
    
    # Transit duration (must range from 1 to 60 months)
    transit_duration = clean_numeric(transit_duration, lo = 1, hi = 60),
    
    # pilot experience (1 to 30 years) 
    pilot_experience = clean_numeric(pilot_experience, lo = 1, hi = 30),
    
    # vessel age (1 to 50 years) 
    vessel_age = clean_numeric(vessel_age, lo = 1, hi = 50),
    
    # container type 
    # In the dataset, the valid types observed are: DeepSpace Haulbox, DockArc Freight Case,
    # HardSeal Transit Crate, LongHaul Vault Canister, QuantumCrate Module
    container_type = clean_charac(coalesce(container_type, "")),
    container_type = if_else(
      container_type %in% c("DeepSpace Haulbox", "DockArc Freight Case",
                            "HardSeal Transit Crate",
                            "LongHaul Vault Canister", "QuantumCrate Module"),
      container_type, NA_character_
    ),
    
    # solar radiation (must range from 0 to 1) 
    solar_radiation = clean_numeric(solar_radiation, lo = 0, hi = 1),
    
    # debris density (must range from 0 to 1) 
    debris_density = clean_numeric(debris_density, lo = 0, hi = 1),
    
    # exposure (0-1, >0) 
    exposure = clean_numeric(exposure, lo = 1e-6, hi = 1),
    
    # claim count (should be integer 0-5), but there are negative values (-3,-2,-1) and floats (5.06, 34.6 etc.)
    claim_count = clean_levels(claim_count, 0:5)
  ) |>
  filter(!is.na(exposure))

summary(cargo_freq)

# ===============================
# cargo loss data - cleaning severity 

cargo_sev_raw <- read_excel(
  'srcsc-2026-claims-cargo.xlsx',
  sheet = 'sev'
)

cargo_sev <- cargo_sev_raw |>
  mutate(
    claim_id    = clean_id(coalesce(claim_id, "")),
    claim_id    = if_else(grepl("^CAR-C-[0-9]+$", claim_id), claim_id, NA_character_),
    policy_id   = clean_id(coalesce(policy_id, "")),
    policy_id   = if_else(grepl("^CL-[0-9]{6}$", policy_id), policy_id, NA_character_),
    shipment_id = clean_id(coalesce(shipment_id, "")),
    shipment_id = if_else(grepl("^S-[0-9]{6}$", shipment_id), shipment_id, NA_character_),
    
    # cargo type 
    cargo_type = clean_charac(coalesce(cargo_type, "")),
    cargo_type = if_else(
      cargo_type %in% c("cobalt","gold","lithium","platinum",
                        "rare earths","supplies","titanium"),
      cargo_type, NA_character_
    ),
    
    # weight, cleaning the negative weights
    #also according to the data dictionary, valid weights are 1.5K – 250K
    weight = clean_numeric(weight, lo = 1500, hi = 250000),
    
    # cargo value: dict states ~50K – 680,000K [FIX 3: replaced positivity check with range bounds]
    cargo_value = clean_numeric(cargo_value, lo = 50000, hi = 680000000),
    
    # route risk 
    route_risk = clean_levels(route_risk, 1:5),
    
    # distance  
    distance = clean_numeric(distance, lo = 0.01, hi = 100),
    
    # transit duration  
    transit_duration = clean_numeric(transit_duration, lo = 1, hi = 60),
    
    # pilot experience  
    pilot_experience = clean_numeric(pilot_experience, lo = 1, hi = 30),
    
    # vessel age  
    vessel_age = clean_numeric(vessel_age, lo = 1, hi = 50),
    
    # container type  
    container_type = clean_charac(coalesce(container_type, "")),
    container_type = if_else(
      container_type %in% c("DeepSpace Haulbox", "DockArc Freight Case",
                            "HardSeal Transit Crate",
                            "LongHaul Vault Canister", "QuantumCrate Module"),
      container_type, NA_character_
    ),
    
    # solar radiation  
    solar_radiation = clean_numeric(solar_radiation, lo = 0, hi = 1),
    
    # exposure (0-1, >0) 
    exposure = clean_numeric(exposure, lo = 0, hi = 1),
    
    # debris density  
    debris_density = clean_numeric(debris_density, lo = 0, hi = 1),
    
    # claim amount: dict states ~31K – 678,000K [FIX 4: replaced positivity check with range bounds]
    claim_amount = clean_numeric(claim_amount, lo = 31000, hi = 678000000)
  ) |>
  filter(!is.na(claim_amount))

summary(cargo_sev)

# =========================================
# Equipment failure - cleaning the frequency data 

ef_freq_raw <- read_excel(
  'srcsc-2026-claims-equipment-failure.xlsx',
  sheet = 'freq'
)

# Standardising equipment type names to match inventory naming
equip_name_map <- c(
  "FexStram Carrier"  = "Fexstram Carrier",
  "Flux Rider"        = "Flux Rider",
  "Graviton Extractor"= "Graviton Extractor",
  "Ion Pulverizer"    = "Ion Pulverizer",
  "Quantum Bore"      = "Quantum Bore",
  "ReglAggregators"   = "ReglAggregators"
)

ef_freq <- ef_freq_raw |>
  mutate(
    policy_id    = clean_id(coalesce(policy_id, "")),
    policy_id    = if_else(grepl("^EF-[0-9]{6}$", policy_id), policy_id, NA_character_),
    equipment_id = clean_id(coalesce(equipment_id, "")),
    equipment_id = if_else(grepl("^EQ-[0-9]{6}$", equipment_id), equipment_id, NA_character_),
    
    # equipment type, cleaning data and there are also name mismatches with inventory 
    equipment_type = clean_charac(coalesce(equipment_type, "")),
    equipment_type = coalesce(equip_name_map[equipment_type], equipment_type),
    equipment_type = if_else(
      equipment_type %in% names(equip_name_map) |
        equipment_type %in% equip_name_map,
      equipment_type, NA_character_
    ),
    
    # equipment age (0-10 Earth years) 
    # however, there are 144 negative values, and 47,387 values above 10 (up to 327)
    # very large data quality issue in the EF dataset
    equipment_age = clean_numeric(equipment_age, lo = 0, hi = 10),
    
    # solar system 
    solar_system = clean_charac(coalesce(solar_system, "")),
    solar_system = if_else(
      solar_system %in% c("Helionis Cluster", "Epsilon", "Zeta"),
      solar_system, NA_character_
    ),
    
    # maintenance interval (100-5000 hours)
    maintenance_int = clean_numeric(maintenance_int, lo = 100, hi = 5000),
    
    # usage intensity (0-24)
    usage_int = clean_numeric(usage_int, lo = 0, hi = 24),
    
    # exposure  
    exposure = clean_numeric(exposure, lo = 0, hi = 1),
    
    # claim count (integer 0-3) --
    claim_count = clean_levels(claim_count, 0:3)
  ) |>
  filter(!is.na(exposure))

# =============================================================================
# equipment failure - cleaning the severity data 

ef_sev_raw <- read_excel(
  'srcsc-2026-claims-equipment-failure.xlsx',
  sheet = 'sev'
)

ef_sev <- ef_sev_raw |>
  mutate(
    claim_id     = clean_id(coalesce(claim_id, "")),
    claim_id     = if_else(grepl("^EF-C-[0-9]+$", claim_id), claim_id, NA_character_),
    policy_id    = clean_id(coalesce(policy_id, "")),
    policy_id    = if_else(grepl("^EF-[0-9]{6}$", policy_id), policy_id, NA_character_),
    equipment_id = clean_id(coalesce(equipment_id, "")),
    equipment_id = if_else(grepl("^EQ-[0-9]{6}$", equipment_id), equipment_id, NA_character_),
    
    # equipment type  
    equipment_type = clean_charac(coalesce(equipment_type, "")),
    equipment_type = coalesce(equip_name_map[equipment_type], equipment_type),
    equipment_type = if_else(
      equipment_type %in% names(equip_name_map) |
        equipment_type %in% equip_name_map,
      equipment_type, NA_character_
    ),
    
    # equipment age 
    equipment_age = clean_numeric(equipment_age, lo = 0, hi = 10),
    
    # solar system 
    solar_system = clean_charac(coalesce(solar_system, "")),
    solar_system = if_else(
      solar_system %in% c("Helionis Cluster", "Epsilon", "Zeta"),
      solar_system, NA_character_
    ),
    
    # maintenance interval 
    maintenance_int = clean_numeric(maintenance_int, lo = 100, hi = 5000),
    
    # usage intensity 
    usage_int = clean_numeric(usage_int, lo = 0, hi = 24),
    
    # exposure  
    exposure = clean_numeric(exposure, lo = 0, hi = 1),
    
    # claim amount: dict states ~11K – 790K [FIX 5: replaced positivity check with range bounds]
    claim_amount = clean_numeric(claim_amount, lo = 11000, hi = 790000)
  ) |>
  filter(!is.na(claim_amount))

# ======================================= 
# workers comp data - cleaning frequency 

wc_freq_raw <- read_excel(
  'srcsc-2026-claims-workers-comp.xlsx',
  sheet = 'freq'
)

valid_occupations <- c("Administrator", "Drill Operator", "Engineer", "Executive",
                       "Maintenance Staff", "Manager", "Planetary Operations",
                       "Safety Officer", "Scientist", "Spacecraft Operator",
                       "Technology Officer")

wc_freq <- wc_freq_raw |>
  mutate(
    policy_id = clean_id(coalesce(policy_id, "")),
    policy_id = if_else(grepl("^WC-[A-Z]+-[0-9]+$", policy_id), policy_id, NA_character_),
    worker_id = clean_id(coalesce(worker_id, "")),
    worker_id = if_else(grepl("^W-[0-9]+$", worker_id), worker_id, NA_character_),
    station_id = clean_charac(coalesce(station_id, "")),
    station_id = if_else(station_id == "", NA_character_, station_id),
    
    solar_system = clean_charac(coalesce(solar_system, "")),
    solar_system = if_else(
      solar_system %in% c("Helionis Cluster", "Epsilon", "Zeta"),
      solar_system, NA_character_
    ),
    
    occupation = clean_charac(coalesce(occupation, "")),
    occupation = if_else(occupation %in% valid_occupations, occupation, NA_character_),
    
    # employemnt type -
    employment_type = clean_charac(coalesce(employment_type, "")),
    employment_type = if_else(
      employment_type %in% c("Full-time", "Contract"),
      employment_type, NA_character_
    ),
    
    # expereince years (from 0.2 to 40)
    experience_yrs = clean_numeric(experience_yrs, lo = 0.1, hi = 40),
    
    # accident history
    accident_history_flag = as.integer(round(accident_history_flag)),
    accident_history_flag = if_else(
      accident_history_flag %in% 0:1,
      accident_history_flag, NA_integer_
    ),
    
    # psych stress index (ordered category 1-5) 
    psych_stress_index = clean_levels(psych_stress_index, 1:5),
    
    # hours per week (20,25,30,35,40)
    hours_per_week = as.integer(round(hours_per_week)),
    hours_per_week = if_else(
      hours_per_week %in% c(20L, 25L, 30L, 35L, 40L),
      hours_per_week, NA_integer_
    ),
    
    # supervision level (ratio 0 to 1) 
    supervision_level = clean_numeric(supervision_level, lo = 0, hi = 1),
    
    # gravity level (0.75 to 1.50) 
    gravity_level = clean_numeric(gravity_level, lo = 0.75, hi = 1.50),
    
    # safety training index (ordered category 1 to 5) --
    safety_training_index = clean_levels(safety_training_index, 1:5),
    
    # protective gear quality (ordered category 1-5)
    protective_gear_quality = clean_levels(protective_gear_quality, 1:5),
    
    # Base salary, according to data dictionary should be from 20k to 130k 
    base_salary = clean_numeric(base_salary, lo = 20000, hi = 130000),
    
    # exposure 
    exposure = clean_numeric(exposure, lo = 1e-6, hi = 1),
    
    # claim count (integer 0-2)
    claim_count = clean_levels(claim_count, 0:2)
  ) |>
  filter(!is.na(exposure))

# ==========================================
# workers compensation - cleaning severity data  

wc_sev_raw <- read_excel(
  'srcsc-2026-claims-workers-comp.xlsx',
  sheet = 'sev'
)

valid_injury_types  <- c("Amputation","Burns","Cut laceration",
                         "Other","Psychological","Sprain, strain","Stress")
valid_injury_causes <- c("Caught in machine","Exposure","Other",
                         "Overexertion","Stress/strain","Vehicle accident","Violence")

wc_sev <- wc_sev_raw |>
  mutate(
    claim_id  = clean_id(coalesce(claim_id, "")),
    claim_id  = if_else(grepl("^C-[0-9]+$", claim_id), claim_id, NA_character_),
    policy_id = clean_id(coalesce(policy_id, "")),
    policy_id = if_else(grepl("^WC-[A-Z]+-[0-9]+$", policy_id), policy_id, NA_character_),
    worker_id = clean_id(coalesce(worker_id, "")),
    worker_id = if_else(grepl("^W-[0-9]+$", worker_id), worker_id, NA_character_),
    station_id = clean_charac(coalesce(station_id, "")),
    station_id = if_else(station_id == "", NA_character_, station_id),
    
    solar_system = clean_charac(coalesce(solar_system, "")),
    solar_system = if_else(
      solar_system %in% c("Helionis Cluster", "Epsilon", "Zeta"),
      solar_system, NA_character_
    ),
    
    occupation = clean_charac(coalesce(occupation, "")),
    occupation = if_else(occupation %in% valid_occupations, occupation, NA_character_),
    
    employment_type = clean_charac(coalesce(employment_type, "")),
    employment_type = if_else(
      employment_type %in% c("Full-time", "Contract"),
      employment_type, NA_character_
    ),
    
    experience_yrs = clean_numeric(experience_yrs, lo = 0.1, hi = 40),
    
    accident_history_flag = as.integer(round(accident_history_flag)),
    accident_history_flag = if_else(
      accident_history_flag %in% 0:1, accident_history_flag, NA_integer_
    ),
    
    psych_stress_index = clean_levels(psych_stress_index, 1:5),
    
    hours_per_week = as.integer(round(hours_per_week)),
    hours_per_week = if_else(
      hours_per_week %in% c(20L, 25L, 30L, 35L, 40L),
      hours_per_week, NA_integer_
    ),
    
    supervision_level = clean_numeric(supervision_level, lo = 0, hi = 1),
    
    gravity_level = clean_numeric(gravity_level, lo = 0.75, hi = 1.50),
    
    safety_training_index = clean_levels(safety_training_index, 1:5),
    
    protective_gear_quality = clean_levels(protective_gear_quality, 1:5),
    
    base_salary = clean_numeric(base_salary, lo = 20000, hi = 130000), # [FIX 6: corrected lo=1000 -> lo=20000, dict states ~20K–130K]
    
    #injury type 
    injury_type = clean_charac(coalesce(injury_type, "")),
    injury_type = if_else(injury_type %in% valid_injury_types, injury_type, NA_character_),
    
    # injury cause
    injury_cause = clean_charac(coalesce(injury_cause, "")),
    injury_cause = if_else(injury_cause %in% valid_injury_causes, injury_cause, NA_character_),
    
    # claim length 
    claim_length = clean_numeric(claim_length, lo = 1, hi = 1000),
    
    # exposure  
    exposure = clean_numeric(exposure, lo = 0, hi = 1),
    
    # claim amount 
    claim_amount = clean_numeric(claim_amount, lo = 0.01, hi = 200000))  |> filter(!is.na(claim_amount))
    
# ============================
# Saving cleaned datasets 
output <- "./outputs/cleaned"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

write.csv(bi_freq,    file.path(output, "bi_freq_clean.csv"))
write.csv(bi_sev,     file.path(output, "bi_sev_clean.csv"))
write.csv(cargo_freq, file.path(output, "cargo_freq_clean.csv"))
write.csv(cargo_sev,  file.path(output, "cargo_sev_clean.csv"))
write.csv(ef_freq,    file.path(output, "ef_freq_clean.csv"))
write.csv(ef_sev,     file.path(output, "ef_sev_clean.csv"))
write.csv(wc_freq,    file.path(output, "wc_freq_clean.csv"))
write.csv(wc_sev,     file.path(output, "wc_sev_clean.csv"))
