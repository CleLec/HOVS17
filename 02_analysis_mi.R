source("helper_functions.R")
load(file.path(dirs$data, "gpdata.Rdata"))
load(file.path(dirs$results, "all_results.Rdata"))





# Code ------------------------------------------------------------------------------


# Load requisite packages
library(tidyverse)
library(furrr)
library(lavaan)
library(broom)
library(parameters)
library(sirt)
library(ggplot2)
library(kableExtra)
library(gt)
library(webshot)
library(scales)
library(readr)
#library(reshape2)

# Setup plan for multi-processor workflows with furrr
plan(multisession, workers = 4)


# Define a function that reorders the items so that the positively worded ones
# come first (facilitates model convergence and interpretation)


# Create lookup tables 

lookup_mi  <- cfa_results %>% 
  filter(!(Dimension == "transcendence_t1"  & Model == "congeneric") ,  
         !(Dimension == "conservation_t1"  & Model == "congeneric"),   
           ) %>% 
  select(-mod_output, -mod_results) %>% 
  crossing( type = c("configural", "metric", "scalar"),
            group = c("gender", "age_group", "edu_group")                                                             
  ) %>% 
select(dimension = Dimension, Model, items = Items, type, form, group)

# Function for MI analysis ------------------------------------------------
# Functions for creating 


mi_mod <- function(items, type, form, group){
  
  if (type == "configural"){
    mod_output <- cfa(model = form, 
                      data = dat, 
                      orthogonal = TRUE,
                      group = group,
                      missing = 'fiml', 
                      estimator = "MLR") 
  }
  else if (type == "metric"){
    mod_output <- cfa(model = form, 
                      data = dat, 
                      orthogonal = TRUE,
                      group = group,
                      group.equal = c("loadings", "regressions"),
                      missing = 'fiml', 
                      estimator = "MLR") 
  } 
  else if (type == "scalar") {
    mod_output <- cfa(model = form, 
                      meanstructure = TRUE,
                      data = dat, 
                      orthogonal = TRUE,
                      group = group,
                      group.equal = c("loadings", "intercepts", "regressions"),
                      missing = 'fiml', 
                      estimator = "MLR") 
  } 
  mod_output
}

# Functions to evaluate model output
mi_evaluate <- function(mod_output) {
  mi_evaluation <- list(
    mi_converged = lavInspect(mod_output, "converged"),
    mi_iterations = lavInspect(mod_output, "iterations"),
    mi_fit = fitmeasures(mod_output) %>% 
      as_tibble_row() ,
    mi_parameters = model_parameters(mod_output)
  )
  mi_evaluation
}

## -------------------------------------------------------------------------- ##
##  Determine Measurement Invariance Level based on Fit Indices - Chen, 2007  ##
## -------------------------------------------------------------------------- ##

### FUNCTIONS config_inv(), metric_inv(), scalar_inv() ###

# Determine fit for configural model
config_inv <- function(type, cfi.scaled, rmsea.scaled, srmr) {
  if (type == "configural") {
    ## RMSEA should be set at <= .080
    ## To get started on MI testing, a fitting configural type is required.
    ## Thus, RMSEA set to <= .0849 to achieve configural
    ## NOTE: There is a safety-switch in metric &scalar MI testing for RMSEA > 0.08
    config.inv <- ifelse(sum(c(cfi.scaled >= 0.9, rmsea.scaled <= 0.0849, srmr <= 0.08)) > 1, # 2 of three
                         TRUE, FALSE)
  } else { 
    config.inv <- FALSE 
  }
  config.inv
}

# Determine fit for metric type
metric_inv <- function(type, lag.config, cfi.diff, rmsea.diff, srmr.diff,
                       cfi.scaled, rmsea.scaled, srmr) {
  if ((type == "metric") & (lag.config == TRUE)) {
    metric.inv <- ifelse(((cfi.diff < -0.01 & (rmsea.diff > 0.015 | srmr.diff > 0.030))
                          | (cfi.scaled < 0.9 | rmsea.scaled > 0.08 | srmr > 0.08)), 
                         FALSE, TRUE)
  } else { 
    metric.inv <- FALSE 
  }
  metric.inv
}

# Determine fit for scalar type
scalar_inv <- function(type, lag.metric, cfi.diff, rmsea.diff, srmr.diff,
                       cfi.scaled, rmsea.scaled, srmr) {
  if (!is.na(lag.metric)) {
    if ((type == "scalar") & (lag.metric == TRUE)) {
      scalar.inv <- ifelse((cfi.diff < -0.01 & (rmsea.diff > 0.015 | srmr.diff > 0.010))
                           | (cfi.scaled < 0.9 | rmsea.scaled > 0.08 | srmr > 0.08), 
                           FALSE, TRUE)
    } else { 
      scalar.inv <- FALSE 
    }
  } else { 
    scalar.inv <- FALSE 
  }
  scalar.inv
}


### FUNCTIONS level_chen(), high_level_chen()

# Integrate levels attained in one column
level_chen <- function(type, config.inv, metric.inv, scalar.inv) {
  if (type == "configural") {
    level.chen <- config.inv
  } else if (type == "metric") {
    level.chen <- metric.inv
  } else if (type == "scalar") {
    level.chen <- scalar.inv
  }
  level.chen
}

# Determine best fitting type
high_level_chen <- function(type, level.chen, lag.level.chen, lag.lag.level) {
  if (type == "scalar" & level.chen == TRUE) {
    high.level <- "scalar"
  } else if (type == "scalar" & level.chen == FALSE & lag.level.chen == TRUE) {
    high.level <- "metric"
  } else if (type == "scalar" & level.chen == FALSE & lag.level.chen == FALSE & 
             lag.lag.level == TRUE) {
    high.level <- "configural"
  } else if (type == "scalar" & level.chen == FALSE & lag.level.chen == FALSE & 
             lag.lag.level == FALSE) {
    high.level <- "-"
  } else {
    high.level <- "-"
  }
  level.chen
}


## -------------------------------------------------------------------------- ##
##  Determine Measurement Invariance Level based on BIC metric - Raftery,1995 ##
## -------------------------------------------------------------------------- ##

### FUNCTIONS bic_metric(), bic_scalar() ###

# Determine parsimony-accuracy for metric type
bic_metric <- function(type, bic.diff) {
  if (type == "metric") {
    bic.metric.inv <- ifelse(bic.diff <= 0, TRUE, FALSE)
  } else {
    bic.metric.inv <- FALSE
  }
}

# Determine parsimony-accuracy for scalar type
bic_scalar <- function(type, bic.diff, lag.bic.metric) {
  if ((type == "scalar") & (lag.bic.metric == TRUE)) {
    bic.scalar.inv <- ifelse(bic.diff <= 0, TRUE, FALSE)
  } else {
    bic.scalar.inv <- FALSE
  }
}

### FUNCTIONS level_bic(), high_level_bic()

# Integrate parsimony-accuracy in one column
level_bic <- function(type, bic.metric, bic.scalar) {
  if (type == "metric") {
    level.bic <- bic.metric
  } else if (type == "scalar") {
    level.bic <- bic.scalar
  } else if (type == "configural") {
    level.bic <- "-"
  } else {
    level.bic <- FALSE
  }
}

# Determine best fitting type
high_level_bic <- function(type, level.bic, lag.level.bic) {
  if (type == "scalar" & level.bic == TRUE) {
    high.level <- "scalar"
  } else if (type == "scalar" & level.bic == FALSE & lag.level.bic == TRUE) {
    high.level <- "metric"
  } else if (type == "scalar" & level.bic == FALSE & lag.level.bic == FALSE) {
    high.level <- "-"
  } else {
    high.level <- "-"
  }
}


## -------------------------------------------------------------------------- ##
##  Determine Measurement Invariance Level Combining Fit Heuristics & BIC     ##
## -------------------------------------------------------------------------- ##

# Arrive at Integrative Decision on Measurement Invariance
int_dec <- function(type, high.level.chen, high.level.bic, config.total) {
  if (type == "scalar" & (high.level.chen == "scalar" | 
                          high.level.bic  == "scalar")) {
    int.dec <- "scalar"
  } else if (type == "scalar" & (high.level.chen == "metric" | 
                                 (high.level.bic  == "metric" & 
                                  config.total == TRUE))) {
    int.dec <- "metric"
  } else if (type == "scalar" & high.level.chen == "configural") {
    int.dec <- "configural"
  } else {
    int.dec <- "-"
  }
}


### FUNCTION extract_fit() ###
##  Create fit.table for Multi-Group CFA results ##
extract_fit <- function(input) {
  fit.table <- input %>% 
    unnest_wider(mod_results) %>% 
    select(-c(items, form, mod_output, mi_converged, 
              mi_iterations, mi_parameters))
  fit.table <- fit.table %>% 
    unnest_wider(mi_fit) %>% 
    select(dimension, group, type, ntotal, chisq.scaled, df.scaled, 
           pvalue.scaled, chisq.scaling.factor, cfi.scaled, 
           rmsea.scaled, srmr, bic)
  fit.table
}

### FUNCTION extract_delta_fit() ###
##  Expand fit.table with info about fit-differences between invariance types
extract_delta_fit <- function(input) {
  fit.table <- input %>% 
    group_by(dimension, group) %>% 
    mutate(  cfi.diff   = c(0, diff(cfi.scaled)),
             rmsea.diff = c(0, diff(rmsea.scaled)), 
             srmr.diff  = c(0, diff(srmr)),
             bic.diff   = c(0, diff(bic)))
  fit.table
}

### Function extract_invariance() ###
## create a fit table with invariance level

extract_invariance <- function(input) {
  fit.table <- input %>% 
    group_by(dimension, group) %>% 
    mutate(config.inv = pmap_lgl(list(type, cfi.scaled, rmsea.scaled, 
                                      srmr), 
                                 config_inv),
           lag.config = lag(config.inv, default = NA),
           metric.inv = pmap_lgl(list(type, lag.config, cfi.diff, 
                                      rmsea.diff, srmr.diff, cfi.scaled, 
                                      rmsea.scaled, srmr),
                                 metric_inv),
           lag.metric = lag(metric.inv, default = NA),
           scalar.inv = pmap_lgl(list(type, lag.metric, cfi.diff, 
                                      rmsea.diff, srmr.diff, cfi.scaled, 
                                      rmsea.scaled, srmr), 
                                 scalar_inv),
           bic.metric = pmap_lgl(list(type, bic.diff), bic_metric),
           lag.bic.metric = lag(bic.metric, default = NA),
           bic.scalar = pmap_lgl(list(type, bic.diff, lag.bic.metric), 
                                 bic_scalar),
           level.chen = pmap_chr(list(type, config.inv, metric.inv, 
                                      scalar.inv), 
                                 level_chen),
           level.bic = pmap_chr(list(type, bic.metric, bic.scalar), 
                                level_bic),
           lag.level.chen = lag(level.chen, default = NA),
           lag.lag.level = lag(level.chen, 2, default = NA),
           high.level.chen = pmap_chr(list(type, level.chen, 
                                           lag.level.chen, lag.lag.level), 
                                      high_level_chen),
           lag.level.bic = lag(level.bic, default = NA),
           high.level.bic = pmap_chr(list(type, level.bic, lag.level.bic), 
                                     high_level_bic),
           config.total    = lag(config.inv, 2, default = NA)) %>% 
    #  replace_na(list(level.chen = FALSE, level.BIC = FALSE)) %>% 
    select(-contains("lag.")) %>% 
    mutate(int.dec = pmap(list(type, high.level.chen, high.level.bic, 
                               config.total),
                          int_dec))
  
  fit.table
}


### FUNCTION extract_mi_level() ###
##  Create table with highest invariance levels attained from fit.table
extract_mi_level <- function(input) {
  highest.level <- input %>% 
    filter(type == "scalar") %>% 
    select(-c(type, config.total)) %>% 
    select(dimension, ntotal, int.dec)
  highest.level
}


### FUNCTION extract_parameters() ###
##  Extract type parameters from mg_cfa tibble with results
#   Note: input used = mi_results
extract_params <- function(input) {
  parameters <- input %>% 
    select(-c(items, group, form, mod_output)) 
  parameters <- parameters %>% 
    mutate(param = map(mod_results, last)) %>% 
    select(-c(mod_results))
  parameters
}



# Run the MI analyses ---------------------------------------------------------------

# Create reduced data
dat <- gpdata[c(t1_value_items, "gender", "age_group", "edu_group")]

# Running MI analysis
mi_results <- lookup_mi %>% 
  mutate(mod_output = future_pmap(list(items, type, form, group), 
                                  mi_mod),
         mod_results = future_map(mod_output, possibly(mi_evaluate, otherwise = NA)) 
  )


# Summarize the MI results ----------------------------------------------------------

# Determine measurement invariance (Multi-Group CFA)
fit.table.mi <- mi_results %>% 
  #group_by(dimension, group) %>% 
  extract_fit() %>% 
  extract_delta_fit()  %>% 
  extract_invariance()

traditionalMI_tbl <- fit.table.mi %>% select(dimension, type, ntotal, chisq.scaled, 
                                             pvalue.scaled, cfi.scaled, rmsea.scaled, 
                                             srmr, bic, level.chen, level.bic) %>%  
  rename(Dimension = dimension, Level_Tested = type, N = ntotal, ChiSquare = chisq.scaled, 
         pvalue = pvalue.scaled, CFI = cfi.scaled, RMSEA = rmsea.scaled, SRMR = srmr, 
         BIC = bic, Level_Chen = level.chen, Level_BIC = level.bic)  %>%  
  mutate(Level_Chen = case_when(Level_Tested == "configural" ~ "", TRUE ~ Level_Chen),
         Level_BIC = case_when(Level_Tested == "configural" ~ "", TRUE ~ Level_BIC),
         Level_achieved = case_when(Level_Tested == "scalar" & Level_Chen == TRUE & Level_BIC == TRUE  ~ "scalar",
                                    Level_Tested == "metric" & Level_Chen == TRUE & Level_BIC == TRUE  ~ "metric",
                                    Level_Tested == "configural" & Dimension == "AGRE" ~ "configural",
                                    Level_Tested == "configural" & Dimension == "CONS" ~ "(configural)",
                                    Level_Tested == "configural" & Dimension == "EMOS" ~ "(configural)",
                                    TRUE ~ "X"))  %>% 
  mutate(N = as.numeric(N), ChiSquare = as.numeric(ChiSquare), pvalue = as.numeric(pvalue), CFI = as.numeric(CFI), 
         RMSEA = as.numeric(RMSEA), SRMR = as.numeric(SRMR), BIC = as.integer(BIC)) %>% 
  mutate(across(c("CFI", "RMSEA", "SRMR"), round, 3)) %>% 
  mutate(across(c("ChiSquare"), round, 2))


# Tables for presentation

mi_fit_tab <- fit.table.mi %>% 
  select(dimension, type, chisq.scaled, df.scaled, pvalue.scaled, cfi.scaled, rmsea.scaled, 
         srmr, bic, high.level.chen, high.level.bic) %>% 
  rename(Dimension = dimension, Group = group, Model = type, ChiSq = chisq.scaled, 
         df = df.scaled, pvalue = pvalue.scaled, CFI = cfi.scaled, RMSEA = rmsea.scaled, 
         SRMR = srmr, BIC = bic, MI_Chen = high.level.chen, 
         MI_BIC = high.level.bic) %>% 
#  mutate(across(where(is.numeric), formatter)) %>% 
  relocate(Dimension, .before = "Group") %>% 
  arrange(Dimension, Group, Model)

# mi_fit_tab  %>%
#   kbl(caption = "Model fit indices of measurement invariance models") %>%
#   kable_classic(full_width = F, html_font = "Cambria") %>% 
#   save_kable("04_reports/mi_fit_tab.pdf")

# Create table for highest measurement invariance level attained
mi.table   <- extract_mi_level(fit.table.mi) %>% 
  mutate(int.dec = flatten_chr(int.dec))

# Tables for presentation

mi_bic = mi_fit_tab %>% 
  filter(MI_BIC != "-") %>% 
  rename("MI BIC" = MI_BIC)

mi_qual <- mi.table %>% 
  select(Dimension = dimension,  Group = group, "MI Chen" = int.dec) %>% 
  arrange(Dimension, Group) %>% 
  left_join(., 
            mi_bic[c("Dimension", "Group", "MI BIC")], 
            by = c("Dimension", "Group"))




# Run analyses with partial scalar invariance -----------------------------


partial_scalar_models <- mi_results %>% 
  filter(type == "scalar") %>% 
  mutate(free_intercept = map_chr(mod_output,
                                  ~modindices(.x, 
                                              sort = TRUE,
                                              free.remove = FALSE) %>%
                                    filter(op == "~1") %>% 
                                    pluck("lhs") %>% 
                                    first())) %>% 
  
  rowwise() %>% 
  mutate(form =  if_else(group == "age_group",
                         str_c(form, str_c(";",
                                           free_intercept, 
                                           "~c(a, b, c)*1")),
                         str_c(form, str_c(";",
                                           free_intercept, 
                                           "~c(a, b)*1")))) %>% 
  select(-mod_output, -mod_results) %>% 
  ungroup()

# Add second free intercept for some models where modindices were very high

partial_scalar_models$form[[12]] <- str_c(partial_scalar_models$form[[12]],
                                          "; bdze019a~c(c, d)*1")
partial_scalar_models$form[[11]] <- str_c(partial_scalar_models$form[[11]],
                                          "; bdze015a~c(c, d)*1")
partial_scalar_models$form[[9]] <- str_c(partial_scalar_models$form[[9]],
                                          "; bdze024a~c(c, d)*1")

#partial_scalar_models$form = list(partial_scalar_models$form)

lookup_mi_partial <- lookup_mi %>% 
  filter(!type == "scalar") %>% 
  unnest(form) %>% 
  bind_rows(partial_scalar_models) %>% 
  select(-free_intercept)


mi_results_partial <- lookup_mi_partial %>% 
  mutate(mod_output = future_pmap(list(items, type, form, group), 
                                  mi_mod),
         mod_results = future_map(mod_output, possibly(mi_evaluate, otherwise = NA)) 
  )



# Summarize the partial scalar MI results ----------------------------------------

# Determine measurement invariance (Multi-Group CFA)
fit.table.mi.partial <- mi_results_partial %>% 
  #group_by(dimension, group) %>% 
  extract_fit() %>% 
  extract_delta_fit()  %>% 
  extract_invariance()

traditionalMI_tbl_partial <- fit.table.mi.partial %>% select(dimension, type, ntotal, chisq.scaled, 
                                             pvalue.scaled, cfi.scaled, rmsea.scaled, 
                                             srmr, bic, level.chen, level.bic) %>%  
  rename(Dimension = dimension, Level_Tested = type, N = ntotal, ChiSquare = chisq.scaled, 
         pvalue = pvalue.scaled, CFI = cfi.scaled, RMSEA = rmsea.scaled, SRMR = srmr, 
         BIC = bic, Level_Chen = level.chen, Level_BIC = level.bic)  %>%  
  mutate(Level_Chen = case_when(Level_Tested == "configural" ~ "", TRUE ~ Level_Chen),
         Level_BIC = case_when(Level_Tested == "configural" ~ "", TRUE ~ Level_BIC),
         Level_achieved = case_when(Level_Tested == "scalar" & Level_Chen == TRUE & Level_BIC == TRUE  ~ "scalar",
                                    Level_Tested == "metric" & Level_Chen == TRUE & Level_BIC == TRUE  ~ "metric",
                                    Level_Tested == "configural" & Dimension == "AGRE" ~ "configural",
                                    Level_Tested == "configural" & Dimension == "CONS" ~ "(configural)",
                                    Level_Tested == "configural" & Dimension == "EMOS" ~ "(configural)",
                                    TRUE ~ "X"))  %>% 
  mutate(N = as.numeric(N), ChiSquare = as.numeric(ChiSquare), pvalue = as.numeric(pvalue), CFI = as.numeric(CFI), 
         RMSEA = as.numeric(RMSEA), SRMR = as.numeric(SRMR), BIC = as.integer(BIC)) %>% 
  mutate(across(c("CFI", "RMSEA", "SRMR"), round, 3)) %>% 
  mutate(across(c("ChiSquare"), round, 2))


# Tables for presentation

mi_fit_tab_partial <- fit.table.mi.partial %>% 
  select(dimension, type, chisq.scaled, df.scaled, pvalue.scaled, cfi.scaled, rmsea.scaled, 
         srmr, bic, high.level.chen, high.level.bic) %>% 
  rename(Dimension = dimension, Group = group, Model = type, ChiSq = chisq.scaled, 
         df = df.scaled, pvalue = pvalue.scaled, CFI = cfi.scaled, RMSEA = rmsea.scaled, 
         SRMR = srmr, BIC = bic, MI_Chen = high.level.chen, 
         MI_BIC = high.level.bic) %>% 
  #  mutate(across(where(is.numeric), formatter)) %>% 
  relocate(Dimension, .before = "Group") %>% 
  arrange(Dimension, Group, Model)

# mi_fit_tab  %>%
#   kbl(caption = "Model fit indices of measurement invariance models") %>%
#   kable_classic(full_width = F, html_font = "Cambria") %>% 
#   save_kable("04_reports/mi_fit_tab.pdf")

# Create table for highest measurement invariance level attained
mi.table.partial   <- extract_mi_level(fit.table.mi.partial) %>% 
  mutate(int.dec = flatten_chr(int.dec))

# Tables for presentation

mi_bic_partial = mi_fit_tab_partial %>% 
  filter(MI_BIC != "-") %>% 
  rename("MI BIC" = MI_BIC)

mi_qual_partial <- mi.table.partial %>% 
  select(Dimension = dimension,  Group = group, "MI Chen" = int.dec) %>% 
  arrange(Dimension, Group) %>% 
  left_join(., 
            mi_bic[c("Dimension", "Group", "MI BIC")], 
            by = c("Dimension", "Group"))


# Save results  ------------------------------------------------------

save(list = c("mi_qual", "mi_fit_tab", "mi_qual_partial", "mi_fit_tab_partial"), 
     file = file.path(dirs$results, "mi_results.Rdata")
)
