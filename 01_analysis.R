# Basic settings ----------------------------------------------------------

# Load libraries required for data wrangling

library(tidyverse)
library(psych)
library(glue)
library(sjlabelled)
library(furrr)
library(EFA.dimensions)
library(lavaan)
library(parameters)


# Setup plan for multi-processor workflows with furrr
plan(multisession, workers = 4)
set.seed(123)

# Create a list of directories for shorter code in subsequent tasks
# For replication, replace the one absolute path to the raw data with the one
# in your file system containing the GESIS Panel data (SPSS *.sav datasets).
# All other paths are relative (to the R project path).


source("helper_functions.R")

# Load the data and helper objects created in 00_data_preparation.R
load(file.path(dirs$data, "gpdata.Rdata"))


# Descriptive statistics --------------------------------------------------

generate_descriptives <- function(data, items) {
  
  data %>% 
    #filter(participation_t1 == 1) %>% 
    summarise(across(all_of(items),
                     list(
                       "mean" = mean, 
                       "sd" = sd,
                       "min" = min,"max" = max,
                       "skewness" = moments::skewness,
                       "kurtosis" = moments::kurtosis,
                       "valid" =  ~sum(!is.na(.x)), 
                       "missing" = ~sum(is.na(.x)) / length(.x) * 100
                     ),
                     na.rm = TRUE
    )
    ) %>% pivot_longer(everything(),
                       names_to = c("variable", "stat"),
                       #names_pattern = "([0-9]{1,3}_[a-z]{2,4}[0-9]{1})_(.*)") %>% 
                       names_pattern = "(.*)_(.*)") %>% 
    pivot_wider(id_cols = c("variable"), names_from = "stat")
}

renamer <- scale_info$`Wording (English)`
names(renamer) <- scale_info$`Item Label (GESIS Panel)`


item_descriptives <- generate_descriptives(gpdata, t1_value_items) %>% 
  mutate(wording = str_replace_all(variable, renamer)) %>% 
  left_join(scale_info[c("Item", "Item Label (GESIS Panel)")], 
            by = c("variable" = "Item Label (GESIS Panel)")) %>% 
  relocate(Item, wording, .before = "mean") 

scale_descriptives <- generate_descriptives(gpdata, 
  items =   c(t1_value_dimensions,
              t1_value_dimensions_ips)
)

rm(renamer)
# Scale intercorrelations -------------------------------------------------

dimension_correlations <- gpdata %>% 
  select(t1_value_dimensions) %>% 
  cor(use = "pairwise")

ips_cors <- gpdata %>% 
  select(t1_value_dimensions_ips) %>% 
  cor(use = "pairwise")

dimension_correlations[lower.tri(dimension_correlations)] <- ips_cors[lower.tri(ips_cors)]

rm(ips_cors)
# Dimensionality ----------------------------------------------------------

lookup_table <- lookup_value_dimensions[1:4] %>% 
  enframe(name = "Dimension", value = "Items") 


test_dimensionality <- function(items, data = gpdata) {
  
  dat <- data[items]
  print(head(dat))
  # Write a function that runs all dimensionality analyses
  
  pa_output <- psych::fa.parallel(
    dat,
    plot = FALSE,
    #n.obs = nrow(data),
    fa = "pc") %>%
    suppressMessages()
  
  
  # Extract number of factors and principal components from PA
  pa_results <- pa_output %>%
    pluck("ncomp") %>% 
    as_tibble_col("ncomp") %>% 
    rename("PA" = ncomp)
  
  # Extract eigenvalue information from PA output
  eigen_values <- pa_output %>%
    magrittr::extract("pc.values")
  
  eigen_info <- tibble(
    first_eigen = eigen_values[["pc.values"]][1],
    second_eigen = eigen_values[["pc.values"]][2],
    eigen_ratio = first_eigen / second_eigen,
  ) %>%
    round(2)
  
  # Run MAP using the EFA.dimensions package and extract revised MAP criterion
  map_results <- MAP(
    data = dat, verbose = FALSE
  ) %>%
    chuck("NfactorsMAP4") %>%
    tibble(MAP = .) %>% 
    suppressMessages()
  
  # Rund EKC using EFA.dimensions and
  ekc_results <- EMPKC(
    data = dat,
    #  Ncases = count(data),
  ) %>% 
    chuck("NfactorsEMPKC") %>%
    tibble(EKC = .) %>% 
    suppressMessages()
  
  # Combine all results to a tibble
  dimensionality_results <- bind_cols(
    pa_results,
    map_results,
    ekc_results,
    eigen_info
  ) 
  
  dimensionality_results
}


# Apply the function to value dimensions 

dimensionality_results <- lookup_table %>% 
  bind_rows(
    tibble(Dimension = "Full item set", Items = list(t1_value_items))
  ) %>% 
  mutate(results = future_map(Items, 
                              .f = safely(test_dimensionality),
                              .options = furrr_options(seed = 123))
  )  %>% 
  select(-Items ) %>% 
  unnest_wider(results, simplify = FALSE) %>% 
  unnest_wider(result) # %>% 
  #mutate(across(c(first_eigen, second_eigen, eigen_ratio), formatter))


test_dimensionality <- function(items, data = gpdata) {
  
  dat <- data[items]
  print(head(dat))
  # Write a function that runs all dimensionality analyses
  
  pa_output <- psych::fa.parallel(
    dat,
    plot = FALSE,
    #n.obs = nrow(data),
    fa = "pc") %>%
    suppressMessages()
  
  
  # Extract number of factors and principal components from PA
  pa_results <- pa_output %>%
    pluck("ncomp") %>% 
    as_tibble_col("ncomp") %>% 
    rename("PA" = ncomp)
  
  # Extract eigenvalue information from PA output
  eigen_values <- pa_output %>%
    magrittr::extract("pc.values")
  
  eigen_info <- tibble(
    first_eigen = eigen_values[["pc.values"]][1],
    second_eigen = eigen_values[["pc.values"]][2],
    eigen_ratio = first_eigen / second_eigen,
  ) %>%
    round(2)
  
  # Run MAP using the EFA.dimensions package and extract revised MAP criterion
  map_results <- MAP(
    data = dat, verbose = FALSE
  ) %>%
    chuck("NfactorsMAP4") %>%
    tibble(MAP = .) %>% 
    suppressMessages()
  
  # Rund EKC using EFA.dimensions and
  ekc_results <- EMPKC(
    data = dat,
    #  Ncases = count(data),
  ) %>% 
    chuck("NfactorsEMPKC") %>%
    tibble(EKC = .) %>% 
    suppressMessages()
  
  # Combine all results to a tibble
  dimensionality_results <- bind_cols(
    pa_results,
    map_results,
    ekc_results,
    eigen_info
  ) 
  
  dimensionality_results
}


# Apply the function to value dimensions 

dimensionality_results <- lookup_table %>% 
  bind_rows(
    tibble(Dimension = "Full item set", Items = list(t1_value_items))
  ) %>% 
  mutate(results = future_map(Items, 
                              .f = safely(test_dimensionality),
                              .options = furrr_options(seed = 123))
  )  %>% 
  select(-Items ) %>% 
  unnest_wider(results, simplify = FALSE) %>% 
  unnest_wider(result) # %>% 


# New dimensionality analyses using the EFA.dimensions package only  ----------------


test_dimensionality2 = function(items, data = gpdata) {
  
  dat <- data[items]

  
  dimensionality_results =  EFA.dimensions::DIMTESTS(dat, 
                                                     tests = c("RAWPAR", "EMPKC", "MAP", "VSS"),
                                                     display = 0) %>% 
    pluck("dimtests") %>% 
    as_tibble(rownames = "test") %>% 
    pivot_wider(names_from = 1, values_from = 2)
  
  dimensionality_results
}



dimensionality_results2 =  lookup_table %>% 
  bind_rows(
    tibble(Dimension = "Full item set", Items = list(t1_value_items))
  ) %>% 
  mutate(results = future_map(Items, 
                              .f = test_dimensionality2,
                              .options = furrr_options(seed = 123))
  )  %>% 
  select(-Items ) %>% 
  unnest_wider(results, simplify = TRUE) 

dimensionality_results2

# Reliability -------------------------------------------------------------

reliability_results <- future_map_dfr(lookup_value_dimensions[1:4],
                                      .id = "Dimension",
                                      .f = ~ psych::reliability(corFiml(gpdata[unlist(.x)])) %>% 
                                        pluck("result.df") %>% 
                                        as_tibble() %>%
                                        select(Alpha = alpha, Omega = omega.tot, max.split, min.split)) %>% 
  mutate("Split-Half" = str_c(sprintf(min.split, fmt = "%.2f"),
                              "-", 
                              sprintf(max.split, fmt = "%.2f"))) %>% 
  select(-max.split, -min.split) 

#Test-retest stability of the facets
rtt <- tibble(
  Dimension = names(lookup_value_dimensions)[1:4],
  Dimension_T2 = names(lookup_value_dimensions)[5:8]
) %>% 
  mutate("Test-Retest" = future_map2_dbl(Dimension, Dimension_T2, 
                           .f = ~ gpdata[gpdata$participation_t1 == 1, c(.x, .y)] %>% 
                             corFiml() %>% .[c(2)]))
  


# Add to reliability table
reliability_results <- left_join(reliability_results,
                                 rtt[c("Dimension", "Test-Retest")],
                                 by = "Dimension") %>% 
    mutate(across(where(is.numeric), ~sprintf(.x, fmt = "%#.2f")))


# CFA model fits ----------------------------------------------------------

cfa_form <- function(dimension, items, fixed = TRUE){
  is_reversed <- str_detect(items, "_r")
  fixed_items <- str_c("1*", items)
  
  
  if (fixed == TRUE) {
    form <- str_c(dimension, " =~ ", 
                  str_c(fixed_items, collapse = " + "))
    
  } else {
    
    form <- str_c(dimension, " =~ ", 
                str_c(items, collapse = " + "))
    
  }
  
 form
}

cfa_mod <- function(items, form, data = gpdata) {
  dat <- data %>% 
  #  filter(participation_t1 == 1) %>% 
    select(all_of(items)) 
  
  if(nrow(dat) == 0){
    mod_output <- NA
  }
  else {
    mod_output <- lavaan::cfa(model = form, 
                              data = dat, 
                              missing = 'FIML', 
                              estimator = "MLR")
  }
  mod_output
}

### FUNCTIONS cfa_converged(), cfa_iterations(), cfa_fit() ###
##  for extracting info from fitted models


cfa_evaluate <- function(mod_output) {
  
  
  stopifnot()
  if (!is(mod_output, "lavaan")) {
    
    stop("No valid input data or model input")
    
  }
  
  cfa_evaluation <- list(
    cfa_converged = lavInspect(mod_output, "converged"),
    cfa_iterations = lavInspect(mod_output, "iterations"),
    cfa_fit = fitmeasures(mod_output) %>% 
      as_tibble_row() ,
    cfa_param = model_parameters(mod_output)
  )
  
  cfa_evaluation
}

# Create tibble for models
cfa_results <- lookup_table |>
  mutate(#tau_equivalent  = future_map2(str_sub(Dimension, 1, 4), Items, cfa_form, fixed = TRUE),
         form = future_map2(str_sub(Dimension, 1, 4), Items, cfa_form, fixed = FALSE),
         Model = "congeneric") 
cfa_results[5,] <- cfa_results[1,]
cfa_results$form[[5]] <- str_c(cfa_results$form[[1]], ";bdze019a~~bdze023a")
cfa_results$Model[[5]] <- "congeneric residual correlation"

cfa_results[6, ] <- cfa_results[4, ]
cfa_results$form[[6]] <- "cons =~ 1*bdze014a + bdze017a + 1*bdze022a"
cfa_results$Model[[6]] <- "1 item fixed"
cfa_results <- cfa_results |> 
  mutate(mod_output = future_map2(Items, form, cfa_mod),
         mod_results = future_map(mod_output, safely(cfa_evaluate)) 
  )


# # CFA additional manual models --------------------------------------------
# 
# ## Add Residual correlation in the Self-Transcendence Model
# map(cfa_results[[2,"mod_output"]], lavaan::modificationindices) 
# cfa_results[9, ] <- cfa_results[2, ]
# cfa_results$form[[9]] <- str_c(cfa_results$form[[9]], ";bdze019a~~bdze023a")
# cfa_results$mod_results[[9]] <- cfa_evaluate(cfa_results$mod_output[[9]])
# cfa_results$mod_output[[9]] <- cfa_mod(cfa_results$Items[[9]], cfa_results$form[[9]])
# cfa_results$Model[[9]] <- str_c(cfa_results$Model[[9]], " (RCV)")
# 
# 
# ## Estimate a model with equal factor loadings for Tradition and Conformity
# map(cfa_results[[2,"mod_output"]], lavaan::modificationindices) 
# cfa_results[10, ] <- cfa_results[8, ]
# cfa_results$form[[10]] <- "cons =~ 1*bdze014a + bdze017a + 1*bdze022a"
# cfa_results$mod_results[[10]] <- cfa_evaluate(cfa_results$mod_output[[10]])
# cfa_results$mod_output[[10]] <- cfa_mod(cfa_results$Items[[10]], cfa_results$form[[10]])
# cfa_results$Model[[10]] <- "partly constrained"


# # Model
# g_factor <- str_c("general =~ ",
#                   str_c("1*", t1_value_items, collapse = " + ") 
# ) %>% 
#   str_c("general ~~ 0*tran; general ~~ 0*open; general ~~ 0*cons;
#         general ~~ 0*enha;", sep = ";")
# 
# #mimic <- str_c(t1_value_items, " ~ wpmean", collapse = "; ")
# 
# joint_model <- cfa_results %>% 
#   filter(Model == "tau_congeneric") %>% 
#   pluck("form") %>% 
#   purrr::reduce(.f = str_c, sep = "; ") %>% 
#   str_c(g_factor, sep = "; ") %>% 
#   cfa(data = gpdata, estimator = "MLR", missing = "FIML", orthogonal = FALSE)


cfa_table <- cfa_results %>% 
  filter(!(Dimension == "conservation_t1" & Model == "congeneric")) |>
 unnest_wider(mod_results) %>% 
  unnest_wider(result) %>% 
  select(Dimension, Model, cfa_fit) %>% 
  unnest(cfa_fit) %>% 
  select(Dimension, Model, 
         "Chi-2" = chisq.scaled, 
         "df" = df.scaled,
         "p value" = pvalue.scaled,
         CFI = cfi.robust, 
         RMSEA = rmsea.robust,
         SRMR = srmr, BIC = bic) %>% 
  mutate(across(where(is.numeric), round, 3))

# Multidimensional scaling (MDS) ------------------------------------------
library(smacof)

# Set the type of MDS to run
mds_type = "interval"

# Compute a distance matrix
dist_mat <- gpdata |>
  mutate(across(all_of(all_value_items[1:17]), ~.x - wpmean)) |>
  select(all_of(all_value_items[1:17])) |>
  psych::corFiml() |>
  sim2diss(to.dist = TRUE) 

# Run th emds

# mds_results <- mds(dist_mat, type = mds_type)
# 
# # Generate a random stress value for comparison for 17 items
# rs_vec <- randomstress(17, ndim = 2, nrep = 500, type = mds_type)
# 
# # Random permutation test (Mair et al., 2016)
# mds_permtest <- permtest(mds_results,
#                          data =  select(gpdata, all_of(all_value_items[1:17])),
#                          method.dat = "euclidean",
#                          nrep = 500, 
#                          verbose = FALSE)


mds_configuration = as_tibble(mds_results$conf, rownames = "Item Label (GESIS Panel)") |>
  left_join(scale_info[c("Item", "Item Label (GESIS Panel)", "Wording (English)", "Higher-order value")],
            by = "Item Label (GESIS Panel)") |>
  rename(c("Dimension 1" = "D1", "Dimension 2" = "D2"))

# Robustness check: Run a PCA with the ipsatized data
gpdata %>% 
  mutate(across(all_of(all_value_items[1:17]), ~.x - wpmean)) %>% 
  select(all_of(all_value_items[1:17])) %>% 
  psych::principal(rotate = "varimax", nfactors = 2) %>% plot()


# EFA / ESEM --------------------------------------------------------------


renamer <- scale_info$Item
names(renamer) <- scale_info$`Item Label (GESIS Panel)`

# EFA
library(GPArotation)

efa_results <- gpdata |>
  #mutate(across(all_of(all_value_items[1:17]), ~.x - wpmean)) |>
  select(all_of(t1_value_items)) |>
  psych::fa(fm = "ml", nfactors = 4, rotate = "none")
  
# Create a target matrix
targets = matrix(0, nrow = 17, ncol = 4)
targets [1:5, 1] <- 1
targets [6:9, 2] <- 1
targets [10:14, 3] <- 1
targets [15:17, 4] <- 1

rotated_solution <- targetQ(efa_results$loadings,
                            Tmat = Random.Start(4),
                            Target = targets,
                            maxit = 10000)

efa_results$loadings <- rotated_solution$loadings
efa_results$phi <- rotated_solution$Phi

domain_loadings <- parameters::model_parameters(efa_results, 
                                                sort = F)  

domain_loadings$Variable = str_replace_all(domain_loadings$Variable, renamer)

print(domain_loadings, threshold = 0.3)


# Test: Run ESEM with target rotation in Mplus

library(MplusAutomation)

mplus_input <- mplusObject(
TITLE =   "Random-intercept EFA for values",
  ANALYSIS = "ROTATION = TARGET;",

  MODEL = "Tran BY bdze011a~1 bdze015a~1 bdze019a~1 bdze023a~1 bdze026a~1 ! Tran
          bdze013a~0 bdze018a~0 bdze020a~0 bdze024a~0 bdze027a~0 !Open
          bdze012a~0 bdze016a~0 bdze021a~0 bdze025a~0 ! Enhancement
         bdze014a~0 bdze017a~0 bdze022a~0 (*1); ! Conservation

Open BY   bdze013a~1 bdze018a~1 bdze020a~1 bdze024a~1 bdze027a~1 !Open
bdze011a~0 bdze015a~0 bdze019a~0 bdze023a~0 bdze026a~0 ! Tran
 bdze012a~0 bdze016a~0 bdze021a~0 bdze025a~0 ! Enhancement
         bdze014a~0 bdze017a~0 bdze022a~0 (*1); ! Conservation

Enha BY  bdze012a~1 bdze016a~1 bdze021a~1 bdze025a~1 ! Enhancement
bdze011a~0 bdze015a~0 bdze019a~0 bdze023a~0 bdze026a~0 ! Tran
          bdze013a~0 bdze018a~0 bdze020a~0 bdze024a~0 bdze027a~0 !Open
     bdze014a~0 bdze017a~0 bdze022a~0 (*1); ! Conservation

Cons BY  bdze014a~1 bdze017a~1 bdze022a~1
        bdze011a~0 bdze015a~0 bdze019a~0 bdze023a~0 bdze026a~0 ! Tran
          bdze013a~0 bdze018a~0 bdze020a~0 bdze024a~0 bdze027a~0 !Open
          bdze012a~0 bdze016a~0 bdze021a~0 bdze025a~0 ! Enhancement
         (*1); !Conservation

g by bdze011a-bdze022a@1;
g WITH tran-cons@0
",
OUTPUT = "stdyx;"
, rdata = gpdata[t1_value_items], autov = TRUE
)

mplus_output <- mplusModeler(mplus_input, 
                             dataout = "01_data/gpdata.dat",
                             run = 1)


# Test: Run ESEM in Lavaan
inverted_renamer = names(renamer)
names(inverted_renamer)= unname(renamer) %>% str_replace_all("([0-9]{1,2})([A-Z]{1,2})", "\\2\\1" )

esem_model <- '
    # efa block
    efa("efa1")*f1 + 
    efa("efa1")*f2 +
    efa("efa1")*f3 +
    efa("efa1")*f4  =~ CO12 + CO7 + CO4 + OC14 + OC3 + OC8 + OC17 + OC10 + 
                      SE2 + SE15 + SE11 + SE6 + ST9 + ST13 + ST16 + ST1 + ST5

    # cfa block
 g =~ 1*CO12 + 1*CO7 + 1*CO4 + 1*OC14 + 1*OC3 + 1*OC8 + 1*OC17 + 1*OC10 + 
      1*SE2 + 1*SE15 + 1*SE11 + 1*SE6 + 1*ST9 + 1*ST13 + 1*ST16 + 1*ST1 + 1*ST5
#
#g =~ CO12 + CO7 + CO4 + OC14 + OC3 + OC8 + OC17 + OC10 + 
#                      SE2 + SE15 + SE11 + SE6 + ST9 + ST13 + ST16 + ST1 + ST5
   g ~~ 0*f1
   g ~~ 0*f2
   g ~~ 0*f3
   g ~~ 0*f4
   
   
'

cfa_model <- '

CO =~ 1*CO12 + start(1)*CO7 + start(1)*CO4
OC =~ 1*OC14 + start(1)*OC3 + start(1)*OC8 + start(1)*OC17 + start(1)*OC10  
SE =~ 1*SE2 + start(1)*SE15 + start(1)*SE11 + start(1)*SE6 
ST =~ 1*ST9 + start(1)*ST13 + start(1)*ST16 + start(1)*ST1 + start(1)*ST5

g =~ 1*CO12 + 1*CO7 + 1*CO4 + 1*OC14 + 1*OC3 + 1*OC8 + 1*OC17 + 1*OC10 + 
                      1*SE2 + 1*SE15 + 1*SE11 + 1*SE6 + 1*ST9 + 1*ST13 + 1*ST16 + 1*ST1 + 1*ST5
   g ~~ 0*CO
   g ~~ 0*OC
   g ~~ 0*ST
   g ~~ 0*SE
   
   
'


fit <- sem(model = esem_model, 
           data = gpdata[t1_value_items] %>% 
             rename(all_of(inverted_renamer)),
           rotation = "target",
           rotation.args = list(target = targets))

fit <- efa(data = gpdata[t1_value_items] %>% 
             rename(all_of(inverted_renamer)),
           nfactors = 4,
           rotation = "target",
           rotation.args = list(target = targets))



fit <- cfa(model = cfa_model, 
           data = gpdata[t1_value_items] %>% 
             rename(all_of(inverted_renamer)))

summary(fit, fit = TRUE, standardized = TRUE)

# Criterion correlations --------------------------------------------------

# Define a vector

library(corrr)

# Vector of correlates/criteria (for selecting and sorting)
# all_criteria <- c(
#   "Age (years)",
#   "Gender",
#   "Education",
#   "Income",
#   "Work",
#   "Finances",
#   "Leisure",
#   "Family",
#   "Friends",
#   "Neighbors",
#   "Health",
#   "Open-Mindedness",
#   "Conscientiousness",
#   "Extraversion",
#   "Agreeableness",
#   "Emotional Stability",
#   "Left–Right Placement",
#   "Life Satisfaction",
#   "Happiness"
#   )


all_criteria = c(
  age = "Age (years)",
  gender = "Gender",
  edu_group = "Education",
  "income" = "Personal net income",
  
  "bfi_open" = "Open-Mindedness",
  "bfi_con" = "Conscientiousness",
  "bfi_ext" = "Extraversion",
  "bfi_agg" = "Agreeableness",
  "bfi_emo" = "Emotional Stability",
  
  importance_work = "Work", # Importance of work
  importance_finance = "Finances", # Importance of financial situation
  importance_leisure = "Leisure", # Importance of leisure
  importance_family = "Family", # Importance of own family
  importance_friends = "Friends", # Importance of friends
  importance_neighbors = "Neighbors", # Importance of neighbors
  importance_health = "Health",
  
  "left_right" = "Left–Right Placement",
  
  "satisfaction" = "Life Satisfaction",
  happiness = "Happiness"
  
)

# Compute the within-person mean for the importance of life domain items 
# such that these items can be centered in the next step
gpdata$wpmean_importance <- gpdata %>% 
  select(starts_with("importance_")) %>% 
  rowMeans(na.rm = FALSE)
  

# Compute correlation with outcomes
# using centered values and importance items
criterion_correlations <- gpdata %>% 
  mutate(across(starts_with("importance"), ~ .x - wpmean_importance)) %>% 
  select(t1_value_dimensions_ips,
         age,
         gender,
         edu_group,
         starts_with("importance"),
         starts_with("bfi"),
         satisfaction,
         happiness,
         left_right,
         income,
         ) %>% 
  correlate() %>% 
  filter(!term %in% names(t1_value_dimensions_ips)) %>% 
  select(term, names(t1_value_dimensions_ips)) %>% 
  #arrange(term, match(c("age", "gender", "edu_group", "importance"))) %>% 
  mutate(term = str_replace_all(term,
                                 all_criteria)) %>% 
  arrange(match(term, all_criteria)) %>% 
  rename_with(.fn = ~str_replace_all(.x, "\\ \\(ctrd.\\)", ""))
  
  adjusted_r_squared <-   crossing(outcome = names(all_criteria), 
                          predictors =  list(t1_value_items, t1_value_dimensions)) %>% 
    mutate(fml = map2(outcome, predictors, ~glue("{.x} ~ {str_c(.y, collapse = '+')}")),
           output = map(fml, ~lm(.x, data = gpdata)),
           tidy_output = map(output, broom::glance)) %>% 
    unnest_wider(tidy_output) %>% 
    mutate(type = map_dbl(predictors, length), preds =  str_replace_all(type, c("4" = "Higher-order values", "17" = "Single items"))) %>% 
    select(outcome, preds, adj.r.squared) %>% 
    mutate(outcome = str_replace_all(outcome, all_criteria)) %>% 
   pivot_wider(id_cols = outcome, 
               names_from = preds, 
               values_from = adj.r.squared) %>% 
    arrange(match(outcome, all_criteria))
  
# Save the results --------------------------------------------------------


save(list = c(
  "all_criteria",
  "item_descriptives",
  "scale_descriptives",
  "mds_results",
  "mds_configuration",
  "dimension_correlations",
  "dimensionality_results",
  "reliability_results",
  "cfa_results",
  "cfa_table",
  "efa_results",
  "domain_loadings",
  "criterion_correlations",
  "adjusted_r_squared"
),
     file = file.path(dirs$results, "all_results.Rdata"))
  
  
  
# ADDITIONAL ANALYSES FOR REVISION

# Write correlatoin table with centered scores for revision ---------------------------------------------

library(apaTables)
criterion_names = names(all_criteria)
names(criterion_names) = unname(all_criteria)

efa_factor_cors = efa_results$phi 
colnames(efa_factor_cors) = c("ST", "SE", "OC", "CO")
rownames(efa_factor_cors) = c("ST", "SE", "OC", "CO")


apa.cor.table(
data =  gpdata %>% 
        mutate(across(starts_with("importance"), ~ .x - wpmean_importance)) %>% 
        select(all_of(c(t1_value_dimensions_ips,
               names(all_criteria)))
    ) %>% 
    rename(all_of(criterion_names)),
  filename="Table1_APA.doc", 
  show.conf.interval = TRUE)

# Write correlatoin table with uncerntered scores for revision ----------------------------

apa.cor.table(
  data =  gpdata %>% 
    select(all_of(c(t1_value_dimensions,
                    names(all_criteria)))
    ) %>% 
    rename(all_of(criterion_names)),
  filename="Table2_APA.doc", 
  show.conf.interval = TRUE)



pca_results <- gpdata |>
  #mutate(across(all_of(all_value_items[1:17]), ~.x - wpmean)) |>
#  mutate(across(all_of(t1_value_items), ~ .x - wpmean_importance)) %>% 
  select(all_of(c(t1_value_items))) |>
  psych::principal(nfactors = 4, rotate = "none")


# Create a target matrix
targets = matrix(0, nrow = 17, ncol = 4)
targets [1:5, 1] <- 1
targets [6:9, 2] <- 1
targets [10:14, 3] <- 1
targets [15:17, 4] <- 1

rotated_solution <- targetQ(pca_results$loadings,
                            Tmat = Random.Start(4),
                            Target = targets,
                            maxit = 10000)

pca_results$loadings <- rotated_solution$loadings
pca_results$phi <- rotated_solution$Phi

domain_loadings <- parameters::model_parameters(pca_results, 
                                                sort = F)  

domain_loadings$Variable = str_replace_all(domain_loadings$Variable, renamer)

print(domain_loadings, threshold = 0.3)
