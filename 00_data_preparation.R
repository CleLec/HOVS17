# Basic settings ----------------------------------------------------------

# Load libraries required for data wrangling

library(tidyverse)
library(psych)
library(glue)
library(sjlabelled)
library(sjmisc)

source("helper_functions.R")

# Create helper objects: Lookup tables, vectors, and lists  ----------------------
# These helper objects will be used for selecting and iterating over variables.

# Vector with all value variables at T1 and T2, respectively 

t1_value_items<- c("bdze011a", "bdze015a", "bdze019a", "bdze023a", "bdze026a",
                    "bdze012a", "bdze016a", "bdze021a", "bdze025a",
                    "bdze013a", "bdze018a", "bdze020a", "bdze024a", "bdze027a",
                    "bdze014a", "bdze017a", "bdze022a")

all_value_items <- c("bdze011a", "bdze015a", "bdze019a", "bdze023a", "bdze026a",
                  "bdze012a", "bdze016a", "bdze021a", "bdze025a",
                  "bdze013a", "bdze018a", "bdze020a", "bdze024a", "bdze027a",
                  "bdze014a", "bdze017a", "bdze022a",
                  "cdze011a", "cdze015a", "cdze019a", "cdze023a", "cdze026a",
                  "cdze012a", "cdze016a", "cdze021a", "cdze025a",
                  "cdze013a", "cdze018a", "cdze020a", "cdze024a", "cdze027a",
                  "cdze014a", "cdze017a", "cdze022a")


t1_value_dimensions <- c("Openness" = "openness_t1",
                          "Conservation" = "conservation_t1",
                          "Self-Transcendence" = "transcendence_t1",
                          "Self-Enhancement" = "enhancement_t1")
                          

t1_value_dimensions_ips <- c("Openness (ctrd.)" = "openness_t1_ips",
                         "Conservation (ctrd.)" = "conservation_t1_ips",
                         "Self-Transcendence (ctrd.)" = "transcendence_t1_ips",
                         "Self-Enhancement (ctrd.)" = "enhancement_t1_ips")

# List of variables that will serve as correlates
lookup_correlates <- c(
  "bazb007a", # Importance of own family
  "bazb008a", # Importance of work
  "bazb009a", # Importance of leisure
  "bazb010a", # Importance of friends
  "bazb011a", # Importance of neighbors
  "bazb012a", # Importance of financial situation
  "bazb013a", # Importance of health
  "bcaj055a", # left-right orientation
  #"a12c010a" # left-right orientation
  "bbak098a", # Life satisfaction 
  "bdan123a" # Happiness - > noch korrigieren
  )

# List of Schwartz value variables
lookup_value_dimensions <- list(
  transcendence_t1 = c("bdze011a", "bdze015a", "bdze019a", "bdze023a", "bdze026a"),
  enhancement_t1 = c("bdze012a", "bdze016a", "bdze021a", "bdze025a"),
  openness_t1 = c("bdze013a", "bdze018a", "bdze020a", "bdze024a", "bdze027a"),
  conservation_t1 = c("bdze014a", "bdze017a", "bdze022a"),
  transcendence_t2 = c("cdze011a", "cdze015a", "cdze019a", "cdze023a", "cdze026a"),
  enhancement_t2 = c("cdze012a", "cdze016a", "cdze021a", "cdze025a"),
  openness_t2 = c("cdze013a", "cdze018a", "cdze020a", "cdze024a", "cdze027a"),
  conservation_t2 = c("cdze014a", "cdze017a", "cdze022a")
)

scale_info <- readxl::read_excel(
  file.path(dirs$material, "scale_info.xlsx")
  ) %>% 
  arrange(`Higher-order value`, `Basic human value`)

# Read the data from SPSS -------------------------------------------------

# Data are in separate files per cohort and wave.
# A function is required to read the separate datasets
read_gpdata <- function(file, vars = NULL) {
  path <- file.path(dirs$raw_data, file)
  data <- sjlabelled::read_spss(path, convert.factors	= FALSE)
  
  if(!is.null(vars))
    data <- data |>
          select(all_of(vars))
  
  data
}

# Data from different time points are read from SPSS using the function.
# They are then joined to a single dataframe (gpdata) 
# Note: 
# Wave 1 is the first measurement of values in Cohort A (henceforth T1)
# Wave 2 is the second (henceforth T2)

gpdata <- list(
  basic = read_gpdata("ZA5665_a1_a11-a12_v41-0-0.sav",
                      vars = c("z000001a", # Person ID
                               "a11d054a", # Sex
                               "a11d056b", # Birth year
                               "a11d082b") # Educational attaiment
                                ), 
  wave_a = read_gpdata("ZA5665_a1_aa-ac_v41-0-0.sav"),
  wave_b = read_gpdata("ZA5665_a1_ba-bf_v41-0-0.sav"),
  wave_c = read_gpdata("ZA5665_a1_ca-cf_v41-0-0.sav")
) %>%
  purrr::reduce(left_join, by = "z000001a")


# Compute scores and other variables --------------------------------------


# Compute and recode relevant variables and select only variables needed
# for the analysis
labeled_data <- gpdata

gpdata <- gpdata %>%
  mutate(participation_t1 = ifelse(bdza003a == 1, 1, 0),
         participation_t2 = ifelse(cdza003a == 1, 1, 0),
         gender = ifelse(a11d054a == 1, 1, 0),
         income = ifelse(bfzh088b %in% c(1:15), bfzh088b, NA),
         ID = z000001a,
         education = a11d082b, # a11d082a ohne offene Angaben
         age = ifelse(a11d056b == -99, NA, 2014 - a11d056b),
         across(all_of(all_value_items), 
                ~ifelse(!.x %in% c(1:6), NA, .x)),
         across(all_of(lookup_correlates),
                ~ifelse(.x %in% c(-111, -99, -77, -33, 98 ), NA, .x))
         ) %>% 
  filter(participation_t1 == 1) %>% 
  mutate(edu_group = sjmisc::rec(education, rec = "2:6=0; 7:8=1; else = NA"),
         age_group = case_when(
           age >= 18 & age  < 40 ~ "18-39",
           age >= 40 & age < 60 ~ "40-59",
           age >= 60 ~ "60-70"
         )
) %>% 
  rename(importance_family =  "bazb007a",
         
    importance_family = "bazb007a", # Importance of own family
    importance_work = "bazb008a", # Importance of work
    importance_leisure = "bazb009a", # Importance of leisure
    importance_friends = "bazb010a", # Importance of friends
    importance_neighbors = "bazb011a", # Importance of neighbors
    importance_finance = "bazb012a", # Importance of financial situation
    importance_health = "bazb013a", # Importance of health
    left_right = "bcaj055a", # left-right orientation
    #"a12c010a" # left-right orientation
    satisfaction = "bbak098a", # Life satisfaction 
    happiness = "bdan123a" # Happiness - > noch korrigieren
  )
  


# Compute scores for the values
gpdata <- scale_builder(gpdata, 
                      score_fun = "mean",
                      item_list = lookup_value_dimensions, 
                      na.rm = FALSE) 

# Compute ipsatiszed scores for the values
gpdata$wpmean <- rowMeans(gpdata[t1_value_items], 
                          na.rm = FALSE)

gpdata <- gpdata |>
  mutate(across(c("transcendence_t1", "enhancement_t1",
                  "openness_t1", "conservation_t1"),
                ~.x - wpmean, .names = "{.col}_ips"))

# Compute Big Five scores

lookup_bfi10 <- list(
  bfi_emo = c("bdze004a", "bdze009a_r"),
  bfi_ext = c("bdze006a", "bdze001a_r"),
  bfi_open = c("bdze010a", "bdze005a_r"),
  bfi_con = c("bdze008a", "bdze003a_r"),
  bfi_agg = c("bdze002a", "bdze007a_r")
  )

gpdata <- gpdata |>
  mutate(
   across(c("bdze009a", "bdze001a",
            "bdze005a", "bdze003a",
            "bdze007a", "bdze004a",
            "bdze006a", "bdze010a",
            "bdze008a", "bdze002a"), ~ifelse(.x %in% c(1:5), .x, NA)),
   across(c("bdze009a", "bdze001a",
                  "bdze005a", "bdze003a",
                   "bdze007a"), 
                ~ 6 - .x, .names = "{.col}_r")
 ) %>% 
 scale_builder(score_fun = "mean",
                item_list = lookup_bfi10, 
                na.rm = FALSE) 





# Restore labels
gpdata <- copy_labels(gpdata, labeled_data)


# select(ID, age, gender, education,
#        participation_t1, participation_t2,
#        all_value_items) %>% 
# Save data and helper objects --------------------------------------------

save(file = file.path(dirs$data, "gpdata.Rdata"),
     list = c("gpdata",
              "scale_info",
              "all_value_items",
              "t1_value_items",
              "t1_value_dimensions",
              "t1_value_dimensions_ips",
              "lookup_correlates",
              "lookup_value_dimensions"))
