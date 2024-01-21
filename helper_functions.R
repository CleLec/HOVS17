
# List of directions ------------------------------------------------------

# This list of subdirectories allows for shorter code

# Create a list of directories for shorter code in subsequent tasks
# All paths are relative to the R project
# For replication, copy the SPSS data to the folder in your file system
# in your file system containing the GESIS Panel data (SPSS *.sav datasets).

dirs = list(
  raw_data = "01_data/spss", # Copy the SPSS data into this folder
  material = "00_material",
  data = "01_data",
  results = "02_results"
)


# Data wrangling functions --------------------------------------------

# Write a function that computes scale scores
# and that takes a lookup table with items and scale labels as input
scale_builder <- function(data, item_list, score_fun, suffix = NULL, ...) {
  if (!is.list(item_list) | is.null(names(item_list))) {
    stop("Supply a named list in the item_list argument")
  }
  
  scale_names <- names(item_list)
  score_fun_name <- deparse(substitute(score_fun))
  
  for (i in 1:length(item_list)) {
    
    #new_name <- str_c(scale_names[i], "_", score_fun_name)
    
    new_name <- scale_names[i]
    
    if (!is.null(suffix)) {
      new_name <- str_c(new_name, suffix)
    }
    
    
    if (score_fun == "mean") {
      
      data[[new_name]] <- rowMeans(data[item_list[[i]]], ...)
      
    }
    
    
    if (score_fun == "sum") {
      
      data[[new_name]] <- rowSums(data[item_list[[i]]], ...)
      
    }
    
    message(
      "Computing ", new_name, " as the ", score_fun_name, " score across ",
      length(item_list[[i]]), " items (", 
      str_flatten(item_list[[i]], collapse = ", "), ")."
    )
  }
  data
}


# Data analysis functions -------------------------------------------------


