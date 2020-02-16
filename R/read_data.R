#' read_data
#'
#' This is the primary fuction for reading CellProfiler data into R with this package.
#' 
#' @param filedir The directory with CellProfiler model data to be read. This directory should also contain metadata file. A design file is optional.
#' @param design Logical parameter, if TRUE then a desing file will be joined to data. If FALSE no design file will be joined.
#' @return A single data frame named raw_data that contains all CellProfiler model output and experimental treatments if a design file is used.
#' @export

read_data <- function(filedir, design = FALSE){
  require(tidyverse);require(readr)
  
  file_list <- list.files(filedir,full.names = T) #create a list of files (model output(s)) from your target directory
  model_names <- list.files(filedir) #list of names for each output file
  
  raw_data <- purrr::map2_dfr(file_list, model_names, ~read_rds(.x) %>% ### !! this is where to change to read_csv !! ###
                                dplyr::mutate(model = .y)) %>% #bind data from each output file
    dplyr::mutate(worm_length_um = 3.2937*Worm_Length)
  
  metadata_file <- readr::read_csv(list.files()[grep("metadata", list.files())])
  
  if(!design) { #if you are not using a design file
    print("NO DESIGN FILE LOADED")
    raw_data <- dplyr::inner_join(raw_data, metadata_file) %>%
      tidyr::separate(model, into="model", sep="[:punct:]", extra = "drop")
  } else { #if you are using a design file
    design_file <- readr::read_csv(list.files()[grep("design", list.files())]) #load design file
    raw_data <- dplyr::inner_join(raw_data, metadata_file) %>%
      dplyr::inner_join(., design_file) %>%
      tidyr::separate(model, into="model", sep="[:punct:]", extra = "drop")
  }
  
  return(raw_data)
}
