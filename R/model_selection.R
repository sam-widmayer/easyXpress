#' model_selection
#'
#' This function will assign the appropraite CellProfiler model to each primary object in a raw_data object output by the read_data function.
#' 
#' @param raw_data The directory with CellProfiler model data to be read. This directory should also contain metadata file. A design file is optional.
#' @param model_num Logical parameter, if TRUE then a desing file will be joined to data. If FALSE no design file will be joined.
#' @return A single data frame named model_selected that contains the best CellProfiler model for detecting worm objects within each primary object detected by CellProfiler.
#' @export

model_selection <- function(raw_data, model_num = 4) {
  require(tidyverse); require(readr)
  
  #loads appropriate file based on number of models selected. 
  #Might be best to keep combo.csv files in one place -- can then specify path
  if (model_num == 2) {
    model_selected <- readr::read_csv("/Users/grad/Documents/GitHub/easyXpress/model_combinations/two_model_combo.csv")
    print("SELECTED TWO MODEL FILE")
  } else if (model_num == 3) {
    model_selected <- readr::read_csv("/Users/grad/Documents/GitHub/easyXpress/model_combinations/three_model_combo.csv")
    print("SELECTED THREE MODEL FILE")
  } else if (model_num == 4) {
    model_selected <- readr::read_csv("/Users/grad/Documents/GitHub/easyXpress/model_combinations/four_model_combo.csv")
    print("SELECTED FOUR MODEL FILE")
  }
  
  #join combination file with raw data
  model_selected <- raw_data %>%
    dplyr::group_by(Metadata_Plate, Metadata_Well, Parent_WormObjects, model) %>%
    dplyr::mutate(num_worms = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Metadata_Plate, Metadata_Well, Parent_WormObjects) %>%
    dplyr::distinct(model, .keep_all = T) %>%   
    dplyr::ungroup() %>%
    dplyr::select(Metadata_Plate, Metadata_Well, Parent_WormObjects, model, num_worms) %>%
    tidyr::spread(model, num_worms) %>%
    dplyr::mutate_at(tail(names(.),model_num), ~replace(.,.>2,2)) %>% ### replace any value >2 with 2
    dplyr::mutate_if(is.numeric , replace_na, replace = 0) %>% ### replaced NA with zero
    dplyr::left_join(., model_selected) %>%
    dplyr::full_join(., raw_data) %>%
    dplyr::filter(model == model_select)
  
  return(model_selected)
}