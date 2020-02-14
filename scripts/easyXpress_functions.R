
### easyXpress functions ####

## be sure to set working directory before beginning (project folder)
setwd("/Users/grad/Documents/GitHub/easyXpress/projects/growth/20191119_growth_HB101")

###############################
#### 1. Load Data and Bind ####
###############################

dirs <- "output_100w_model" #folder containing output data from cell profiler pipeline run.
# it is also important to be sure that your output.csv files are named as your model
# ie. mine are named: L1.RDS, L2L3.RDS, L4.RDS, Adult.RDS

## !!note: I have my files saved as .RDS. One may want to change "read_rds" function to "read_csv" if yours are in .csv format!! ##
## option to forgo design file
read_data <- function(filedir, design = FALSE){
  require(tidyverse); require(readr)
  
  file_list <- list.files(filedir,full.names = T) #create a list of files (model output) from your target directory
  model_names <- list.files(filedir) #list of names for each output file
  
  raw_data <- purrr::map2_dfr(file_list, model_names, ~read_rds(.x) %>% ### !! this is where to change to read_csv !! ###
                         dplyr::mutate(model = .y)) %>% #bind data from each output file
    dplyr::mutate(worm_length_um = 3.2937*Worm_Length) %>%
    dplyr::select(-starts_with("Distance"), -starts_with("Intensity"), -starts_with("Location"), -starts_with("Worm_Angle"), -starts_with("Worm_ControlPoint"))
  
  metadata_file <- readr::read_csv(list.files()[grep("metadata", list.files())]) %>% #load metadata file from working directory
    dplyr::select(Image_FileName_RawBF, Metadata_Plate, Metadata_Well)
  
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

#run
raw <- read_data(dirs)

##################################################
##### 2. Model Selection and Flag Assignment #####
##################################################

## option to select number of models
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

#run
mselect <- model_selection(raw, model_num = 4)

#############################
##### 3. Edge Flag Data #####
#############################

## option to set your own radius
edge_flag <- function(raw_data, radius = 825) {
  edge_flag_center_x <- 1024; edge_flag_center_y <- 1024; edge_flag_radius <- radius #set well edge flag parameters for plotting and flagging
  
  edge_flagged <- raw_data %>%
    dplyr::mutate(well_edge_flag = ifelse(sqrt((AreaShape_Center_X - edge_flag_center_x)^2 + (AreaShape_Center_Y - edge_flag_center_y)^2) <= edge_flag_radius, FALSE, TRUE))
  
  return(edge_flagged)
}

#run
edge_flg <- edge_flag(mselect, radius = 825)

########################
##### 4. Flag Data #####    # still need to add bimodality identification
########################

## this function will calculate outliers based on flags of interest (cluster_flag/well_edge_flag)
## it will then record which flags were removed to calculate well outlier (flags_removed column)
set_flags <- function(data, cluster_flag = TRUE, well_edge_flag = TRUE) {
  require(tidyverse)
  
  # removes outliers function (tukey's fences)
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }
  
  if(!cluster_flag & !well_edge_flag) {
    print("NO FLAGS SELECTED FOR FILTERING")
    flag_data <- data %>%
      dplyr::group_by(Metadata_Plate, Metadata_Well) %>%
      dplyr::mutate(well_outlier_flag = (remove_outliers(worm_length_um)),
                    well_outlier_flag = ifelse(is.na(well_outlier_flag), TRUE, FALSE)) %>%
      dplyr::ungroup() %>%
      dplyr::full_join(.,data)
    
  } else if(!well_edge_flag) {
    print("FILTERING CLUSTER FLAGS ONLY")
    flag_data <- data %>%
      dplyr::filter(cluster_flag != T, well_edge_flag !=F) %>% 
      dplyr::mutate(flag_removed = TRUE) %>% # record the flags that were removed
      dplyr::group_by(Metadata_Plate, Metadata_Well) %>%
      dplyr::mutate(well_outlier_flag = (remove_outliers(worm_length_um)),
                    well_outlier_flag = ifelse(is.na(well_outlier_flag), TRUE, FALSE)) %>%
      dplyr::ungroup() %>%
      dplyr::full_join(.,data) %>%
      dplyr::mutate(flags_removed= "rm_cluster_flag") # record the flags that were removed
    
  } else if(!cluster_flag) {
    print("FILTERING WELL EDGE FLAGS ONLY")
    flag_data <- data %>%
      dplyr::filter(cluster_flag != F, well_edge_flag !=T) %>% 
      dplyr::mutate(flag_removed = TRUE) %>%
      dplyr::group_by(Metadata_Plate, Metadata_Well) %>%
      dplyr::mutate(well_outlier_flag = (remove_outliers(worm_length_um)),
                    well_outlier_flag = ifelse(is.na(well_outlier_flag), TRUE, FALSE)) %>%
      dplyr::ungroup() %>%
      dplyr::full_join(.,data) %>%
      dplyr::mutate(flags_removed= "rm_well_edge_flag") 
    
  } else {
    print("FILTERING BOTH CLUSTER AND WELL EDGE FLAGS")
    flag_data <- data %>%
      dplyr::filter(cluster_flag != T, well_edge_flag !=T) %>% 
      dplyr::mutate(flag_removed = TRUE) %>%
      dplyr::group_by(Metadata_Plate, Metadata_Well) %>%
      dplyr::mutate(well_outlier_flag = (remove_outliers(worm_length_um)),
                    well_outlier_flag = ifelse(is.na(well_outlier_flag), TRUE, FALSE)) %>%
      dplyr::ungroup() %>%
      dplyr::full_join(.,data) %>%
      dplyr::mutate(flags_removed= "rm_cluster_flag, rm_well_edge_flag") 
  }
  
  return(flag_data)
}

#run
flag_labeled <- set_flags(edge_flg, cluster_flag = T, well_edge_flag = T)


##########################
##### 5. Filter Data #####
##########################

# this function creates 4 dataframes
## !!! be sure to specify what you are summarizing by in "..." !!! ##
process <- function(flag_data,  ...) {
  require(tidyverse)
  
  # 1 -- data as is, after above steps. Nothing removed
  raw_data <- flag_data
  
  # 2 -- removing rows where flag = TRUE
  processed_data <- flag_data %>%
    dplyr::filter(flag_removed == TRUE) %>% #indicates rows in which flagged data are appropriately removed. See "set_flags()" function
    dplyr::filter(well_outlier_flag == FALSE) #data in this row is NOT an outlier
  
  # 3 -- summarizing #1 by "..."
  summarized_raw <- raw_data %>%
    dplyr::group_by(...) %>%
    dplyr::summarize(mean_wormlength_um = mean(worm_length_um, na.rm = TRUE),
                     min_wormlength_um = as.numeric(quantile(worm_length_um, na.rm = TRUE)[1]),
                     q10_wormlength_um = as.numeric(quantile(worm_length_um, probs = 0.1, na.rm = TRUE)[1]),
                     q25_wormlength_um = as.numeric(quantile(worm_length_um, probs = 0.25, na.rm = TRUE)[1]),
                     median_wormlength_um = median(worm_length_um, na.rm = T),
                     sd_wormlength_um = sd(worm_length_um, na.rm = T),
                     q75_wormlength_um = as.numeric(quantile(worm_length_um, probs = 0.75, na.rm = TRUE)[1]),
                     q90_wormlength_um = as.numeric(quantile(worm_length_um, probs = 0.90, na.rm = TRUE)[1]),
                     max_wormlength_um = as.numeric(quantile(worm_length_um, na.rm = TRUE)[5]),
                     cv_wormlength_um = (sd_wormlength_um/mean_wormlength_um))
  
  # 4 -- summarizing #2 by "..."
  summarized_processed <- processed_data %>%
    dplyr::group_by(...) %>%
    dplyr::summarize(mean_wormlength_um = mean(worm_length_um, na.rm = TRUE),
                     min_wormlength_um = as.numeric(quantile(worm_length_um, na.rm = TRUE)[1]),
                     q10_wormlength_um = as.numeric(quantile(worm_length_um, probs = 0.1, na.rm = TRUE)[1]),
                     q25_wormlength_um = as.numeric(quantile(worm_length_um, probs = 0.25, na.rm = TRUE)[1]),
                     median_wormlength_um = median(worm_length_um, na.rm = T),
                     sd_wormlength_um = sd(worm_length_um, na.rm = T),
                     q75_wormlength_um = as.numeric(quantile(worm_length_um, probs = 0.75, na.rm = TRUE)[1]),
                     q90_wormlength_umF = as.numeric(quantile(worm_length_um, probs = 0.90, na.rm = TRUE)[1]),
                     max_wormlength_um = as.numeric(quantile(worm_length_um, na.rm = TRUE)[5]),
                     cv_wormlength_um = (sd_wormlength_um/mean_wormlength_um))
  
  sum_by <- raw_data %>% 
    dplyr::select(...)
  print(paste("SUMMARIZED BY", names(sum_by))) #altering you what is bring summarized

  return(list("raw_data" = raw_data, "processed_data" = processed_data, 
              "summarized_raw" = summarized_raw, "summarized_processed" = summarized_processed))
}

#run
final <- process(flag_labeled, Metadata_Plate, Metadata_Hour, Metadata_Well)

#pulling summarized_processed dataframe
processed_sum <- all$summarized_processed

#######################################
##### easyXpress wrapper function #####   # I though a wrapper function would be cool. Not sure how packages work but
#######################################   # this function would need to access all prior function in order to work properly 


### dont forget to load ALL prior functions from above... ####
#step 1. load and bind data
#step 2. model selection
#step 3. assign edge flags
#step 4. set all flags
#step 5. filter flags and process

#####
# this function will basically run each of the above functions in order in a single step
## just be sure to specify any additional arguments!!
Xpress <- function(filedir, ..., design = FALSE, model_num = 4, radius = 825, cluster_flag = TRUE, well_edge_flag = TRUE) {

  output <- read_data(filedir) %>%
    model_selection(raw_data = .) %>%
    edge_flag(raw_data = .) %>%
    set_flags(data = .) %>%
    process(flag_data = ., ...)
    
  return(output)
}


dirs <- "output_100w_model" #folder containing output data from cell profiler pipeline run.
total_run <- Xpress(dirs, Metadata_Plate, Metadata_Well)









