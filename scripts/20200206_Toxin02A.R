#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(png)
library(viridis)
library(pals)
library(rebus)

# remove outliers function (tukey's fences)
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#########################################
### 1: making get_data function       ###
#########################################

# load data: function might need to take a list of directories...maybe another function can make this list based on experiment name?
# load cell profiler data from full dose worm model (can load all .csv files with NonOverlappingWorms?)
worm_dat_full <- data.table::fread("data/20200206_toxin02A_Analysis01/20200206_toxin02_all_full.csv") %>%
  dplyr::mutate(model = "full_dose")

worm_dat_control <- data.table::fread("data/20200206_toxin02A_Analysis01/20200206_toxin02_all_control.csv") %>%
  dplyr::mutate(model = "control")

worm_dat_high_dose<- data.table::fread("data/20200206_toxin02A_Analysis01/20200206_toxin02_all_high.csv") %>%
  dplyr::mutate(model = "high_dose")

# load tratements
design <- data.table::fread("data/20200206_toxin02A_Analysis01/20200206_Toxin02A_design.csv", header = T)

# load image file names
file_df <- data.table::fread('data/20200206_toxin02A_Analysis01/20200206_toxin02A_metadata.csv') %>%
  dplyr::select(FileName_RawBF = Image_FileName_RawBF, Metadata_Plate, Metadata_Well) %>%
  tidyr::separate(Metadata_Well, into = c("Metadata_Row", "Metadata_Column"), sep = 1, remove = FALSE) # IF no row and column data are supplied

# join all
raw_dat <- full_join(worm_dat_full, worm_dat_control) %>%
  full_join(worm_dat_high_dose) %>%
  full_join(design) %>%
  full_join(file_df) %>%
  dplyr::group_by(Metadata_Well, Metadata_Plate) %>%
  dplyr::mutate(worm_length_um = 3.2937*Worm_Length) %>%
  dplyr::ungroup() %>%
  dplyr::select(model, names(design), FileName_RawBF, everything()) 

############################################
### 2: flag_dat: Flagging raw data       ###
############################################
# identify which model to choose for a particular parent worm object and create variables to tell if each model produces clusters.
# NOTE: "conflict file" could be used to specify each of the 27 permutations rather than this long list of shit below.
model_selection <- raw_dat %>%
  group_by(Metadata_Plate, Metadata_Well, Parent_WormObjects, model) %>%
  dplyr::mutate(num_worms = n()) %>%
  group_by(Metadata_Plate, Metadata_Well, Parent_WormObjects) %>%
  dplyr::distinct(model, .keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::select(Metadata_Plate, Metadata_Well, Metadata_Column, Metadata_Row, Parent_WormObjects, model, num_worms) %>%
  tidyr::spread(model, num_worms) %>%
  dplyr::select(-`<NA>`) %>%
  dplyr::mutate(model_code = case_when(
    control == 1 & full_dose == 1&   high_dose == 1 ~ "111",   control >= 2& full_dose == 1&   high_dose == 1 ~ "211",   is.na(control)& full_dose == 1& high_dose == 1 ~ "011",
    control == 1 & full_dose == 1&   high_dose >= 2 ~ "112",   control >= 2& full_dose == 1&   high_dose >= 2 ~ "212",   is.na(control)& full_dose == 1& high_dose >= 2 ~ "012",
    control == 1 & full_dose == 1&   is.na(high_dose) ~ "110", control >= 2& full_dose == 1&   is.na(high_dose) ~ "210", is.na(control)& full_dose == 1& is.na(high_dose) ~ "010",
    control == 1 & full_dose >= 2&   high_dose == 1 ~ "121",   control >= 2& full_dose >= 2&   high_dose == 1 ~ "221",   is.na(control)& full_dose >= 2& high_dose == 1 ~ "021",
    control == 1 & full_dose >= 2&   high_dose >= 2 ~ "122",   control >= 2& full_dose >= 2&   high_dose >= 2 ~ "222",   is.na(control)& full_dose >= 2& high_dose >= 2 ~ "022",
    control == 1 & full_dose >= 2&   is.na(high_dose) ~ "120", control >= 2& full_dose >= 2&   is.na(high_dose) ~ "220", is.na(control)& full_dose >= 2& is.na(high_dose) ~ "020",
    control == 1 & is.na(full_dose)& high_dose == 1 ~ "101",   control >= 2& is.na(full_dose)& high_dose == 1 ~ "201",   is.na(control)& is.na(full_dose)& high_dose == 1 ~ "001",
    control == 1 & is.na(full_dose)& high_dose >= 2 ~ "102",   control >= 2& is.na(full_dose)& high_dose >= 2 ~ "202",   is.na(control)& is.na(full_dose)& high_dose >= 2 ~ "002",
    control == 1 & is.na(full_dose)& is.na(high_dose) ~ "100", control >= 2& is.na(full_dose)& is.na(high_dose) ~ "200", is.na(control)& is.na(full_dose)& is.na(high_dose) ~ "000")) %>% #assign model code for all possible permutations of models
  dplyr::mutate(model_select = case_when(model_code == "001" ~ "high_dose",
                                         model_code == "002" ~ "high_dose",
                                         model_code == "010" ~ "full_dose",
                                         model_code == "011" ~ "full_dose",
                                         model_code == "012" ~ "full_dose",
                                         model_code == "020" ~ "full_dose",
                                         model_code == "021" ~ "full_dose",
                                         model_code == "022" ~ "full_dose",
                                         model_code == "100" ~ "control",
                                         model_code == "101" ~ "control",
                                         model_code == "102" ~ "control",
                                         model_code == "110" ~ "control",
                                         model_code == "111" ~ "control",
                                         model_code == "112" ~ "control",
                                         model_code == "120" ~ "control",
                                         model_code == "121" ~ "control",
                                         model_code == "122" ~ "control",
                                         model_code == "200" ~ "control",
                                         model_code == "201" ~ "control",
                                         model_code == "202" ~ "control",
                                         model_code == "210" ~ "control",
                                         model_code == "211" ~ "control",
                                         model_code == "212" ~ "control",
                                         model_code == "220" ~ "control",
                                         model_code == "221" ~ "control",
                                         model_code == "222" ~ "control")) %>%
  dplyr::mutate(model_flag = case_when(model_code == "002" ~ "high_dose_cluster",
                                       model_code == "012" ~ "high_dose_cluster",
                                       model_code == "020" ~ "full_dose_cluster",
                                       model_code == "021" ~ "full_dose_cluster",
                                       model_code == "022" ~ "full_dose_cluster",
                                       model_code == "102" ~ "high_dose_cluster",
                                       model_code == "112" ~ "high_dose_cluster",
                                       model_code == "120" ~ "full_dose_cluster",
                                       model_code == "121" ~ "full_dose_cluster",
                                       model_code == "122" ~ "full_dose_cluster",
                                       model_code == "200" ~ "control_cluster",
                                       model_code == "201" ~ "control_cluster",
                                       model_code == "202" ~ "control_cluster",
                                       model_code == "210" ~ "control_cluster",
                                       model_code == "211" ~ "control_cluster",
                                       model_code == "212" ~ "control_cluster",
                                       model_code == "220" ~ "control_cluster",
                                       model_code == "221" ~ "control_cluster",
                                       model_code == "222" ~ "control_cluster")) %>%
  dplyr::mutate(cluster_flag = case_when(model_code == "001" ~ F,
                                         model_code == "002" ~ T,
                                         model_code == "010" ~ F,
                                         model_code == "011" ~ F,
                                         model_code == "012" ~ F,
                                         model_code == "020" ~ T,
                                         model_code == "021" ~ T,
                                         model_code == "022" ~ T,
                                         model_code == "100" ~ F,
                                         model_code == "101" ~ F,
                                         model_code == "102" ~ F,
                                         model_code == "110" ~ F,
                                         model_code == "111" ~ F,
                                         model_code == "112" ~ F,
                                         model_code == "120" ~ F,
                                         model_code == "121" ~ F,
                                         model_code == "122" ~ F,
                                         model_code == "200" ~ T,
                                         model_code == "201" ~ T,
                                         model_code == "202" ~ T,
                                         model_code == "210" ~ T,
                                         model_code == "211" ~ T,
                                         model_code == "212" ~ T,
                                         model_code == "220" ~ T,
                                         model_code == "221" ~ T,
                                         model_code == "222" ~ T)) %>%
  dplyr::rename(num_worms_control = control, num_worms_full_dose = full_dose, num_worms_high_dose = high_dose)

# join model selction dataframe to full data frame and filter by model selection.
# set well edge flag parameters for ploting and flagging
edge_flag_center_x <- 1024
edge_flag_center_y <- 1024
edge_flag_radius <- 825

flag_dat <- full_join(raw_dat, model_selection) %>%
  dplyr::mutate(Metadata_Date = as.factor(Metadata_Date)) %>%
  dplyr::filter(model == model_select) %>%
  dplyr::mutate(well_edge_flag = ifelse(sqrt((AreaShape_Center_X - edge_flag_center_x)^2 + (AreaShape_Center_Y - edge_flag_center_y)^2) <= edge_flag_radius, FALSE, TRUE)) %>%
  dplyr::select(model_select, model_code, cluster_flag, model_flag, well_edge_flag, model:FileName_RawBF, ObjectNumber, Parent_WormObjects,	Metadata_Column,	Metadata_Date, Metadata_Plate,
                Metadata_Row,	Metadata_Well, AreaShape_Area,	AreaShape_Center_X,
                AreaShape_Center_Y, Worm_Length,  worm_length_um) %>%
  dplyr::group_by(Metadata_Date, Metadata_Plate, Metadata_Well) %>%
  dplyr::mutate(raw_well_mean_CP_length_um = mean(worm_length_um, na.rm = T),
                raw_well_sd_CP_length_um = sd(worm_length_um, na.rm = T),
                raw_well_cv_CP_length_um = (raw_well_sd_CP_length_um/raw_well_mean_CP_length_um)*100,
                raw_well_worm_num = n()) %>%
  dplyr::ungroup()

##############################################################################################################
### 3: flagout_dat: calculate well outliers after filtering flags from above and calculate new well stats  ###
##############################################################################################################
# this should be turned into easyxpress::flag_out() function. Default is to remove cluster_flag only.

flag_rm_dat <- flag_dat %>%
  dplyr::filter(cluster_flag != T, well_edge_flag !=T) %>% # set flags to remove here
  dplyr::mutate(flag_removed = TRUE) %>% # record the flags that were removed for calculating flag removed phenos
  dplyr::group_by(Metadata_Date, Metadata_Plate, Metadata_Well) %>%
  dplyr::mutate(flag_rm_well_mean_CP_length_um = mean(worm_length_um, na.rm = T),
                flag_rm_well_sd_CP_length_um = sd(worm_length_um, na.rm = T),
                flag_rm_well_cv_CP_length_um = (flag_rm_well_sd_CP_length_um/flag_rm_well_mean_CP_length_um)*100,
                flag_rm_well_worm_num = n()) %>%
  dplyr::ungroup() %>%
  dplyr::full_join(flag_dat) %>%
  dplyr::mutate(flags_removed= "rm_cluster_flag, rm_well_edge_flag") %>% # record the flags that were removed for calculating flag removed phenos
  dplyr::mutate(flag_rm = ifelse(is.na(flag_rm_well_mean_CP_length_um), T, F)) # logical for if worm object pheno was removed or not

##############################################################################################################
### 4: well_outlier_rm_dat:  Identify well outliers for raw and flag_rm phenotypes if present              ###
##############################################################################################################
# Calculate well outliers for each pheno
well_outlier_rm_dat <- flag_rm_dat %>%
  dplyr::filter(flag_rm == F) %>%
  dplyr::group_by(Metadata_Date, Metadata_Plate, Metadata_Well) %>%
  dplyr::mutate(flag_rm_well_outlier_flag = (remove_outliers(worm_length_um)),
                flag_rm_well_outlier_flag = ifelse(is.na(flag_rm_well_outlier_flag), TRUE, FALSE)) %>%
  dplyr::filter(!flag_rm_well_outlier_flag == TRUE) %>% # Joy edit to handle outliers properly
  dplyr::mutate(well_outlier_flag_rm_well_mean_CP_length_um = mean(worm_length_um, na.rm = T),
                well_outlier_flag_rm_well_sd_CP_length_um = sd(worm_length_um, na.rm = T),
                well_outlier_flag_rm_well_cv_CP_length_um = (well_outlier_flag_rm_well_sd_CP_length_um/well_outlier_flag_rm_well_mean_CP_length_um)*100,
                well_outlier_flag_rm_well_worm_num = n()) %>%
  dplyr::ungroup() %>%
  dplyr::full_join(flag_rm_dat) %>%
  dplyr::group_by(Metadata_Date, Metadata_Plate, Metadata_Well) %>%
  dplyr::mutate(raw_well_outlier_flag = (remove_outliers(worm_length_um)),
                raw_well_outlier_flag = ifelse(is.na(raw_well_outlier_flag), TRUE, FALSE)) %>%
  dplyr::mutate(well_outlier_rm_raw_well_mean_CP_length_um = mean(worm_length_um, na.rm = T),
                well_outlier_rm_raw_well_sd_CP_length_um = sd(worm_length_um, na.rm = T),
                well_outlier_rm_raw_well_cv_CP_length_um = (well_outlier_rm_raw_well_sd_CP_length_um/well_outlier_rm_raw_well_mean_CP_length_um)*100,
                well_outlier_rm_raw_well_worm_num = n()) %>%
  dplyr::ungroup() %>%  #naming objects for plotting. Note, the flag_rm_well_outlier_flag is only an outlier for the well given the other flags are removed
  dplyr::mutate(object_type = case_when(
    cluster_flag == T & well_edge_flag == T & flag_rm_well_outlier_flag == T ~ "all_flags",
    cluster_flag == T & well_edge_flag == T & flag_rm_well_outlier_flag == F ~ "cluster/edge_flag",
    cluster_flag == T & well_edge_flag == F & flag_rm_well_outlier_flag == T ~ "cluster/outlier_flag",
    cluster_flag == T & well_edge_flag == F & flag_rm_well_outlier_flag == F ~ "cluster_flag",
    cluster_flag == F & well_edge_flag == T & flag_rm_well_outlier_flag == T ~ "edge/outlier_flag",
    cluster_flag == F & well_edge_flag == T & flag_rm_well_outlier_flag == F ~ "edge_flag",
    cluster_flag == F & well_edge_flag == F & flag_rm_well_outlier_flag == T ~ "outlier_flag",
    cluster_flag == T & well_edge_flag == T & is.na(flag_rm_well_outlier_flag) ~ "cluster/edge_flag",
    cluster_flag == T & well_edge_flag == F & is.na(flag_rm_well_outlier_flag) ~ "cluster_flag",
    cluster_flag == F & well_edge_flag == T & is.na(flag_rm_well_outlier_flag) ~ "edge_flag",
    cluster_flag == F & well_edge_flag == F & is.na(flag_rm_well_outlier_flag) & model_select == "control" ~ "control_model",
    cluster_flag == F & well_edge_flag == F & is.na(flag_rm_well_outlier_flag) & model_select == "full_dose" ~ "full_model",
    cluster_flag == F & well_edge_flag == F & is.na(flag_rm_well_outlier_flag) & model_select == "high_dose" ~ "high_model",
    cluster_flag == F & well_edge_flag == F & flag_rm_well_outlier_flag == F & model_select == "control" ~ "control_model",
    cluster_flag == F & well_edge_flag == F & flag_rm_well_outlier_flag == F & model_select == "full_dose" ~ "full_model",
    cluster_flag == F & well_edge_flag == F & flag_rm_well_outlier_flag == F & model_select == "high_dose" ~ "high_model"))

##############################################################################################################
### 4: Plotting dose responses                                                                             ###
##############################################################################################################

# set color palette
strain_colors <- c(PD1074 = "orange", CB4856 = "blue")

# plot dose responses raw data
raw_dose <- ggplot(well_outlier_rm_dat %>% dplyr::distinct(Metadata_Date, Metadata_Plate, Metadata_Well, .keep_all =T)) +
  aes(x = factor(concentration_um), y = raw_well_mean_CP_length_um, fill = strain, alpha = 0.85) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = strain_colors) +
  guides(alpha = FALSE) +
  ylim(0,850) +
  geom_point(aes(fill = strain), position = position_jitterdodge(jitter.width = 0.25), shape = 21) +
  facet_wrap(~drug, scales = "free") +
  theme_bw() +
  labs(x = "drug dose (uM)", title = "20200206_raw_dose")
raw_dose
ggsave('plots/20200206_raw_dose_reponse_trimodel.png', width = 10, height = 5)

well_outlier_removed_dose <- ggplot(well_outlier_rm_dat %>% dplyr::distinct(Metadata_Date, Metadata_Plate, Metadata_Well, .keep_all =T)) +
  aes(x = factor(concentration_um), y = well_outlier_flag_rm_well_mean_CP_length_um, fill = strain, alpha = 0.85) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = strain_colors) +
  guides(alpha = FALSE) +
  ylim(0,850) +
  geom_point(aes(fill = strain), position = position_jitterdodge(jitter.width = 0.25), shape = 21) +
  facet_wrap(~drug, scales = "free") +
  theme_bw() +
  labs(x = "drug dose (uM)", title = "20200131_well_outlier_removed_dose")
well_outlier_removed_dose
ggsave('plots/20200206_well_outlier_removed_dose_reponse_trimodel.png', width = 10, height = 5)

##############################################################################################################
### 5: view_well function                                                                                  ###
##############################################################################################################
# set image directory
proc_img_dir <- "data/20200206_toxin02A_Analysis01/ProcessedImages"

#define view_well function
view_well <- function(df, plate, well, proc_img_dir) {
  # use stringr to find middle part of directory that contains the date and experiment name.
  date <- DGT %R% DGT %R% DGT %R% DGT %R% DGT %R% DGT %R% DGT%R% DGT
  separator <- char_class("/_-")
  run <- "RUN" %R% DGT%R% DGT
  name1 <- separator %R% date %R% separator %R% one_or_more(WRD) %R% "_"
  name2 <- date %R% separator %R% one_or_more(WRD) %R% "_"
  dir_date_name <- str_extract(proc_img_dir, pattern = name1)
  date_name <- str_extract(proc_img_dir, pattern = name2)
  
  # define object type palette
  object_type_palette <- c("control_model" = "#2b66f8", #blue
                           "full_model" = "#BE0032", #red
                           "high_model" = "#36802d", #green
                           "all_flags" = "#C2B280", #tan
                           "outlier/edge_flag" = "#848482", #grey
                           "cluster/edge_flag" = "#F38400", #orange
                           "cluster/outlier_flag" = "#875692", #purple
                           "edge_flag" = "#F3C300", #yellow
                           "outlier_flag" = "#222222", #dary_grey
                           "cluster_flag" = "#F2F3F4") #white
  # subset dataframe
  plot_dat <- df %>%
    dplyr::mutate(Metadata_Plate = case_when(
      Metadata_Plate == 1 ~ "01",
      Metadata_Plate == 2 ~ "02",
      Metadata_Plate == 3 ~ "03",
      Metadata_Plate == 4 ~ "04",
      Metadata_Plate == 5 ~ "05",
      Metadata_Plate == 6 ~ "06",
      Metadata_Plate == 7 ~ "07",
      Metadata_Plate == 8 ~ "08",
      Metadata_Plate == 9 ~ "09")) %>%
    dplyr::filter(Metadata_Plate == glue::glue({plate}), Metadata_Well == glue::glue({well}))
  
  img <- readPNG(glue::glue("{proc_img_dir}{dir_date_name}{plate}_2x_{well}_overlay.png"))
  
  h<-dim(img)[1] # image height
  w<-dim(img)[2] # image width
  
  well_img <- ggplot(plot_dat) +
    aes(x = AreaShape_Center_X, y = AreaShape_Center_Y, fill = object_type) +
    scale_fill_manual(values = object_type_palette) +
    annotation_custom(grid::rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc")), 0, w, 0, -h) + # The minus is needed to get the y scale reversed. these are pixel dim of image
    scale_x_continuous(expand=c(0,0),limits=c(0,w)) +
    scale_y_reverse(expand=c(0,0),limits=c(h,0)) +  # The y scale is reversed because in image the vertical positive direction is typically downward.Also note the limits where h>0 is the first parameter.
    labs(x = "", y = "", title = glue::glue("{plot_dat %>% pull(Metadata_Date)}_Plate{plate}_{well}_{plot_dat %>% pull(drug)}_{plot_dat %>% pull(concentration_um)}uM_{plot_dat %>% pull(strain)}"), fill = "Object Class") +
    geom_point(shape = 21, alpha = 0.75) +
    coord_equal() +
    theme(legend.position = "none") +
    annotate("path",
             x=edge_flag_center_x+edge_flag_radius*cos(seq(0,2*pi,length.out=100)),
             y=edge_flag_center_y+edge_flag_radius*sin(seq(0,2*pi,length.out=100)), color = "red", alpha = 0.25)
  
  well_img_box_plot <- ggplot(plot_dat) +
    geom_boxplot(aes(x = Metadata_Well, y = worm_length_um), outlier.shape = NA) +
    geom_jitter(shape = 21, width = 0.25, size = 3, aes(x = Metadata_Well, y = worm_length_um, fill = object_type)) +
    scale_fill_manual(values = object_type_palette) +
    labs(x ="") +
    ylim(0,1200) +
    theme_bw() +
    theme(legend.position = "right")
  
  
  full_diagnostic <- cowplot::plot_grid(well_img, well_img_box_plot, nrow = 1, rel_widths = c(1,0.25), align = "hv", axis = "tb")
  full_diagnostic
  
  ggsave(glue::glue('plots/{date_name}_Plate{plate}_{well}_{plot_dat %>% distinct(drug) %>% pull(drug)}_{plot_dat %>% distinct(concentration_um) %>% pull(concentration_um)}uM_{plot_dat %>% distinct(strain) %>% pull(strain)}_diagnostic.png'), width = 16, height = 10)
  full_diagnostic
}

# testing view_well function
view_well(well_outlier_rm_dat, "02", "E02", "data/20200206_toxin02A_Analysis01/ProcessedImages")

##############################################################################################################
### 5: view_dose function   NOT WORKING YET                                                                ###
##############################################################################################################

view_dose <- function(df, drug, proc_img_dir) {
  # use stringr to find middle part of directory that contains the date and experiment name.
  date <- DGT %R% DGT %R% DGT %R% DGT %R% DGT %R% DGT %R% DGT%R% DGT
  separator <- char_class("/_-")
  run <- "RUN" %R% DGT%R% DGT
  name1 <- separator %R% date %R% separator %R% one_or_more(WRD) %R% "_"
  name2 <- date %R% separator %R% one_or_more(WRD) %R% "_"
  dir_date_name <- str_extract(proc_img_dir, pattern = name1)
  date_name <- str_extract(proc_img_dir, pattern = name2)
  
  # define object type palette
  object_type_palette <- c("control_model" = "#2b66f8", #blue
                           "full_model" = "#BE0032", #red
                           "high_model" = "#36802d", #green
                           "all_flags" = "#C2B280", #tan
                           "outlier/edge_flag" = "#848482", #grey
                           "cluster/edge_flag" = "#F38400", #orange
                           "cluster/outlier_flag" = "#875692", #purple
                           "edge_flag" = "#F3C300", #yellow
                           "outlier_flag" = "#222222", #dary_grey
                           "cluster_flag" = "#F2F3F4") #white
  # subset dataframe
  plot_dat <- df %>%
    dplyr::filter(drug == glue::glue({drug}))
  dplyr::group_by(concentration_um) %>%
    dplyr::distinct(concentration_um, .keep_all = T) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Metadata_Plate = case_when(
      Metadata_Plate == 1 ~ "01",
      Metadata_Plate == 2 ~ "02",
      Metadata_Plate == 3 ~ "03",
      Metadata_Plate == 4 ~ "04",
      Metadata_Plate == 5 ~ "05",
      Metadata_Plate == 6 ~ "06",
      Metadata_Plate == 7 ~ "07",
      Metadata_Plate == 8 ~ "08",
      Metadata_Plate == 9 ~ "09",
      Metadata_Plate >= 10 ~ as.character(Metadata_Plate))) %>%
    dplyr::arrange(concentration_um)
  
  # Pull image names
  img_names <- plot_dat %>%
    lpull(FileName_RawBF) %>%
    stringr::str_replace_all(., pattern = "-", replacement = "_") %>%
    stringr::str_replace(., pattern = ".TIF", replacement = "_overlay.png")
  
  # make list of images for dose response 
  dose_img_list <- list()
  for(i in 1:length(unique(img_names))){
    img <- readPNG(glue::glue("{proc_img_dir}/{img_names[i]}"))
    dose_img_list[[i]] <- img
  }
  
  # find dimensions of images
  h<-dim(dose_img_list[[1]])[1] # image height
  w<-dim(dose_img_list[[1]])[2] # image width
  
  
  testdf <- well_outlier_rm_dat %>% dplyr::filter(Metadata_Plate == "7", Metadata_Well=="A01")
  #  plot well images for doses
  test_0 <- ggplot(well_outlier_rm_dat %>% dplyr::filter(Metadata_Well ==  test$Metadata_Well[1] & Metadata_Plate %in% test$Metadata_Plate[1])) +
    aes(x = AreaShape_Center_X, y = AreaShape_Center_Y, fill = object_type) +
    scale_fill_manual(values = object_type_palette) +
    annotation_custom(grid::rasterGrob(dose_img_list[[1]], width=unit(1,"npc"), height=unit(1,"npc")), 0, w, 0, -h) + # The minus is needed to get the y scale reversed. these are pixel dim of image
    scale_x_continuous(expand=c(0,0),limits=c(0,w)) +
    scale_y_reverse(expand=c(0,0),limits=c(h,0)) +  # The y scale is reversed because in image the vertical positive direction is typically downward.Also note the limits where h>0 is the first parameter.
    labs(x = "", y = "", title = glue::glue("{test %>% pull(Metadata_Date)}_Plate{test %>% pull(Metadata_Plate)}_{test %>% pull(Metadata_Well)}_{test %>% pull(drug)}_{test %>% pull(concentration_um)}uM_{test %>% pull(strain)}"), fill = "Object Class") +
    geom_point(shape = 21, alpha = 0.75) +
    coord_equal() +
    theme(legend.position = "none") +
    annotate("path",
             x=edge_flag_center_x+edge_flag_radius*cos(seq(0,2*pi,length.out=100)),
             y=edge_flag_center_y+edge_flag_radius*sin(seq(0,2*pi,length.out=100)), color = "red", alpha = 0.25)
  test_0

##############################################################################################################
### 6: view_plate function                                                                                 ###
##############################################################################################################

# plot the plate
plate = "02"
plot_dat <- well_outlier_rm_dat %>%
  dplyr::mutate(Metadata_Plate = case_when(
    Metadata_Plate == 1 ~ "01",
    Metadata_Plate == 2 ~ "02",
    Metadata_Plate == 3 ~ "03",
    Metadata_Plate == 4 ~ "04",
    Metadata_Plate == 5 ~ "05",
    Metadata_Plate == 6 ~ "06",
    Metadata_Plate == 7 ~ "07",
    Metadata_Plate == 8 ~ "08",
    Metadata_Plate == 9 ~ "09",
    Metadata_Plate >= 10 ~ as.character(Metadata_Plate))) %>%
  dplyr::filter(Metadata_Plate == glue::glue({plate})) %>%
  dplyr::distinct(Metadata_Well, .keep_all = T)

require(plotly)
ggplotly(ggplot(plot_dat) +
           aes(x = factor(Metadata_Column, levels = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")), y = factor(Metadata_Row, levels = c("H", "G", "F", "E", "D", "C", "B", "A")),
               color = well_outlier_flag_rm_well_mean_CP_length_um, text = Metadata_Well, key = well_outlier_flag_rm_well_mean_CP_length_um) +
           geom_point( size = 2) +
           theme_bw() +
           theme(legend.position = "right") + 
           labs(x = "Column", y = "Row", title = glue::glue("Plate{plate}")) +
           scale_color_viridis(),
         tooltip = "key")
