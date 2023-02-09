### About -----------------------------------------------------------------------------------'

# This script attaches metadata to the data produced by the Norra Midsjöbank Project 
# It is modified to run on the data produced from the Norra Midsjöbank Project
# The accompanying script -- 5.0_metadata.R -- is needed to run this script

#by: Oscar Törnqvist
#modified by: Francis Freire

####

###########################################################################################
rm(list = ls())

library(raster)
library(tidyverse)
# Include metadata script
source("M:/marin/nkp16/nmi18/arbetskataloger/habmod/r-project/R/5.0_metadata.R")

use_packages<-function(package_list) {
  for (p in package_list){
    if(!suppressMessages(require(p, character.only = TRUE))){
      install.packages(p)
      suppressMessages(library(p, character.only = TRUE))
    }
  }
} 

#All tags... execept those generated from the raster, i e. extent and cell size

# These five  parameters are originally set by Melody:
ref_syst <- curr_ref_syst <- "+init=epsg:32633"
#product_title <- "Geological/Substrate maps of Norra Midsjöbank"
friendly_title <- "Long description of product"
model_area <- "Baltic sea"
product_version <- "1"

product_language_code <- "eng"
product_character_set <- "utf8"
product_metadata_standard_name <- "SS-EN ISO 19115:2005-geodata.se version 4.0"
product_metadata_standard_version <- "4.0"

#product_abstract <- "Gridded geological substrate map produced for Norra Midsjöbank project. The substrate map was created using data derived from bio-geophysical data processed using machine learning algorithms (i.e., Boosted regression tree) to produce the raster substrate maps. The surface is gridded at 5, 10, 25, 50, 100, 250 m resolution. Data collection and processing details are found at the project."
product_purpose <- "Maps for conservation and management"

product_base_url <- "http://www.sgu.se"
product_sub_url <- "/produkter/"
product_org <- "Geological Survey of Sweden (SGU)"
product_postal1 <- "Box 670"
product_postal2 <- "751 28"
product_city <- "Uppsala"
product_country <- "Sweden"

product_phone <- "+4618 17 90 00"
product_contact <- "Customer support"
product_email <- "kundservice@sgu.se"
product_access_constraint = "Open data"
product_credits = c("SGU")

product_distribution_link <- "http://somelink"
#product_distribution_name <- "Substrate Maps for Norra Midsjöbank"
product_distribution_description <- "description"
product_distribution_protocol <- "WWW:LINK-1.0-http--link"

product_topic_category_codes<-c("environment","geoscientificInformation","oceans")
product_update_frequency <- "Final Product"
product_use_limitation  <- c("not for use as navigation")
product_access_constraints <-"Open access data"
product_secrecy = "Open acces"
product_secrecy_notes = 
product_secrecy_classification = "Förordning (2016:320) om skydd för geografisk information"
product_secrecy_handling = "Använd endast med säkerhetsklassad IT-utrustning"
#product_theme_keywords = c("Oceanographic geographical features","Geology","Habitat mapping")
product_place_keywords = "Baltic"

product_data_type = "grid"

product_supplemental_info <- "Habitat maps developed for Swedish Agency for Marine and Water Management (HAVS- OCH VATTEN MYNDIGHETEN) as part of the National Marine Mapping Project (NMK) objectives"

product_data_quality_check <- "Data Quality check" #Method, protocol, standard
product_data_quality_check_report <- "QC passed" #Result
product_quality_conformance <- "Model quality evaluation supplied in project report"
#product_lineage_statement <- "The data was created by using predictor variables derived from various geo-physical methods and analyzed using machine learning algorithms to predict and create the raster surfaces."

product_legal_constraints_url <-  "http://inspire.ec.europa.eu/metadata-codelist/ConditionsApplyingToAccessAndUse/noConditionsApply"
product_legal_constraints_desc <- "No conditions apply to access and use."
inspire_spec1_title <- "Commission Regulation (EU) No 1089/2010 of 23 November 2010 implementing Directive 2007/2/EC of the European Parliament and of the Council as regards interoperability of spatial data sets and services"
inspire_spec2_title <- "COMMISSION REGULATION (EC) No 1205/2008 of 3 December 2008 implementing Directive 2007/2/EC of the European Parliament and of the Council as regards metadata"
inspire_explanation <- "See the referenced specification"

# "Product name" and final raster to ge given XML metadata

#pred_name <-"Hoburg HUB level X or the like" #predicted parameter, e.g. "soft bottom" of "Substrate" product

# set paths to directory ------------------------------------------------------------------------------------------------------------

# set directory for each model theme and version
version <- "v11" #v00, v01 ...
prj_names <- c("sea18", "nmi18", "hob16") # list projects to combine
i=2 #nmi18
prj <- prj_names[i]

#### add metadata to the group outputs directory ------------------------------------------------------------------------------------------------------------

# A.) To add metadata to all the geology/substrate final outputs

input_dir_sub <- "/data/outputs/predictions/geology/substrate/"
input_dir_sub <- paste0(getwd(), input_dir_sub, prj, "/", version, "/", "delivery", sep="")

all_files_sub <- input_dir_sub %>%
  list.files(pattern = ".tif$", full.names = TRUE, recursive = TRUE)

for (each_file in all_files_sub) {
  #these metadata values were specific to the type of map, (i.e. geology maps)
  product_title <- "Geological/Substrate maps of Norra Midsjöbank"
  product_abstract <- "Gridded geological substrate map produced for Norra Midsjöbank project. The substrate map was created using data derived from bio-geophysical data and then processed using machine learning algorithms (i.e., Boosted regression tree) to produce the raster the substrate maps. The surface is gridded at 5, 10, 25, 50, 100, 250 m resolution. Data collection and processing details are found in the project report."
  product_distribution_name <- "Substrate Maps for Norra Midsjöbank"
  product_theme_keywords = c("Oceanographic geographical features","Geology","Habitat mapping")
  product_lineage_statement <- "The data was created by using predictor variables derived from various geo-physical methods and analyzed using machine learning algorithms to predict and create the raster surfaces."
  
  pred_name <- basename(each_file)
  pred_name1 <- paste("Predicted ", str_to_title(gsub("\\d.*|sgu.*", "", (gsub("_", " ", pred_name)))), "areas (@ ", gsub(".tif","", sub(".*_", "", pred_name)), " resolution)", " in Norra Midsjöbank", sep = "")
  xml_file_path <<- paste0(each_file, ".xml")
  generate_metadata(each_file, pred_name)
}


# B.) To add metadata to all the biology final outputs

input_dir_bio <- "/data/outputs/predictions/biology/"
input_dir_bio <- paste0(getwd(), input_dir_bio, prj, "/", version, "/", "delivery", sep="")

all_files_bio <- input_dir_bio %>%
  list.files(pattern = ".tif$", full.names = TRUE, recursive = TRUE)

for (each_file in all_files_bio) {
  #these metadata values were specific to the type of map, (i.e. biology maps)
  product_title <- "Biological maps of Norra Midsjöbank"
  product_abstract <- "Gridded biological distribution map of habitat classes for Norra Midsjöbank project. The map was created using data derived from bio-geophysical data and then processed using machine learning algorithms (i.e., Boosted regression tree) to produce the raster the substrate maps. The surface is gridded at 5, 10, 25, 50, 100, 250 m resolution. Data collection and processing details are found in the project report."
  product_distribution_name <- "Biological Maps for Norra Midsjöbank"
  product_theme_keywords = c("Oceanographic geographical features","Biological distribution","Habitat mapping")
  product_lineage_statement <- "The data was created by using predictor variables derived from various geo-physical methods and analyzed using machine learning algorithms to predict and create the raster surfaces."
  
  pred_name <- basename(each_file)
  pred_name1 <- paste("Predicted ", str_to_title(gsub("\\d.*|sgu.*", "", (gsub("_", " ", pred_name)))), "areas (@ ", gsub(".tif","", sub(".*_", "", pred_name)), " resolution)", " in Norra Midsjöbank", sep = "")
  xml_file_path <<- paste0(each_file, ".xml")
  generate_metadata(each_file, pred_name)
}


# C.) To add metadata to all the HELCOM HUB final outputs

input_dir_classifications_hub <- "/data/outputs/predictions/classifications/hub/"
input_dir_classifications_hub <- paste0(getwd(), input_dir_classifications_hub, prj, "/", version, sep="")

all_files_classifications_hub <- input_dir_classifications_hub %>%
  list.files(pattern = ".tif$", full.names = TRUE, recursive = TRUE)

for (each_file in all_files_classifications_hub) {
  #these metadata values were specific to the type of map, (i.e. HELCOM HUB maps)
  product_title <- "HELCOM HUB maps of Norra Midsjöbank"
  product_abstract <- "Gridded HELCOM HUB map of Norra Midsjöbank. The map was created using data derived from bio-geophysical data and then processed using machine learning algorithms (i.e., Boosted regression tree) to produce the raster the substrate maps. The surface map represents various levels of the HELCOM HUB classification and is gridded at 5, 10, 25, 50, 100, 250 m resolution. Data collection and processing details are found in the project report."
  product_distribution_name <- "HELCOM HUB Maps of Norra Midsjöbank"
  product_theme_keywords = c("Oceanographic geographical features","HELCOM classification","Habitat mapping")
  product_lineage_statement <- "The data was created by using predictor variables derived from various geo-physical methods and analyzed using machine learning algorithms to predict and create the raster surfaces."
  
  pred_name <- basename(each_file)
  pred_name1 <- paste("Predicted HELCOM ", str_to_title(gsub("sgu.*", "", (gsub("_", " ", pred_name)))), "(@ ", gsub(".tif","", sub(".*_", "", pred_name)), " resolution)", " in Norra Midsjöbank", sep = "")
  xml_file_path <<- paste0(each_file, ".xml")
  generate_metadata(each_file, pred_name)
}


# D.) To add metadata to all the Natura 2000 final outputs

input_dir_classifications_natura2000 <- "/data/outputs/predictions/classifications/natura2000/"
input_dir_classifications_natura2000 <- paste0(getwd(), input_dir_classifications_natura2000, prj, "/", version, sep="")

all_files_classifications_natura2000 <- input_dir_classifications_natura2000 %>%
  list.files(pattern = ".tif$", full.names = TRUE, recursive = TRUE)

for (each_file in all_files_classifications_natura2000) {
  #these metadata values were specific to the type of map, (i.e. natura200) 
  product_title <- "Natura 2000 maps of Norra Midsjöbank"
  product_abstract <- "Gridded Natura 2000 map of Norra Midsjöbank. The map was created using data derived from bio-geophysical data and then processed using machine learning algorithms (i.e., Boosted regression tree) to produce the raster the substrate maps. The surface is gridded at 5, 10, 25, 50, 100, 250 m resolution. Data collection and processing details are found in the project report."
  product_distribution_name <- "Natura 2000 Maps of Norra Midsjöbank"
  product_theme_keywords = c("Oceanographic geographical features","Natura 2000 classification","Habitat mapping")
  product_lineage_statement <- "The data was created by using predictor variables derived from various geo-physical methods and analyzed using machine learning algorithms to predict and create the raster surfaces."
  
  pred_name <- basename(each_file)
  pred_name1 <- paste("Predicted ", str_to_title(gsub("sgu.*", "", (gsub("_", " ", pred_name)))), "map (@ ", gsub(".tif","", sub(".*_", "", pred_name)), " resolution)", " in Norra Midsjöbank", sep = "")
  xml_file_path <<- paste0(each_file, ".xml")
  generate_metadata(each_file, pred_name)
  
}


# E.) To add metadata to Backscatter

input_dir_bs <- paste0(getwd(), "/data/inputs/predictors/nmi18/bs/", sep="")

all_files_bs <- input_dir_bs %>%
  list.files(pattern = ".tif$", full.names = TRUE, recursive = TRUE)

for (each_file in all_files_bs) {
  product_title <- "Backscatter maps of Norra Midsjöbank"
  product_abstract <- "Gridded Backscatter map of Norra Midsjöbank. The map was created from raw Multibeam data (EM2040 sytem) that was processed using QPS-Fledermaus software to produce Backscatter products and Angular range (ARA) derivatives. Data processing details are found in the project report."
  product_distribution_name <- "Backscatter Maps of Norra Midsjöbank"
  product_theme_keywords = c("Oceanographic geographical features","Multibeam Backscatter","Habitat mapping")
  pred_name <- basename(each_file)
  pred_name1 <- paste("Predicted ", str_to_title(sub("bs1m |bs5m ", "Backscatter product ",(gsub(".tif", "", (gsub("_", " ", pred_name)))))), " resolution)", " in Norra Midsjöbank", sep = "")
  xml_file_path <<- paste0(each_file, ".xml")
  generate_metadata(each_file, pred_name)
  
}

# F.) To add metadata to Depth

input_dir_bs <- paste0(getwd(), "/data/inputs/predictors/nmi18/depth/", sep="")

all_files_bs <- input_dir_bs %>%
  list.files(pattern = ".tif$", full.names = TRUE, recursive = TRUE)

for (each_file in all_files_bs) {
  product_title <- "Bathymetry map of Norra Midsjöbank"
  product_abstract <- "Gridded Bathymetry map of Norra Midsjöbank. The map was created from raw Multibeam data (EM2040 sytem) that was processed using Caris software. Data collection and processing details are found in the project report."
  product_distribution_name <- "Backscatter Maps of Norra Midsjöbank"
  product_theme_keywords = c("Oceanographic geographical features","Bathymetry","Habitat mapping")
  pred_name <- basename(each_file)
  pred_name1 <- paste("Bathymetry Map of Norra Midsjöbank")
  xml_file_path <<- paste0(each_file, ".xml")
  generate_metadata(each_file, pred_name)
  
}

# F.) To add metadata to Points (.csv file)

input_dir_points <- paste0(getwd(), "/data/inputs/points/", sep="")

all_files_points <- input_dir_points %>%
  list.files(pattern = ".csv$", full.names = TRUE, recursive = TRUE)




for (each_file in all_files_bs) {
  product_title <- "Ground Validation points for the Norra Midsjöbank Project"
  product_abstract <- "The file contain the ground validation (gv) points used in the Norra Midsjöbank Project. The data was collected during the nmi18 Fieldwork using a modified Point intercept transect method. Analysis was done in the office. For details, please refer to the project report."
  product_distribution_name <- "Ground validation (GV) point file for Norra Midsjöbank"
  product_theme_keywords = c("Oceanographic geographical features","Point Intercept","Ground Validation", "Habitat mapping")
  pred_name <- basename(each_file)
  pred_name1 <- paste("Ground Validation point file for Norra Midsjöbank")
  xml_file_path <<- paste0(each_file, ".xml")
  generate_metadata(each_file, pred_name)
  
}

####add shape, creative common license cc0