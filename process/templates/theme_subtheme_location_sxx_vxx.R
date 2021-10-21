# ========================================================================================================= #
#### ABOUT ####
# ========================================================================================================= #
# Brief description: 
# Script by: Initial, org
# Updated: Initial, org, comments
# Developer check: Initial, org, yymmdd
# External check: Initial, org, yymmdd, comments
# ========================================================================================================= #
#### INSTRUCTIONS ####
# ========================================================================================================= #    
# Style
  # Write code that is easy to follow and repeat.
  # Suggested R style is tidyverse. Style guide here https://style.tidyverse.org/index.html

# Metadata
  # All indata from "data_raw" must have a unique ID with *_metasym.txt file to go along, check before you add data.
  # All indata from and outputs to "data" must have a *_sourcesym.txt file with all (accumulated) sourceIDs used to create the layer.
  # The metasym files and tracking of those in the sourcesym files ensures correct sources for products and ISO metadata, this template helps you get it right..
  
# Naming and folders
  # Predefined standard names for main/sub folders and final products in two text files that can be added to in need.
  # ..\shiny_data_upload\modify_txt_files_v01.2.xlsx and ..\process\templates\modify_txt_files_data_v01.0.xlsx
  # Main version (e.g. v01 for 2021) follows dev cycles for all WIOSYm. For version control within v01 use V01.1, V01.2 etc..
  # More code/instructions to standardize folders and names under "set destinations" section
  # Each product dir have a "proc" folder. It's the temp space where all workfiles are stored. No final products here.  
  # All final products are stored under the root folder (e.g. ..\data\reg\eco\bh_coral\v01) together with 
  # a nameofmyfile_sourcesym.txt file with all source IDs used. No other files here. move products to _archive\ dir if defunct 

# ========================================================================================================= # 
#### PREPARATIONS ####
# ========================================================================================================= #

# Read packages ---------------------------------------------------------------------------------------------
# list all packages that this script is dependent on, remove any redundant ones
x <- c("sf", "raster", "rgdal", "tidyverse")

#install.packages(x, dependencies=TRUE)
                                            
library(sf)
library(raster)
library(rgdal)
library(tidyverse)

# Set main destinations ----------------------------------------------------------------------------------------

# set version
# INPUT
version = "v01" # relates to major iteration of WIOSym (v01 for grid 2021, v02 for grid 2022)
d_version = ".0"  # detailed version (you can add another level e.g. 0.0  -> v01.0.0)

# set geographic scope (regional for whole WIO area)
read_tsv("./shiny_data_upload/locations.txt")  # prints available locations to choose from and their abbreviations
# INPUT
location = "reg" # choose applicable 3 letter abbreviation from "location_val" (regional wiosym = "reg")

# Set names and working folder using standard themes, if available
folders_raw <- read_tsv("./shiny_data_upload/folders.txt")
folders_data <- read_tsv("./shiny_data_upload/folders_data.txt")
folders <- folders_raw %>% 
  add_row(folders_data)

# select main theme
unique(folders["theme_folder"])
# INPUT
theme <- "theme" # input your main theme e.g. "eco" (choose from list in "theme_folder")

# select sub theme / component
folders %>%
  filter(theme_folder == theme) %>% 
  dplyr::select(subtheme_folder, subtheme_name) %>% 
  print(n = Inf)
# INPUT
subtheme <- "subtheme_folder" # input your subtheme based on the most appropriate option from the list. If your sub theme is missing, please add new sub theme to folders.txt first (and make sure to sync with shared location for WIOSym if applicable)


# Create destination path and folder ------------------------------------------------------------------------------------------------
dest_path <- paste("./data/", location, "/", theme, "/", subtheme, "/", version, "/", sep="")  # location of final products
dest_path

if (!dir.exists(dest_path)){
  dir.create(dest_path, recursive = TRUE)
}


# create work directory for any temporary files
work_path <- paste(dest_path, "proc/", sep="") # in this folder you can organise any files you use in the process using your preferred structure. No "final products" here, they go under root (path_dest1)
work_path

if (!dir.exists(work_path)){
  dir.create(work_path, recursive = TRUE)
}

# process log directory, for saving process logs and the exact script (or mxd etc) used to create the output files
proc_path <- paste(dest_path, "proc_log/", sep="") # in this folder you can organize any files you use in the process using your prefferred structure. No "final products" here, they go under root (path_dest1)
proc_path

if (!dir.exists(proc_path)){
  dir.create(proc_path, recursive = TRUE)
}

# archive directory for any redundant products when this script is updated...
archive_path <- paste(dest_path, "_archive/", sep="") # in this folder you can organize any files you use in the process using your prefferred structure. No "final products" here, they go under root (path_dest1)
archive_path
if (!dir.exists(archive_path)){
  dir.create(archive_path, recursive = TRUE)
}


# Save your R script ---------------------------------------------------------------------------------------------------------------------------------------
#location for all active R scripts (latest version) according to theme, subtheme
script_path <- "./process/r"  # main path
script_path
if (!dir.exists(script_path)){
  dir.create(script_path, recursive = TRUE)
}

# run lines below, copy path and name and use File/Save As in R studio to save your .R workfile in the correct workspace - 
# Note: in the end of the script (e.g. "write to file section") you are also asked to save a copy of your final code in the destination proc_log directory

script_path_windows <- paste(getwd(), "/process/r/", sep="")  

#script_path_windows <- gsub("/", "\\\\", script_path_windows)
script_path_windows  # filepath

# name of our script
script_name <- paste(theme, "_", subtheme,"_", location, "_s01_",version, ".R", sep="")  # s  (e.g. "_s01_") stands for sequential order incase you have divided the work in several scripts
script_name # file name
# ========================================================================================================= #
# INDATA
# ========================================================================================================= #
# Set "data_raw" sources ----------------------------------------------------------------------------------------
# Make sure all data_raw files you use have a corresponding metasym.txt file that look ok, code to help below. 
# use xxx.r script to summarize latest sources in data_raw according to your theme.. ED WE NEED TO INSERT LINK HERE

# SOURCE 1 data_raw:
source1_dir <- source_dir <-  "./path/uniqueid/"
source1_id <- source_id <- "uniqueid" # cut and past your source id here
source1_metasym <- read_tsv(paste(source_dir, source_id, "_metasym.txt", sep="")) # reads the associated metasym file
source1_metasym # check that metasym file is populated, and fix if not

# check what data is in data_raw directory and read in all files to use from the directory:
dir(source_dir, recursive=T) 

# File 1:
dir(source_dir) # check what data is in data_raw directory
source_file <- "internalpath/name.ext"  # your data file(s)
source1_path1 <- source_path <- paste(source_dir, source_file, sep="")
source_path
descriptive_name <- indata1 <- st_read(source_path)

# File 2: ...

# SOURCE 2 data_raw:
# file 1 ...
  
  
# Set "data" sources -------------------------------------------------------------------------------------------------
# Products in data are based on information in data_raw. Make sure they are all accompanied by sourcesym.txt to track source IDs

# DATA 1: standard_grid 1km
data1_dir <- data_dir <-  "./data/reg/grid/grid/v01/"  # input your data directory

dir(data_dir) # check what data is in this directory
data_file <- "grid_1km_v01.1.tif" # your data file
data1_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data1_sourcesym # check so source IDs exists and make sense
data1_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
grid_1km <- raster(data_path) # read in your data object

# DATA 2: standard_grid 1km na values
data2_dir <- data_dir <-  "./data/reg/grid/grid/v01/"  # input your data directory
dir(data_dir) 
data_file <- "grid_1km_na_v01.1.tif" # your data file
data2_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data2_sourcesym # check so source IDs exists and make sense
data2_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
grid_1km_na <- raster(data_path) # read in your data object 


# DATA 3: standard_grid 250m
data3_dir <- data_dir <-  "./data/reg/grid/grid/v01/"  # input your data directory
dir(data_dir) 
data_file <- "grid_250m_v01.1.tif" # your data file
data3_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data3_sourcesym # check so source IDs exists and make sense
data3_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
grid_250m <- raster(data_path) # read in your data object 

# DATA 4: standard_grid 250m NA values
data4_dir <- data_dir <-  "./data/reg/grid/grid/v01/"  # input your data directory
dir(data_dir) # check what data is in this directory
data_file <- "grid_250m_na_v01.1.tif" # your data file
data4_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data4_sourcesym # check so source IDs exists and make sense
data4_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
grid_250m_na <- raster(data_path) # read in your data object 
grid_250m_na


# DATA 5: bounding box shapefile
data5_dir <- data_dir <-  "./data/reg/grid/scope/v01/"  # input your data directory
dir(data_dir) # check what data is in this directory
data_file <- "wiosym_grid_bounding_box_v01.shp" # your data file
data5_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data5_sourcesym # check so source IDs exists and make sense
data5_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
grid_poly <- st_read(data_path) # read in your data object if appropriate and perhaps change name to more informative...
grid_poly
plot(grid_poly)

# DATA 6: keep adding as needed, or remove unwanted data above....


# Create sourcesym file for all sources --------------------------------------------------------------------
# adding all unique ID together for output sourcesym files

# Summary of data raw source directorys
(data_raw_sources <- tibble(id = c(source1_id, ...))) # REPLACE - keep adding all data_raw souces e.g. c( source1_id, source2_id, ...)


# Summary of data sourcesym files
data_sources <- data_raw_sources %>% 
  add_row(data1_sourcesym)%>%  # add any "sourcesym" files applicable to this product, keep adding rows with data sourcesym files as needed
  add_row(data2_sourcesym)%>% 
  add_row(data3_sourcesym)%>% 
  add_row(data4_sourcesym)%>% 
  add_row(data5_sourcesym)%>% 
  #add_row(data6_sourcesym)%>% 
  unique() %>% 
  print()

# Sources combined
(all_sources <- data_raw_sources %>% add_row(data_sources))

# if some output files have less sources, copy this section and create individual files below
# ========================================================================================================= #
#### PROCESS ####
# ========================================================================================================= #

# from here and on you are free to write your own process, annotate well and use the style guide ref above when applicable 

# example to help get started (delete if easier to write your own way..)
# name of first section ... -----------------------------------------------------------------------------------------------------
# well annotated code... annotate why you do something, less focus on what you do (the code will speak for that)

# brief description why
object1_name <- somefunction(somedata)

# write to file in work directory (if needed)
workfile1 <- workfile <- paste(work_path, subtheme, "_name_", version, d_version, ".tif", sep="")
workfile

writeRaster(object_name, workfile, overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')
object1_name <- raster(workfile)
object1_name

# name of second section ... --------------------------------------------------------------------------------

# ========================================================================================================= #
#### EXPORT ####
# ========================================================================================================= #
# Write products to root directory, accompanied by yourfilename_sourcesym.txt to track source IDs
# Temp files (../proc/) generated in "PROCESS Section will lack metadata (don't use from elsewhere). Save files of value to root with correct naming & *_sourcesym.txt. 
# Check outputs - If needed, check your outputs in /proc/ before writing to file (Make note of manual check elsewhere (e.g. qgis). perhaps integrate leaflet here?
# Note on looped products with same sources, make still sure to save *_sourcesym.txt for each individual product.

# SOURCESYM ---------------------------------------------------------------------------------------------------------------
# Create sourcesym file WITH all source IDs used (data and data_raw) 

# data raw sources (dir)
(data_raw_sources <- tibble(id = c(source1_id, ...))) # REPLACE - keep adding all data_raw souces e.g. c( source1_id, source2_id, ...)

# data sources (files)
data_sources <- data1_sourcesym %>% 
  add_row(data2_sourcesym)%>%  # add any "sourcesym" files applicable to this product (check indata section), keep adding sourcesym files as needed
  #add_row(data3_sourcesym)%>% 
  #add_row(data4_sourcesym)%>% 
  #add_row(data5_sourcesym)%>% 
  #add_row(data6_sourcesym)%>% 
  unique() %>% 
  print()

# Sources combined
(all_sources <- data_raw_sources %>% add_row(data_sources))

# if a product only use a subset of the total sources, copy the section above and create individual files below (e.g. product1_sources <- ...)

# PRODUCT 1:       -------------------------------------------------------------------------------------------------------
# read objects
product_orig <- REPLACE # input and check your object to be written to file
product_norm01 <- REPLACE # input and check your object to be written to file
product_uncertainty <- REPLACE # input and check your object to be written to file

# check products, if done in other places (e.g. qgis) write comment
plot(c(product_orig, product_norm01, product_uncertainty)) # use leaflet here instead?

# set names
scale <- "REPLACE" # "1km" / "250m"
unit <- "REPLACE"  # add SI unit for product with original values (use short abbreviation e.g. meter = m, kilometer = km, percent = perc) 
component_names <- read_tsv("./shiny_data_upload/component_names.txt") 
component_names_selected <- component_names %>%
  dplyr::filter(theme_folder == theme) %>% 
  dplyr::select(group, file_name, description) %>% 
  print(n = Inf)

descriptive_name <- "REPLACE"  # name of product, for pressure/ecosystem components select from name list above

# paths

(product_orig_path <- paste(dest_path,"/", descriptive_name, "_orig_", unit, "_", scale, "_", version, d_version, ".tif", sep="" ))
(product_norm01_path <- paste(dest_path,"/", descriptive_name, "_norm01_", scale, "_", version, d_version, ".tif", sep="" ))
(product_uncertainty_path <- paste(dest_path,"/", descriptive_name, "_uncertainty_", scale, "_", version, d_version, ".tif", sep="" ))


# Write to file
writeRaster(product_orig, product_orig_path, COMPRESS=LZW, datatype = raster_data_type, overwrite=TRUE)
writeRaster(product_norm01, product_norm01_path, COMPRESS=LZW, datatype = 'FLT4S', overwrite=TRUE)
writeRaster(product_uncertainty, product_uncertainty_path, COMPRESS=LZW, datatype = 'INT1U', overwrite=TRUE)
write_tsv(all_sources, paste(product_orig_path, "_sourcesym.txt", sep=""))  # writes unique sourcesym file for all products, change "all_sources" to "product1_sources" if subset only


# SAVE SCRIPT ------------------------------------------------------------------------------------------------------------
# 1. save your current R script (File/Save)
# 2. Save a carbon copy of the R script to proc_log (identical script to what is used to create the products exported
# Do not change the *.Rmd carbon copy in proc_log/ unless you change this version of the products (new versions are developed under process/r)


# 2. run code below and go to File/Save As and paste link and name:
dest_path_script <- paste(getwd(), proc_path)
print(dest_path_script <- gsub(" .", "", dest_path_script)) # path
print(script_name_copy <- paste(gsub(".Rmd", "", script_name), d_version, ".Rmd", sep="")) # script name


#==========================================================================================================================#
#### FINAL CHECK ####
#==========================================================================================================================#
# Developer check
# Do your own final check on the products and double check sourcesym file exists and are complete for all files in root 
# -> Sign at top of script in about section (# Approved by:)

# External check 
# To ensure repeatability and quality, make sure one colleage can run the script
# When script is proven repeatable/understandable and products/metadata look ok 
# -> Sign at top of script in about section (# Approved by:), and make any comments if needed
# There is also a work_log document in onedrive / gitHUB, check current arrangement with Swedish dev team

