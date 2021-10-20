# ========================================================================================================= #
#### ABOUT ####
# ========================================================================================================= #
# Brief description: 
# Script by: Initial, org
# Updated: Initial, org
# Checked by: Initial, org

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
folders_data <- read_tsv("./process/templates/folders_data.txt")
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

# IMPORTANT only write final files to the root directory accompanied by yourfilename_sourcesym.txt to track all IDs used
# temp files (../proc/) from "PROCESS Section will lack metadata, don't use from elsewhere. Any file of value resave to root with correct naming & *_sourcesym.txt. 

# Check outputs 
# If needed, check your outputs in /proc/ before writing to file. In this case a manual check in ArcGIS was done.

# PRODUCT 1 to file:   ------------------------------------------------------------------------------------ 
# if multiple products are produced in a loop with same sources they can be written to file in one chunk (i.e. this section), 
# just make sure to save one sourcesym file for each individual produc. 

# INPUT
scale <- "REPLACE" # "1km" / "250m"
descriptive_name <- "REPLACE"  # "name_of_product" for pressure/ecosystem components used standard names in list

# original values product
unit_name <- "orig"  # add "orig_unit" when appropriate (e.g. orig_dB)
raster_data_type <- 'FLT4S' # If integer use 'INT1U' or 'INT4S' to save space

product_path <- paste(dest_path, theme, "_", subtheme, "_", descriptive_name, "_", "unit_name", "_", scale, "_", version, d_version, ".tif", sep="" )
product_path
product <- REPLACE # input and check your object to be written to file
writeRaster(product, product_path, COMPRESS=LZW, datatype = raster_data_type, overwrite=TRUE)

# normalised value product
unit_name <- "norm01"
raster_data_type <- 'FLT4S'

product_path <- paste(dest_path, theme, "_", subtheme, "_", descriptive_name, "_", "unit_name", "_", scale, "_", version, d_version, ".tif", sep="" )
product_path
product <- REPLACE # input and check your object to be written to file
writeRaster(product, product_path, COMPRESS=LZW, datatype = raster_data_type, overwrite=TRUE)

# uncertainty value product
unit_name <- "norm01"
raster_data_type <- 'INT1U'

product_path <- paste(dest_path, theme, "_", subtheme, "_", descriptive_name, "_", "unit_name", "_", scale, "_", version, d_version, ".tif", sep="" )
product_path
product <- REPLACE # input and check your object to be written to file
writeRaster(product, product_path, COMPRESS=LZW, datatype = raster_data_type, overwrite=TRUE)


# Sourcesym file ------------------------------------------------------------------------------------
# help track all sources used in each product

# INPUT - manually list all data_raw source_IDs
file1_data_raw_sources <- file_data_raw_sources <- tibble(id = c("uniqueid1", "uniqueid2", "uniqueid3"))

# INPUT - manually list all sourcesym files related to this product, keep adding rows as needed
file1_sourcesym <- file_sourcesym <- file_data_raw_sources %>% 
#  add_row(data1_sourcesym)%>%
#  add_row(data2_sourcesym)%>% 
  unique() %>% 
  print()

file_sourcesym # double check that that all uniqueIDs for our data product sources are included 

# write to file 
write_tsv(file_sourcesym, paste(product_path, "_sourcesym.txt", sep=""))  # write unique sourcesym files for all products

# PRODUCT 2 to file:   -------------------------------------------------------------------------------- 
#copy from above..


# Sourcesym file: if sources same as previous products use this.., if other sources update using above code. 
file_sourcesym 
# write to file 
write_tsv(file_sourcesym, paste(product_path, "_sourcesym.txt", sep=""))  # write unique sourcesym files for all products


# Save copy of R scripts ---------------------------------------------------------------------------------------------------

# 1. save your current R script (File/Save)

# 2. Save a copy of the R script here (if you update the files, make sure to redo the save):

#   run code below and go to File/Save As and paste link and name:

dest_path_script <- paste(getwd(), proc_path)
dest_path_script <- gsub(" .", "", dest_path_script) 

#dest_path_script <- gsub("/", "\\\\", dest_path_script) # code needs work...

dest_path_script # path

script_name_copy <- paste(gsub(".R", "", script_name), d_version, ".R", sep="")
script_name_copy # name


# OBS! If you make more edits after this "final" save make sure to go back to your original R script under /process/r again, then re save this one when final product is complete


#### FINAL CHECK ####
# ========================================================================================================= #
# Do your own final check on the products and double check sourcesym file exists and are complete for all files in root
# checked by: initial, org, date

# External check -------------------------------------------------------------------------------------
# To ensure repeatability and quality, make sure one colleage can run the script
# When script is proven repeatable/understandable and products/metadata look ok 
# -> Sign at top of script in about section (# Approved by:), and make any comments if needed
# There is also a work_log document in onedrive / gitHUB, check current arrangement with Swedish dev team





