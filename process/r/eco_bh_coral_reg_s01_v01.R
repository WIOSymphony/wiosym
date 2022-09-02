
#### ABOUT ####
# ========================================================================================================= #
# Purpose: Script to create WIOSym regional standard grid
#
# Brief description: 

# Suggestions on future imporvements: 
#
# Created: Initial, organisation, yymmdd
# Updated: Initial, organisation, yymmdd, comments

# Script and products checked by: (initial, organisation), yymmdd, 
  # comments: 

# ========================================================================================================= #

#### Note for developers ####
# ========================================================================================================= #    
# Style
  # Design the code for any process within WIOSym to be as human readable and repeatable as possibel.
  # This includes using tidyverse coding to the extent the programmer is comfortable with. An excellent style guide 
  # is found here https://style.tidyverse.org/index.html, use it to your best ability (but dont worry about PerfeCtion!) 

# Metadata
  # To track data flow in WIOSym its important to make sure all indata have a unique ID and a _metasym.txt file to go along
  # and that the are all logged in the final products, this ensures we know what sources where used.
  # Note that there is a difference how we ingest this information if the indata comes from "data" or from "data_raw", 
  # instructions found in "set destinations"

# Naming and folders
  # To standardise how we name products as well as create folders there are instructions under "set destinations" 
  # in short, we try to use predfined standard names for any main folders and final products we create. 
  # They are stored in a seperate textfile that we can update when need be. 
  # Naming conventions in "proc" folder is up to developer. see it as a skratch space blizzfully free of conventiones! 
  # Main version is synked with the overall development for WIOSym, 2021 we are working with V1, version controll might include git or simply V01.1 etc..
  # Layers that are ready to be used by other scripts/processes are stored under the root destination folder together with 
  # a nameofmyfile_datasym.txt files with all ids used. No other files here.
# ========================================================================================================= # 

# PREPERATIONS ----------------------------------------------------------------------------------------------
# List dependencis - any code or processes that needs to be run for this script to work
# Make sure the following processes are run first: 1. ..., 2. ..., 


# Read packages ---------------------------------------------------------------------------------------------
# list all packages that this script is dependent on, remove any redundant ones
x <- c("tidyverse", "sf", "raster", "rgdal", "fasterize", "labelled", "gdalUtils", "foreign")

#install.packages(x, dependencies=TRUE)
                                            
library(tidyverse)
library(sf)
library(raster)
library(rgdal)
library(fasterize)
library(labelled)
library(gdalUtils)
library(foreign) # write dbf 

# Set main destinations ----------------------------------------------------------------------------------------

# set version
version = "v01" # (1, 2, 3... relates to major releases for WIOSym, use v01.1, or v01.1.1 if you need additional version control use "detailed version" below)
d_version = ".0"  # detailed version (you can add another level e.g. 0.0  -> v01.0.0)

# set main location (i.e. regional for products for whole WIO area, or individual countries for local WIOSym products)
read_tsv("./shiny_data_upload/locations.txt")  # prints available locations to choose from and their abbreviations

location = "reg" # choose applicable 3 letter abbreviation from "location_val"

# Set names and working folder using standard themes, if available

folders_raw <- read_tsv("./shiny_data_upload/folders.txt")
folders_data <- read_tsv("./process/templates/folders_data.txt") # additional folders for data not applicable for data_raw 
folders <- folders_raw %>% 
  add_row(folders_data)

# select main theme
unique(folders["theme_folder"])

theme1 <- theme <- "eco" # input your main theme e.g. "eco" (choose from list in "theme_folder")

# select subtheme / component

folders %>%
  filter(theme_folder == theme) %>% 
  dplyr::select(subtheme_folder, subtheme_name) %>% 
  print(n = Inf)

subtheme1 <- subtheme <- "bh_coral" # input your subtheme based on the most appropriate option from the list. If your subtheme is missing, please add new subtheme to folders.txt first (and make sure to sync with shared location for WIOSym if applicable)

# Create destination path and folder

dest1_path <- dest_path <- paste("./data/", location, "/", theme, "/", subtheme, "/", version, "/", sep="")  # location for any final products that can be used elsehwere along with datasym.txt file which logs all source IDs. 
dest_path

if (!dir.exists(dest_path)){
  dir.create(dest_path, recursive = TRUE)
}


# work directory for any temporary files
work1_path <- work_path <- paste(dest_path, "proc/", sep="") # in this folder you can organise any files you use in the process using your prefferred structure. No "final products" here, they go under root (path_dest1)
work_path

if (!dir.exists(work_path)){
  dir.create(work_path, recursive = TRUE)
}

# Process log directory, for saving process logs and the exact script (or mxd etc) used to create the output files
proc1_path <- proc_path <- paste(dest_path, "proc_log/", sep="") # in this folder you can organise any files you use in the process using your prefferred structure. No "final products" here, they go under root (path_dest1)
proc_path

if (!dir.exists(proc_path)){
  dir.create(proc_path, recursive = TRUE)
}

# ... if multiple themes keep adding destinations here, and not futher down in script


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
script_name <- paste(theme1, "_", subtheme1,"_", location, "_s01_",version, ".R", sep="")  # s  (e.g. "_s01_") stands for sequential order incase you have divided the work in several scripts

script_name # file name


# Set "data_raw" sources ----------------------------------------------------------------------------------------

# data1: coral data from allen atlas, full cover of wio area (except somalia)
source1_dir <- source_dir <-  "./data_raw/reg/eco/bhab/allen/gk2104281131/"
source1_id <- source_id <- "gk2104281131" # cut and past your source id here
source1_metasym <- read_tsv(paste(source_dir, source_id, "_metasym.txt", sep="")) # reads the associated metasym file
source1_metasym # check that metasym file is populated, and fix if not
# check what data is in data_raw directory and read in all files to use from the directory:
dir(source_dir) 
#File1: Allen coral data
dir(source_dir) # check what data is in data_raw directory
source_file <- "Benthic-Map/benthic.gpkg" # your data file(s)
source1_path1 <- source_path <- paste(source_dir, source_file, sep="")
source_path
coral_allen <- st_read(source_path)
#coral_sf_full


# data2: wcmc
source2_dir <- source_dir <-  "./data_raw/glo/eco/bh_coral/wcmc/gk2105272129/"
source2_id <- source_id <- "gk2105272129" # cut and past your source id here
source2_metasym <- read_tsv(paste(source_dir, source_id, "_metasym.txt", sep="")) # reads the associated metasym file
source2_metasym # check that metasym file is populated, and fix if not
# check what data is in data_raw directory and read in all files to use from the directory:
dir(source_dir, recursive=T) 
#File1: wcmc coral data
source_file <- "01_Data/WCMC008_CoralReef2018_Py_v4_1.shp"  # your data file(s)
source2_path1 <- source_path <- paste(source_dir, source_file, sep="")
source_path
coral_wcmc <- st_read(source_path)
coral_wcmc



# Set "data" sources -------------------------------------------------------------------------------------------------
# declare paths to all sources from "data" directory which are products created from information in data_raw. make sure that all sources are accompanied by the a datasym.txt file, if not fix...

# data1: standard_grid 1km
data1_dir <- data_dir <-  "./data/reg/grid/grid/v01/"  # input your data directory

dir(data_dir) # check what data is in this directory
data_file <- "grid_1km_v01.1.tif" # your data file
data1_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data1_sourcesym # check so source IDs exists and make sense
data1_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
grid_1km <- raster(data_path) # read in your data object if appropriate and perhaps change name to more informative...

# data2: standard_grid 1km na values
data2_dir <- data_dir <-  "./data/reg/grid/grid/v01/"  # input your data directory
dir(data_dir) 
data_file <- "grid_1km_na_v01.1.tif" # your data file
data2_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data2_sourcesym # check so source IDs exists and make sense
data2_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
grid_1km_na <- raster(data_path) # read in your data object if appropriate and perhaps change name to more informative...


# data3: standard_grid 250m
data3_dir <- data_dir <-  "./data/reg/grid/grid/v01/"  # input your data directory
dir(data_dir) 
data_file <- "grid_250m_v01.1.tif" # your data file
data3_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data3_sourcesym # check so source IDs exists and make sense
data3_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
grid_250m <- raster(data_path) # read in your data object if appropriate and perhaps change name to more informative...

# data4: standard_grid 250m NA values
data4_dir <- data_dir <-  "./data/reg/grid/grid/v01/"  # input your data directory
dir(data_dir) # check what data is in this directory
data_file <- "grid_250m_na_v01.1.tif" # your data file
data4_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data4_sourcesym # check so source IDs exists and make sense
data4_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
grid_250m_na <- raster(data_path) # read in your data object if appropriate and perhaps change name to more informative...
grid_250m_na


# data5: bounding box shapefile
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


# data6: mangrove
data6_dir <- data_dir <-  "./data/reg/eco/ch_mangrove/v01/"  # input your data directory
dir(data_dir) # check what data is in this directory
data_file <- "eco_ch_mangrove_1km_combined_presence_v01.0.tif" # your data file
data6_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data6_sourcesym # check so source IDs exists and make sense
data6_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
Mangrove_1km <- raster(data_path) # read in your data object if appropriate and perhaps change name to more informative...
Mangrove_1km

# data7: depth
data7_dir <- data_dir <-  "./data/reg/env/topo/v01/"  # input your data directory
dir(data_dir) # check what data is in this directory
data_file <- "env_topo_depth_mean_1km_v01.0.tif" # your data file
data7_sourcesym <- read_tsv(paste(data_dir, data_file, "_sourcesym.txt", sep=""))
data7_sourcesym # check so source IDs exists and make sense
data7_path <- data_path <- paste(data_dir, data_file, sep="")
data_path
depth_1km <- raster(data_path) # read in your data object if appropriate and perhaps change name to more informative...
depth_1km


# combined sourcesym file ------------------------------------------------------------------------------------------------------------------
# create sourcesym file to track all sources used in each output file (manualy list all "source_id" relevant objects, make sure to double check all your IDs are in there)
data_raw_sources <- tibble(id = c(source1_id, source2_id))

file_sourcesym <- data_raw_sources %>% 
    add_row(data1_sourcesym)%>%  # add any "sourcesym" files applicable to this product, keep adding rows with data sourcesym files as needed
    add_row(data2_sourcesym)%>%
    add_row(data3_sourcesym)%>% 
    add_row(data4_sourcesym)%>% 
    add_row(data5_sourcesym)%>% 
    #add_row(data6_sourcesym)%>% 
    #add_row(data7_sourcesym)%>% 
    unique() %>% 
    print()

file_sourcesym # double check that that all uniqueIDs for our data product sources are included 

# sourcesym for uncertainty layer
file_sourcesym2 <- data_raw_sources %>% 
  add_row(data1_sourcesym)%>%  # add any "sourcesym" files applicable to this product, keep adding rows with data sourcesym files as needed
  add_row(data2_sourcesym)%>%
  add_row(data3_sourcesym)%>% 
  add_row(data4_sourcesym)%>% 
  add_row(data5_sourcesym)%>% 
  add_row(data6_sourcesym)%>% 
  add_row(data7_sourcesym)%>% 
  unique() %>% 
  print()

file_sourcesym2



# PROCESS ----------------------------------------------------------------------------------------------------------------

# example -----------------------------------------------------------------------------
# some processing that is briefly described here

# brief description why
object1_name <- somefunction(somedata)

# write to file in work directory (if needed)
workfile1 <- workfile <- paste(work1_path, subtheme1, "_name_", version, d_version, ".tif", sep="")
workfile

writeRaster(object_name, workfile, overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')
object1_name <- raster(workfile)
object1_name


# Map Allen coral area (Allen coral atlas) from shape to raster grid --------------------------------------------------------------------

coral_sf <- coral_allen
glimpse(coral_sf)
unique(coral_sf$class)
class(coral_sf)

# Select coral polygons

coral_sf_sel1 <- coral_sf %>%
  filter(class=="Coral/Algae") %>% 
  #  mutate(class_code=class) %>% 
  mutate(class_code = recode(class, "Coral/Algae" = "1")) %>% 
  mutate(class_code = as.numeric(class_code))
#  recode(class_code, "Coral/Algae" = 1)

glimpse(coral_sf_sel1)


# convertion from coral shape to raster 

#write coral file to work directory for gdalutil
workfile1 <- workfile <- paste(work1_path, "coral_r_sel1_full.shp", sep="")
workfile

st_write(coral_sf_sel1,workfile) 

# write empty grid for gdalutil work
workfile2 <- workfile <- paste(work1_path, "grid_1km_na_coral_allen_full.tif", sep= "")
workfile
writeRaster(grid_1km_na, workfile,  overwrite=T, COMPRESS=LZW)
workfile3 <- workfile <- paste(work1_path, "grid_250m_na_coral_allen_full.tif", sep= "")
workfile
writeRaster(grid_250m_na, workfile,  overwrite=T, COMPRESS=LZW)


# 1km grid mapping
coral_warp <- gdalUtils::gdal_rasterize(src_datasource = workfile1,
                                        dst_filename = workfile2,
                                        b = 1,
                                        at = T,
                                        a = "class_code",
                                        output_Raster = TRUE,
)


# 250m grid mapping
coral_warp <- gdalUtils::gdal_rasterize(src_datasource = workfile1,
                                        dst_filename = workfile3,
                                        b = 1,
                                        at = T,
                                        a = "class_code",
                                        output_Raster = TRUE,
)



# Map all habitat in Allen coral atlas from shape to raster grid --------------------------------------------------------------------

coral_sf <- coral_allen
glimpse(coral_sf)
unique(coral_sf$class)
class(coral_sf)

#Select only coral polygons

coral_sf_sel2 <- coral_sf %>%
  filter(class!="Unknown") %>% 
  #  mutate(class_code=class) %>% 
  mutate(class_code = recode(class, "Coral/Algae" = "0", "Rubble" = "0", "Microalgal Mats" = "0", "Seagrass" = "0", "Rock"="0", "Sand"="0")) %>% 
  mutate(class_code = as.numeric(class_code))
#  recode(class_code, "Coral/Algae" = 1)

glimpse(coral_sf_sel2)


# convertion from coral shape to raster 

#write habitat file to work directory for gdalutil
workfile4 <- workfile <- paste(work1_path, "coral_r_sel2_full.shp", sep= "")

st_write(coral_sf_sel2,workfile, delete_layer = TRUE) 

# write empty grid for gdalutil work
workfile5 <- workfile <- paste(work1_path, "grid_1km_na_all_habitat_allen_full.tif", sep= "")
writeRaster(grid_1km_na, workfile,  overwrite=T, COMPRESS=LZW)
workfile6 <- workfile <- paste(work1_path, "grid_250m_na_all_habitat_allen_full.tif", sep= "")
writeRaster(grid_250m_na, workfile,  overwrite=T, COMPRESS=LZW)


# 1km grid mapping
#coral_warp <- gdalUtils::gdal_rasterize(src_datasource = paste(path_work, "coral_r_sel2_full.shp", sep= ""),
#                                        dst_filename = paste(path_work, "grid_1km_na_all_habitat_allen_full.tif", sep= ""),
#                                        b = 1,
#                                        at = T,
#                                        a = "class_code",
#                                        output_Raster = TRUE,
#)


# 250m grid mapping
coral_warp <- gdalUtils::gdal_rasterize(src_datasource = workfile4,
                                        dst_filename = workfile6,
                                        b = 1,
                                        at = T,
                                        a = "class_code",
                                        output_Raster = TRUE,
)

coral_warp
summary(coral_warp)




# Map wcmc coral area (global data) from shape to raster grid --------------------------------------------------------------------


coral_sf <- coral_wcmc
glimpse(coral_sf)
unique(coral_sf$METADATA_I)
class(coral_sf)

coral_sf_sel1 <- coral_sf %>%
  #  filter(class=="Coral/Algae") %>% 
  #  mutate(class_code=class) %>% 
  #  mutate(class_code = recode(class, "Coral/Algae" = "1")) %>% 
  mutate(class_code = 1)
#  recode(class_code, "Coral/Algae" = 1)

glimpse(coral_sf_sel1)


# convertion from coral shape to raster 

#write coral file to work directory for gdalutil
workfile7 <- workfile <- paste(work1_path, "coral_sf_wcmc_sel1_full.shp", sep= "")
st_write(coral_sf_sel1, workfile)


# write empty grid for gdalutil work
workfile8 <- workfile <- paste(work1_path, "grid_1km_na_coral_wcmc_full.tif", sep= "")
writeRaster(grid_1km_na, workfile,  overwrite=T, COMPRESS=LZW)
workfile9 <- workfile <- paste(work1_path, "grid_250m_na_coral_wcmc_full.tif", sep= "")
writeRaster(grid_250m_na, workfile,  overwrite=T, COMPRESS=LZW)


# 1km grid mapping
coral_warp <- gdalUtils::gdal_rasterize(src_datasource = workfile7,
                                        dst_filename = workfile8,
                                        b = 1,
                                        at = T,
                                        a = "class_code",
                                        output_Raster = TRUE,
)


# 250m grid mapping
coral_warp <- gdalUtils::gdal_rasterize(src_datasource = workfile7,
                                        dst_filename = workfile9,
                                        b = 1,
                                        at = T,
                                        a = "class_code",
                                        output_Raster = TRUE,
)


# Combine all corals into one layer -------------------------------------------------------------------------------------
# paths to input rasters
path_work <- work1_path  # adjusting to fit old code...
wcmc_coral_250m <- paste(path_work, "grid_250m_na_coral_wcmc_full.tif", sep= "")

allen_coral_250m <- paste(path_work, "grid_250m_na_coral_allen_full.tif", sep= "")

allen_habitat_250m <- paste(path_work, "grid_250m_na_all_habitat_allen_full.tif", sep= "")

#allen_habitat_250m <- raster(paste(path_work, "grid_250m_na_all_habitat_allen_full.tif", sep= ""))
#summary(allen_habitat_250m)



# combine layers to see overlapping areas
wcmc_coral_250m_r <- raster(paste(path_work, "grid_250m_na_coral_wcmc_full.tif", sep= ""))

allen_coral_250m_r <- raster(paste(path_work, "grid_250m_na_coral_allen_full.tif", sep= ""))

combined_r <- allen_coral_250m_r + wcmc_coral_250m_r +wcmc_coral_250m_r
writeRaster(combined_r, paste(path_work, "grid_250m_na_coral_combined_overlap_1-3.tif", sep= ""), overwrite=T, COMPRESS=LZW)


# reclassify wcmc cells which are not mapped as coral in Allen atlas, due to lower spatial precision in wcmc data


# write empty grid for gdalutil output

writeRaster(grid_250m_na, paste(path_work, "grid_250m_na_coral_combined_2.tif", sep= ""), overwrite=T, COMPRESS=LZW)

coral_250m_combined <- paste(path_work, "grid_250m_na_coral_combined_2.tif", sep= "")

coral_250m_combined <- gdalUtils::mosaic_rasters(
  gdalfile = c(wcmc_coral_250m, allen_habitat_250m, allen_coral_250m),
  dst_dataset = coral_250m_combined,
  output.vrt = NULL,
  output_Raster = FALSE,
  separate = FALSE,
  trim_margins = NULL,
  gdalwarp_index = 1,
  gdalwarp_params = list(r = "near"),
  #force_ot = "Int16",
  verbose = FALSE,
  COMPRESS=LZW,
)

#-co COMPRESS=DEFLATE
coral_250m_combined <- raster(paste(path_work, "grid_250m_na_coral_combined_2.tif", sep= ""))
writeRaster(coral_250m_combined, paste(path_work, "grid_250m_na_coral_combined_2_lzw.tif", sep= ""), overwrite=T, COMPRESS=LZW, datatype = 'INT4S')

combined_coral_250m_agg1km <- aggregate(coral_250m_combined, fact=4, fun=sum) / 16 #4*4=16, i.e. divide by the maximum cell number to get proportions
writeRaster(combined_coral_250m_agg1km, paste(path_work, "grid_1km_na_coral_combined_2_proportion250m.tif", sep= ""),  overwrite=T, COMPRESS=LZW)

combined_coral_250m_agg1km_binary <- aggregate(coral_250m_combined, fact=4, fun=max) 
writeRaster(combined_coral_250m_agg1km_binary, paste(path_work, "grid_1km_na_coral_combined_2.tif", sep= ""),  overwrite=T, COMPRESS=LZW)


# alternative route, aggregate the two coral layers individually
#allen_coral_250m_agg1km <- aggregate(allen_coral_250m, fact=4, fun=sum) / 16 #4*4=16, i.e. divide by the maximum cell number to get proportions
#writeRaster(allen_coral_250m_agg1km, paste(path_work, "grid_1km_na_coral_allen_proportion250m.tif", sep= ""),  overwrite=T, COMPRESS=LZW)

#wcmc_coral_250m_agg1km <- aggregate(wcmc_coral_250m, fact=4, fun=sum) / 16 #4*4=16, i.e. divide by the maximum cell number to get proportions
#writeRaster(wcmc_coral_250m_agg1km, paste(path_work, "grid_1km_na_coral_wcmc_proportion250m.tif", sep= ""),  overwrite=T, COMPRESS=LZW)




# Map Uncertainty ----------------------------------------------------------------------------------------------------

#steps:
# 1. identify well mapped areas without corals - (coral atlas without NA) as well as mangroves could be a start

# mapped_habitat <-  
coral_habitat_1km <- raster(paste(path_work, "grid_1km_na_coral_combined_2.tif", sep= "")) #(0 = mapped no coral, 1= coral)
coral_habitat_1km

coral_habitat_1km_wcmc <- raster(paste(path_work, "grid_1km_na_coral_wcmc_full.tif", sep= ""))
coral_habitat_1km_wcmc[coral_habitat_1km_wcmc <2] <- 1
coral_habitat_1km_allen <- raster(paste(path_work, "grid_1km_na_all_habitat_allen_full.tif", sep= ""))
coral_habitat_1km_allen[coral_habitat_1km_allen <2] <- 1

Mangrove_1km

# Potential habitat / no potential habitat
depth_1km

reclass_df <- c(-Inf, -200, 0,
                -200, Inf, 1)

reclass_df
#Reorder into matrix
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)

reclass_m

depth_1km_potential <- reclassify(depth_1km,
                                  reclass_m, include.lowest=TRUE)



# identify max min latitude for corals using combined dataset --------------------------------------------------
r1 <- coral_habitat_1km
r2 <- depth_1km_potential # raster to limit extent on
r1NaM <- is.na(as.matrix(r1))
# Find the columns and rows that are not completely filled by NAs
colNotNA <- which(colSums(r1NaM) != nrow(r1))
rowNotNA <- which(rowSums(r1NaM) != ncol(r1))
# Find the extent of the new raster by using the first ans last columns and rows that are not completely filled by NAs. Use crop() to crop the new raster.
# r3Extent <- extent(r1, rowNotNA[1], rowNotNA[length(rowNotNA)],
#                    colNotNA[1], colNotNA[length(colNotNA)])

# find min/max latitude with corals in data and use to drop potential habitat
r3Extent <- extent(r1, r1=rowNotNA[1], r2=rowNotNA[length(rowNotNA)])



r3 <- crop(r2, r3Extent)
# Plot the rasters for comparison.
layout(matrix(1:2, nrow=1))
plot(r1)
plot(r3)

potential_habitat <- raster::extend(r3, extent(r1))
grid_1km[grid_1km==1] <- 0
potential_habitat <- mosaic(potential_habitat, grid_1km, fun=max)
plot(potential_habitat)
#NAvalue(potential_habitat) <- 3

writeRaster(potential_habitat, paste(path_work, "coral_potential_habitat_v00.tif", sep= ""),  overwrite=T, COMPRESS=LZW)

potential_habiat <- raster(paste(path_work, "coral_potential_habitat_v00.tif", sep= ""))

r_1 <- coral_habitat_1km_allen
r_2 <- mosaic(coral_habitat_1km_wcmc, grid_1km, fun=max)
r_3 <- mosaic(Mangrove_1km, grid_1km, fun=max)
r_4 <- potential_habitat
#r_5 <- grid_1km
#r_5 <- r_5[r_5==0] <- 5


# uncertainty classification---------------------------------------------------

coral_habitat_1km_allen <- raster(paste(path_work, "grid_1km_na_all_habitat_allen_full.tif", sep= ""))
coral_habitat_1km_allen[coral_habitat_1km_allen <2] <- 2
coral_habitat_1km_allen

coral_habitat_1km_wcmc <- raster(paste(path_work, "grid_1km_na_coral_wcmc_full.tif", sep= ""))
coral_habitat_1km_wcmc[coral_habitat_1km_wcmc <2] <- 3
coral_habitat_1km_wcmc

Mangrove_1km
Mangrove_1km[Mangrove_1km < 2] <- 2


potential_habitat[potential_habitat == 1] <- 5
potential_habitat

r_mosaic <- merge(coral_habitat_1km_allen, coral_habitat_1km_wcmc, Mangrove_1km, potential_habitat)

# calculate raster statistics and ratify -----------------------------------------------------
df <- as.data.frame(r_mosaic)
df <- as_tibble(df)   %>%  
  mutate_all(as.integer) %>% 
  drop_na()
resolution = 1000
stat <- df %>% 
  group_by(layer) %>% 
  summarize(n = n()) %>%
  mutate(km2 = round(n*resolution*resolution/1000000, 2)) %>% 
  mutate(percent = round(100 * n/sum(n), 2)) %>% 
  print()


r_mosaic <- ratify(r_mosaic)
r_mosaic_rat <- levels(r_mosaic)[[1]]

value <- c(0, 1, 2, 3, 4, 5)
class_name <- c("outside range", "confirmed presence", "very good model", "good model", "poor model", "no data")
count <- as.factor(c(stat[["n"]][0], stat[["n"]][1],	stat[["n"]][2],	stat[["n"]][3],	stat[["n"]][4],	stat[["n"]][5]))
class_km2 <- as.factor(c(stat[["km2"]][0], stat[["km2"]][1],	stat[["km2"]][2],	stat[["km2"]][3],	stat[["km2"]][4],	stat[["km2"]][5]))

stat <- stat %>%
  rename(value = Value)

rat_table <- data.frame(value, class_name) %>%
  as_tibble() %>% 
  mutate(value = as.integer(value))%>% 
  left_join(stat, by="value")
rat_table

rat_table <- rat_table %>% 
  filter(n>=0) %>% 
  print()

r_mosaic_rat$Value <- rat_table$value
r_mosaic_rat$class_name <- rat_table$class_name
r_mosaic_rat$count <- rat_table$n
r_mosaic_rat$class_km2 <- rat_table$km2
r_mosaic_rat$class_perc <- rat_table$percent
levels(r_mosaic) <- r_mosaic_rat
r_mosaic_rat

write.dbf(r_mosaic_rat, file = paste(path_work, "grid_1km_coral_combined2_uncertainty_", version, "m.tif.vat.dbf", sep=''), factor2char = TRUE, max_nchar = 254)

writeRaster(r_mosaic, filename = paste(path_work, "grid_1km_coral_combined2_uncertainty_", version, ".tif", sep= ""), overwrite=TRUE, COMPRESS=LZW, datatype = 'INT2S')




#suggested cathegories
#outside_range <- 0 # outside range/habitat of ecosystem component
#confirmed_presence <- 1
#very_good_model <- 2 #remote sensing assisted
#good_model <- 3
#poor_model <- 4
#no_data <- 5






# OUTPUTS -----------------------------------------------------------------------------------------
# OBS only write final files to the root directory using the exampels from here and on, any files generated in the 
# "PROCESS Section will lake proper metadata (i.e. ..._sourcesym.txt) should be written only to ../proc/

# Check outputs  --------------------------------------------------------------------------------
# If needed, check your outputs in /proc/ before writing to file. In this case a manual check in ArcGIS was done.

# check sourcesym
file_sourcesym
# check sourcesym uncertainty
file_sourcesym2

# write to product file: product 1  -------------------------------------------------------------------------------------- 
product1_path <- product_path <- paste(dest1_path, theme1, "_", subtheme1, "_reef_presence", "_1km_", version, d_version, ".tif", sep="" )
product_path



# write to file 
write_tsv(file_sourcesym, paste(product_path, "_sourcesym.txt", sep=""))  # write unique sourcesym files for all products
product <- combined_coral_250m_agg1km_binary # input and check your object to be written to file
writeRaster(product, product_path, COMPRESS=LZW, datatype = 'INT4S', overwrite=TRUE)


# write to product file: product 2 1km reef proportion -------------------------------------------------------------------------------------- 

product2_path <- product_path <- paste(dest1_path, theme1, "_", subtheme1, "_reef_proportion", "_1km_", version, d_version, ".tif", sep="" )
product_path

# write to file 
write_tsv(file_sourcesym, paste(product_path, "_sourcesym.txt", sep=""))  # write unique sourcesym files for all products
product <- combined_coral_250m_agg1km # input and check your object to be written to file
writeRaster(product, product_path, COMPRESS=LZW, datatype = 'INT4S', overwrite=TRUE)


# write to product file: product 3 250m reef combined-------------------------------------------------------------------------------------- 

product3_path <- product_path <- paste(dest1_path, theme1, "_", subtheme1, "_reef_presence", "_250m_", version, d_version, ".tif", sep="" )
product_path

# write to file 
write_tsv(file_sourcesym, paste(product_path, "_sourcesym.txt", sep=""))  # write unique sourcesym files for all products
product <- coral_250m_combined # input and check your object to be written to file
writeRaster(product, product_path, COMPRESS=LZW, datatype = 'INT4S', overwrite=TRUE)


# write to product file: product 4 uncertainty 1km  -------------------------------------------------------------------------------------- 

product4_path <- product_path <- paste(dest1_path, theme1, "_", subtheme1, "_reef_uncertainty", "_1km_", version, d_version, ".tif", sep="" )
product4_path_rat <- product_path_rat <- paste(dest1_path, theme1, "_", subtheme1, "_reef_uncertainty", "_1km_", version, d_version, ".rat", sep="" )

product_path

file_sourcesym2

# write to file 
write_tsv(file_sourcesym2, paste(product_path, "_sourcesym.txt", sep=""))  # write unique sourcesym files for all products
product <- r_mosaic # input and check your object to be written to file
product_rat <- r_mosaic_rat
writeRaster(product, product_path, COMPRESS=LZW, datatype = 'INT2S', overwrite=TRUE)
write.dbf(product_rat, product_path_rat, factor2char = TRUE, max_nchar = 254)

# Save copy of R scripts ---------------------------------------------------------------------------------------------------

# 1. save your current R script (File/Save)

# 2. Save a copy of the R script here (if you update the files, make sure to redo the save):

#   run code below and go to File/Save As and paste link and name:

dest1_path_script <- paste(getwd(), proc1_path)
dest1_path_script <- gsub(" .", "", dest1_path_script) 

#dest1_path_script <- gsub("/", "\\\\", dest1_path_script) # code needs work...

dest1_path_script # path

script_name_copy <- paste(gsub(".R", "", script_name), d_version, ".R", sep="")
script_name_copy # name


# OBS! If you make more edits after this "final" save make sure to go back to your original R script under /process/r again, then resave this one in the end


# FINAL CHECK  -------------------------------------------------------------------------------------------
# Do your own final check on the products and doublecheck sourcesym file exists and are complete for all files in root
# checked by: initial, org, date

# External check -------------------------------------------------------------------------------------
# To ensure repeatability and quality, make sure atleast one more colleage can run the script and make any final edits.
# When script is proven repeatable/understandable and products/metadata look ok 
# -> Sign at top of script in about section (# Approved by:), and make any comments if needed





