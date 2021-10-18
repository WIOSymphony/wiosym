
#### ABOUT ####
# ========================================================================================================= #
# Purpose: Script to create WIOSym depth grid from gebco data
#
# Brief description: standardised depth data to 250m and 1km grid 

# Suggestions on future improvements: 
#
# Created: GK, sgu, 210528
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
x <- c("tidyverse", "sf", "raster", "rasterVis", "rgdal") #, "fasterize", "labelled", "gdalUtils")

install.packages(x, dependencies=TRUE)


library(tidyverse)
library(sf)          # working with sf objects in R
library(raster)      # working with rasters in R
library(rasterVis)   # raster visualization
library(rgdal)       # links to the GDAL library

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

theme1 <- theme <- "act" # input your main theme e.g. "eco" (choose from list in "theme_folder")

# select subtheme / component

folders %>%
  filter(theme_folder == theme) %>% 
  dplyr::select(subtheme_folder, subtheme_name) %>% 
  print(n = Inf)

subtheme1 <- subtheme <- "ship" # input your subtheme based on the most appropriate option from the list. If your subtheme is missing, please add new subtheme to folders.txt first (and make sure to sync with shared location for WIOSym if applicable)

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
# locate available mangrove data

# data_raw1: description
source1_dir <- source_dir <-  "./data_raw/glo/act/ship/lit/gk2105311532/"
source1_id <- source_id <- "gk2105311532" # cut and past your source id here
source1_metasym <- read_tsv(paste(source_dir, source_id, "_metasym.txt", sep="")) # reads the associated metasym file
source1_metasym # check that metasym file is populated, and fix if not
# check what data is in data_raw directory and read in all files to use from the directory:
dir(source_dir, recursive=T) 
#File1: description
dir(source_dir) # check what data is in data_raw directory
source_file <- "L1_World_Grids_History_Yearly_10_Min_2015_public version.csv"  # your data file(s)
source1_path1 <- source_path <- paste(source_dir, source_file, sep="")
source_path

ais_data <- indata1 <- read_csv(source_path)

ais_data_v00_1 <- read_csv("./_archive/version0/data_raw/activity/shipping/reg/linwu/20200317/L1_World_Grids_History_Yearly_10_Min_2015_EastAfrica_Categories.csv")
ais_data_v00_2 <- read_csv("./_archive/version0/data_raw/activity/shipping/reg/linwu/20200324/L1_World_Grids_History_Yearly_10_Min_2015——EastAfrica.csv")

# data_raw2: description
# file 1 ...
  
  
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

# create sourcesym file for all sources -------------------------------------------------------------
# adding all unique ID together for output sourcesym files
data_raw_sources <- tibble(id = c(source1_id))


file_sourcesym <- data_raw_sources %>% 
  add_row(data1_sourcesym)%>%  # add any "sourcesym" files applicable to this product, keep adding rows with data sourcesym files as needed
  add_row(data2_sourcesym)%>% 
  add_row(data3_sourcesym)%>% 
  add_row(data4_sourcesym)%>% 
  add_row(data5_sourcesym)%>% 
  unique() %>% 
  print()

# if some output files have less sources, copy this chunk and create indiviual files

# PROCESS ----------------------------------------------------------------------------------------------------------------

# convert old script names to new template --------------------------------------------------
path_data <- work1_path

# compute lat long center coordinate ----------------------------------------------------------------------
# from csv file cell boundaries (divided by 600 000 for unknown reasons...) 

ais_data_tbl <- ais_data_v00_2 %>%
  mutate(long=((.[[1]] + ( .[[3]]- .[[1]])))/600000) %>%
  mutate(lat=((.[[2]] - ( .[[2]] - .[[4]])))/600000) %>%
  dplyr::select(lat, long, everything()) %>%
  print()


# convert to raster and stack ---------------------------------------------------------------------

# main themes to test... all ais time, all ais crossing
# "CrossingTimeSum" 
r1 <- rasterFromXYZ(ais_data_tbl[, c(2,1,7)])
crs(r1) <- crs("+proj=longlat +datum=WGS84")
r1
plot(r1)
# "Vessel_Count"
r2 <- rasterFromXYZ(ais_data_tbl[, c(2,1,9)])
crs(r2) <- crs("+proj=longlat +datum=WGS84")
r2
plot(r2)


r <- stack(rasterFromXYZ(ais_data_tbl[, c(2,1,7)]))
for (i in 8:length(ais_data_tbl)){
  r_stack <- rasterFromXYZ(ais_data_tbl[, c(2,1,i)])
  r <- addLayer(r, r_stack)
  #rname <- names(ais_data[i])
  #writeRaster(r, paste(path_data, "ais_worldgrid_10min_2015_global_", rname, ".tif", sep= ""),  overwrite=T, COMPRESS=LZW)
}

#r_path <- paste(path_data, "ais_worldgrid_10min_2015_global_", "stack", ".tif", sep= "")
#writeRaster(r, r_path,  overwrite=T, COMPRESS=LZW)
#r_stacked <- raster(r_path)
#grid_1km
#r_res <- resample(r, grid_1km, method='bilinear')

#dir(path_data)

#r <- raster(paste(work1_path,"ais_worldgrid_10min_2015_global_Vessel_Count.tif", sep=""))
#crs(r) <- crs("+proj=longlat +datum=WGS84")
#plot(r)

# mapping raster to grid
#projectRaster(from = r1, to = r2, alignOnly = TRUE) 


# interpolate rasters ----------------------------------------------------------------------------------
# interpolate using focal statistcis

#grid bounding box in lat/long
grid_na <- raster(xmn=6, xmx=80, ymn=-42, ymx=18, 
                     crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", resolution = res(r1), vals=NA)

# r1 ------------------------------------------------------------------------------------------------------

#Exand cells with NA to datacollection area
r1_ext <- extend(r1, extent(grid_na), val=NA)

# function to fill na cells
fill.na <- function(x) {
  center = 0.5 + (width*width/2) 
  if( is.na(x)[center] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[center],0) )
  }
}  

width = 3

r <- r1_ext
it1 = 50
for (i in c(1:it1)){
r <- focal(r, w = matrix(1,width,width), fun = fill.na, 
            pad = TRUE, na.rm = F)
}

#r1_proj <- projectRaster(r1, grid_1km)
#r1_mask <- mask(r1_proj, grid_1km)
#r1_mosaic <- merge(r1_mask, r1_50mask_mask)
#r1_mosaic <- merge(r1_50mask, r1_50mask_mask)

r1_50mask <- projectRaster(r, grid_1km)
r1_50mask_mask <- mask(r1_50mask, grid_1km)


outfile <- paste(work1_path, subtheme1, "_", names(r1), "_interp", "_1km_", version, d_version,".tif", sep="")
outfile
writeRaster(r1_50mask_mask, outfile, COMPRESS=LZW, overwrite=T)


r_transform <- log(r1_50mask_mask)
r_transform <- r_transform/maxValue(r_transform)

outfile <- paste(work1_path, subtheme1, "_", names(r1), "_interp_log", "_1km_", version, d_version,".tif", sep="")
outfile
writeRaster(r_transform, outfile, COMPRESS=LZW, overwrite=T)




# r2 ------------------------------------------------------------------------------------------------------

#Exand cells with NA to datacollection area
r2_ext <- extend(r2, extent(grid_na), val=NA)

# function to fill na cells
fill.na <- function(x) {
  center = 0.5 + (width*width/2) 
  if( is.na(x)[center] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[center],0) )
  }
}  

width = 3

r <- r2_ext
it1 = 50
for (i in c(1:it1)){
  r <- focal(r, w = matrix(1,width,width), fun = fill.na, 
             pad = TRUE, na.rm = F)
}


r2_50mask <- projectRaster(r, grid_1km)
r2_50mask_mask <- mask(r2_50mask, grid_1km)


outfile <- paste(work1_path, subtheme1, "_", names(r2), "_interp", "_1km_", version, d_version,".tif", sep="")
outfile
writeRaster(r2_50mask_mask, outfile, COMPRESS=LZW, overwrite=T)

r_transform <- log(r2_50mask_mask)
r_transform <- r_transform/maxValue(r_transform)

outfile <- paste(work1_path, subtheme1, "_", names(r2), "_interp_log", "_1km_", version, d_version,".tif", sep="")
outfile
writeRaster(r_transform, outfile, COMPRESS=LZW, overwrite=T)






#grid_na <- resample(grid_na, r1)
#r1_focal <- mosaic(r1, grid_na, fun="max")
r1_focal <- r1
r2_focal <- mosaic(r2, grid_na, fun="max")

gf <- focalWeight(r1_focal, d=0.2, type='Gauss')


r1_f <- focal(r1_f, gf, fun="mean", NAonly = F, pad=T, na.rm=T)
#test <- focal(r1_focal, w=matrix(1/9,nrow=3,ncol=3)) 






# additional resource - create raster example ----------------------------------------------

# in your case (no need to compute rows/colums):
#x <- raster(xmn=-125, xmx=-65, ymn=25.5, ymx=50.5, res=0.75, crs="+proj=longlat +datum=WGS84")

#us_fire <- rasterize(temp[, c('lon', 'lat')], x, temp[, 'temp'], fun=mean)
#plot(us_fire)



# OUTPUTS -----------------------------------------------------------------------------------------
# OBS only write final files to the root directory using the exampels from here and on, any files generated in the 
# "PROCESS Section will lake proper metadata (i.e. ..._sourcesym.txt) should be written only to ../proc/

# Check outputs  --------------------------------------------------------------------------------
# If needed, check your outputs in /proc/ before writing to file. In this case a manual check in ArcGIS was done.


# write to product file: product 1  -------------------------------------------------------------------------------------- 
# if multiple products are produced in a loop with same sources they can be written to file in one chunk (i.e. this section), 
# just make sure to save one sourcesym file for each individual produc. 

scale <- "1km"
descriptive_name <- "nameofproduct"
product1_path <- product_path <- paste(dest1_path, theme1, "_", subtheme1, "_", descriptive_name, "_", scale, "_", version, d_version, ".tif", sep="" )
product_path

# for normalised products, name like this
#product1_path_norm <- product_path <- paste(dest1_path, theme1, "_", subtheme1, "_", descriptive_name, "_norm01_", scale, "_", version, d_version, ".tif", sep="" )
#product1_path_norm



# create sourcesym file to track all sources used in each output file (manualy list all "source_id" relevant objects, make sure to double check all your IDs are in there)
file1_data_raw_sources <- file_data_raw_sources <- tibble(id = c("uniqueid1", "uniqueid2", "uniqueid3"))

file1_sourcesym <- file_sourcesym <- file_data_raw_sources %>% 
#  add_row(data1_sourcesym)%>%  # add any "sourcesym" files applicable to this product, keep adding rows with data sourcesym files as needed
#  add_row(data2_sourcesym)%>% 
  unique() %>% 
  print()

file_sourcesym # double check that that all uniqueIDs for our data product sources are included 

# write to file 
write_tsv(file_sourcesym, paste(product_path, "_sourcesym.txt", sep=""))  # write unique sourcesym files for all products

product <- nameofobject # input and check your object to be written to file

writeRaster(product, product_path, COMPRESS=LZW, datatype = 'INT4S', overwrite=TRUE)


# write to product file: product 2  -------------------------------------------------------------------------------------- 

product2_path <- product_path <- paste(dest1_path, theme1, "_", subtheme1, "_additional_description", "_", version, d_version, ".tif", sep="" )
product_path

# create sourcesym file, not needed if its the same as product1
file_sourcesym # check so still sourcesym still, applies (e.g. identical sources as product 1), if not copy and adjust code from file 1 to create new file. 

# write to file 
write_tsv(file_sourcesym, paste(product_path, "_sourcesym.txt", sep=""))  # write unique sourcesym files for all products
product <- nameofobject # input and check your object to be written to file
writeRaster(product, product_path, COMPRESS=LZW, datatype = 'INT4S', overwrite=TRUE)

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





