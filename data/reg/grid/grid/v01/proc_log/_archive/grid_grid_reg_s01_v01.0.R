
#### ABOUT ####
# ========================================================================================================= #
# Purpose: Script to create WIOSym regional standard grid
#
# Brief description: v01.0 uses srtm data that has already been modified in ARcGIS to represent the marine areas as "1" 
#   and land areas as NA. The purpose of this script is to exactly define the extent, coordinate system and resolution 
#   of the standard grid which will be used for all products in version 1, as well as provide multiple resolutions
#
# Suggestions on future imporvements: Add more detailed data describing estuaries, suc as described here: https://rpubs.com/MRufino/World_Estuaries
#   Also, need to include country specific coastline definitions, such as SA has already proveded but its not yet implemented here. Most important 
#   is to not exclude areas in analysis grid, so a maximum fotprint of approved files would be one way of combining sources
#
# Created: GK, SGU, 20210430
# Latest update: GK, SGU, 20210506

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
# Read packages ---------------------------------------------------------------------------------------------
# list all packages that this script is dependent on, remove any redundant ones
#x <- c("tidyverse", "sf", "raster", "rasterVis", "rgdal", "stars", "parallel", "gdalUtils")

#install.packages(x, dependencies=TRUE)
install.packages("snow")
                                            
library(tidyverse)
library(sf)          # working with sf objects in R
library(raster)      # working with rasters in R
#library(rasterVis)   # raster visualization
library(rgdal)       # links to the GDAL library
library(parallel)
library(gdalUtils)
library(fasterize)
#library(stars) 
library(snow)

# Set data_raw sources ----------------------------------------------------------------------------------------

# set all sources and check id/metadata
  # OBS! only use relative filepaths e.g.  ./data/...  
  # The workdirectory is automatically set to.../wiosym/ when you use start Rstudio by clicking on the *.rproj file under .../wiosym/ 
  # this ensure that all paths are relative so that your data and code can be moved to a new location and still function as is. 
                                            
# declare paths to all source directories (data_raw) and check that all sources have "metasym" files

# path 1 (repeat same code chunk for any additional source)
source1_dir <-  "./data_raw/reg/env/topo/nasa/gk2104301347/"
source1_id <- "gk2104301347"
source1_metasym <- read_tsv(paste(source1_dir, source1_id, "_metasym.txt", sep=""))
source1_metasym # check that metasym file is populated, and fix if not - before finalising the processing

# source1_object1 <- read()
# source1_object2 <- read()


# Set data sources -------------------------------------------------------------------------------------------------
# declare paths to all sources from "data" directory. make sure that all sources are accompanied by the a datasym.txt file, if not fix...

# data1: watermask
data1_dir <- "./data/reg/grid/mask/v00/"
dir(data1_dir)
data1_file <- "srtm_0m_marine.tif" # this data was created through through Global Mapper using srtm data and some manual work, future update of this script could include direct access of srtm data/mask from R
data1_sourcesym <- read_tsv(paste(data1_dir, data1_file, "_sourcesym.txt", sep=""))
data1_sourcesym # check so source IDs exists and make sense
data1_path <- paste(data1_dir, data1_file, sep="")
data1_path # note: this data is not read in as an object since the process in this script happens in GDAL outside R
data1_object <- raster(data1_path)

# data2: bounding box for the grid (created by SGU in GIS so no data_raw sources involved)
data2_dir <- "./data/reg/grid/scope/v01/" 
dir(data2_dir)
data2_file <- "wiosym_grid_bounding_box_v01.shp" 
data2_sourcesym <- read_tsv(paste(data2_dir, data2_file, "_sourcesym.txt", sep=""))
data2_sourcesym # check so source IDs exists and make sense
data2_path <- paste(data2_dir, data2_file, sep="")
data2_path # note: this data is not read in as an object since the process in this script happens in GDAL outside R

data2_bounding_box <- st_read(data2_path)   # "gridbox" in v00. Bounding box for processing area created in ArcGIS by SGU, and approved through twg workshop 1 Nov 2020 


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
theme1 <- "grid" # e.g. theme_folder

# select subtheme / component

folders %>%
  filter(theme_folder == theme1) %>% 
  print(n = Inf)

subtheme1 <- "grid" # e.g. subtheme_folder. If an appropriete subtheme is missing, please add new subtheme to folders.txt first (and make sure to sync with shared location for WIOSym if applicable)
  
dest1_path <- paste("./data/", location, "/", theme1, "/", subtheme1, "/", version, "/", sep="")  # location for any final products that can be used elsehwere along with datasym.txt file which logs all source IDs. 
dest1_path

if (!dir.exists(dest1_path)){
  dir.create(dest1_path, recursive = TRUE)
}

# work directory for any temporary files
work1_path <- paste(path_dest1, "proc/", sep="") # in this folder you can organise any files you use in the process using your prefferred structure. No "final products" here, they go under root (path_dest1)
work1_path

if (!dir.exists(work1_path)){
  dir.create(work1_path, recursive = TRUE)
}

# Process log directory, for saving process logs and the exact script (or mxd etc) used to create the output files
proc1_path <- paste(path_dest1, "proc_log/", sep="") # in this folder you can organise any files you use in the process using your prefferred structure. No "final products" here, they go under root (path_dest1)
proc1_path

if (!dir.exists(proc1_path)){
  dir.create(proc1_path, recursive = TRUE)
}

# ... if multiple themes keep adding destinations here, and not futher down in script


# Save your R script ---------------------------------------------------------------------------------------------------------------------------------------
#location for all active R scripts (latest version) according to theme, subtheme
script_path <- "./process/r"  # main path
script_path
if (!dir.exists(script_path)){
  dir.create(script_path, recursive = TRUE)
}

# name of our script
script_name <- paste(theme1, "_", subtheme1,"_", location, "_s01_",version, sep="")  # s  (e.g. "_s01_") stands for sequential order incase you have divided the work in several scripts

# run lines below, copy path and name and use File/Save As in R studio to save your .R workfile in the correct workspace - 
# Note: in the end of the script (e.g. "write to file section") you are also asked to save a copy of your final code in the destination proc_log directory

script_path_windows <- paste(getwd(), "/process/r/", sep="")  # filepath

# script_path_windows <- gsub("/", "\\\\", script_path_windows)
# script_path_windows

paste(script_name, ".R", sep="") # file name



# PROCESS ----------------------------------------------------------------------------------------------------------------
# from here and on you are free to write your own process, annotate well and use the style guide ref above when applicable 

# Set coordinates --------------------------------------------------------------
# cylindircal equal area projection (same as lambert cylindrical equal area), with center shifted to mid lattitude for wiosym, and start longitude at min x (i.e. 0m in the far left of the grid)
sr <- "+proj=cea +lat_ts=-12 +lon_0=12 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"    


# Define bounding box in projected coordinates ----------------------------------------------------------------------------
res_30arcsec = 30/3600 # 30 arc seconds is roughly the resolution we are aiming for on regional basis

#grid bounding box in lat/long
grid_30sec <- raster(xmn=8, xmx=78, ymn=-40, ymx=16, 
                     crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", resolution = res_30arcsec, vals=NA)

#reproject lat long to equal area projection defined above
grid_30sec_proj <- projectRaster(grid_30sec, crs = sr)

#reproject grid bounding box
grid_box_proj <- st_transform(data2_bounding_box, crs=sr)


#plot grid extent in the projected coordinate system
extent(grid_30sec_proj)

# adjust extents slightly manually to stay on whole numbers etc, visual check in arcgis was performed as well.
extent_grid <- extent(c(-441000, 7192000, -4175000, 1790000))

# 250m grid -----------------------------------------------------------------------------------------------------------------

# reclassify NA values in watermask so that they can be summed by gdalutils
data1_object

#m <- c(1, NA, 0)
#m <- matrix(m, ncol=3, byrow=TRUE)
#data1_object_nona <- reclassify(data1_object, m)
#data1_object_nona


# optional reclassify
replaceNA <- function(x, na.rm, ...){ 
  if(is.na(x[1]))
    return(0)
  else
    return(x)
} 



beginCluster()
data1_object_nona <- clusterR(data1_object, calc, args=list(fun=replaceNA))#, export='a')
endCluster()


data1_object_nona <- calc(data1_object, fun = replaceNA)
data1_path
data1_work_path <- paste("./data/reg/grid/mask/v00/proc/", sep="") 

if (!dir.exists(data1_work_path)){
  dir.create(data1_work_path, recursive = TRUE)
}
data1_file
data1_work_path_nona <- paste(data1_work_path, "srtm_0m_marine_nona.tif")
# write to file
writeRaster(data1_object_nona, data1_work_path_nona, overwrite=TRUE)

# create grid frame for 250m  resolution
grid_proj_250m <- raster(extent_grid,
                    crs=sr, resolution = c(250,250), vals=1)

# write to file 
work1_path_250m_1 <- paste(work1_path, "grid_cea_250m.tif", sep="")
work1_path_250m_1
writeRaster(grid_proj_250m, work1_path_250m_1, overwrite=TRUE)


# 250m grid transformation from high resolution water mask
sr <- "+proj=cea +lat_ts=-12 +lon_0=12 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"    

bat_250m <- gdalUtils::gdalwarp(srcfile = data1_work_path_nona,
                           dstfile = work1_path_250m_1,
                           s_srs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                           t_srs = sr,
                           tr = c(250, 250),
                           tap = TRUE,
                           output_Raster = TRUE,
                           overwrite = TRUE,
                           r = "average",
                           multi = TRUE,
                           srcnodata = "NA"
)

bat_250m

# due to missing extent in the landmask on the west side we extend bat_250m to match grid area (its only missing ocean pixels on west side within area of intrest, so no need to update landmask source file)

extent_bat250m_extend <- extent(c(extent(extent_grid)[1], extent(bat_250m)[2], extent(bat_250m)[3], extent(bat_250m)[4]))
extent_bat250m_extend

bat_250m_extend <- raster(extent_bat250m_extend,
                            crs=sr, resolution = c(250,250), vals=255)

# new extended mask of the waterbody for the analysis grid
bat_250m_ext <- merge(bat_250m, bat_250m_extend, overlap=F)

# check so new updated file look ok
plot(bat_250m_ext)
bat_250m
bat_250m_ext

# overwrite old file with the extended

bat_250m <- bat_250m_ext

# turn grid box to raster with same extent
grid_box_250m_r <- fasterize(grid_box_proj[1], grid_proj_250m, field=NULL, fun="sum")

#make sure extents are identical
grid_box_250m_r2 <- crop(grid_box_250m_r, extent(bat_250m))

# crop tranformed file (two steps, first extent, then a mask)
bat_250m_crop <- crop(bat_250m, grid_box_250m_r2)
bat_250m_mask <- mask(bat_250m_crop, grid_box_250m_r2)

bat_250m_mask
plot(bat_250m_mask)

# reclassify values (for unknown reasons it converts to 255 during the transformation...)
m <- c(0, 255, 1)
#m <- matrix(m, ncol=3, byrow=TRUE)
bat_250m_mask <- reclassify(bat_250m_mask, m)

# write 250m grid to (work) file

#outfile1 <- paste(work1_path, subtheme1, "_250m_", version, d_version, ".tif", sep="")
outfile1 <- paste(work1_path, subtheme1, "_250m_rev1",".tif", sep="")
outfile1
writeRaster(bat_250m_mask, outfile1, overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S') 

# create na raster with same cells / extents
grid_250m <- raster(outfile1)
values(grid_250m) <- NA
outfile2 <- paste(work1_path, subtheme1, "_250m_na", ".tif", sep="")
writeRaster(grid_250m, outfile2, overwrite=F, COMPRESS=LZW, datatype = 'INT4S')



# 1km grid -----------------------------------------------------------------------------

# aggregate to 1 km cells using the 250m grid

grid_1km <- raster::aggregate(bat_250m_mask, fact=4, fun=max, na.rm=T)
grid_1km
plot(grid_1km)
# write to file
outfile3 <- paste(work1_path, subtheme1, "_1km", ".tif", sep="")
outfile3
writeRaster(grid_1km, outfile3, overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')


# create na raster with same cells / extents
grid_1km <- raster(outfile3)
grid_1km
values(grid_1km) <- NA
outfile4 <- paste(work1_path, subtheme1, "_1km_na", ".tif", sep="")
outfile4
writeRaster(grid_1km, outfile4, overwrite=F, COMPRESS=LZW, datatype = 'INT4S', owerwrite=T)

# Alternative way - diabled - create grid frame for 1km  "regional scale"
#grid_proj_1km <- raster(extent_grid,
#                        crs=sr, resolution = c(1000,1000), vals=1)
#dstfile2 <- paste(path_data, "proc/grid_cea_1km.tif", sep="")

#writeRaster(grid_proj_1km, dstfile2, overwrite=TRUE)


# 5km grid -----------------------------------------------------------------------------

# aggregate to 5 km cells using the 1km grid
grid_1km <- raster(outfile3)
grid_5km <- raster::aggregate(grid_1km, fact=5, fun=max, na.rm=T)
grid_5km
plot(grid_5km)
# write to file
outfile5 <- paste(work1_path, subtheme1, "_5km", ".tif", sep="")
outfile5
writeRaster(grid_5km, outfile5, overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')

# o
ptional : add eez to 250m grid  - disabled for now -----------------------------------------------

# reproject eez
#eez_proj <- st_transform(eez, crs=sr) 

# turn grid box to raster with same extent
#eez_proj_r <- fasterize(eez_proj, bat_250m_mask, field="MRGID", fun="first")

# make sure extents are identical
#grid_box_250m_r2 <- crop(grid_box_250m_r, extent(bat_250m))

# crop tranformed file (two steps, first extent, then a mask)
#bat_250m_crop <- crop(bat_250m, grid_box_250m_r2)
#bat_250m_mask <- mask(bat_250m_crop, grid_box_250m_r2)

#bat_250m_mask
#plot(bat_250m_mask)



# OUTPUTS -----------------------------------------------------------------------------------------
# OBS only write final files to the root directory using the exampels from here and on, any files generated in the 
# "PROCESS Section will lake proper metadata (i.e. ..._sourcesym.txt) should be written only to ../proc/

# Check outputs  --------------------------------------------------------------------------------
# If neeed, check your outputs in /proc/ before writing to file. In this case a manual check in ArcGIS was done.

# Write to file ---------------------------------------------------------------------------------

# naming template: dest1_path_file1 <-paste(dest1_path, theme1, "_", subtheme1, "_additional_description", "_", version, d_version, ".format" )
dest1_path_file1 <-paste(dest1_path, subtheme1, "_250m", "_", version, d_version, ".tif", sep="")

# create sourcesym file to track all sources used in each output file (manualy list all "source_id" relevant objects, make sure to double check all your IDs are in there)
# test: # file1_data_raw_sources <- tibble(id = c("gk2104301346", "gk2104301348", "gk2104301349"))
file1_data_raw_sources <- tibble(id = c("gk2104301347"))

file1_sourcesym <- file1_data_raw_sources %>% 
#  add_row(data1_sourcesym)%>% 
#  add_row(data2_sourcesym)%>% # keep adding rows with data sourcesym files as needed
  unique() %>% 
  print()

file1_sourcesym

write_tsv(file1_sourcesym, paste(dest1_path_file1, "_sourcesym.txt", sep=""))

writeRaster(grid_250m, dest1_path_file1, overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')

# Save copy of R scripts ----------------------------------------------------------------------------------

# 1. save your current R script (File/Save)

# 2. Save a copy of the R script here (if you update the files, make sure to redo the save):

#   run code below and go to File/Save As and paste link and name:

dest1_path_script <- paste(getwd(), proc1_path)
dest1_path_script <- gsub(" .", "", dest1_path_script) 

#dest1_path_script <- gsub("/", "\\\\", dest1_path_script) # code needs work...

dest1_path_script # path

script_name_copy <- paste(script_name, d_version, ".R", sep="")
script_name_copy # name


# OBS! If you make more edits after this "final" save make sure to go back to your original R script under /process/r again, then resave this one in the end


# FINAL CHECK  -------------------------------------------------------------------------------------------
# Do your own final check on the products and doublecheck sourcesym file exists and are complete for all files in root
# checked by: initial, org, date

# External check -------------------------------------------------------------------------------------
# To ensure repeatability and quality, make sure atleast one more colleage can run the script and make any final edits.
# When script is proven repeatable/understandable and products/metadata look ok 
# -> Sign at top of script in about section (# Approved by:), and make any comments if needed





