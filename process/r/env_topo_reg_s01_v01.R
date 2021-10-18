
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
x <- c("tidyverse", "sf", "raster", "rgdal", "fasterize", "labelled", "gdalUtils")

install.packages(x, dependencies=TRUE)
                                            
library(tidyverse)
library(sf)
library(raster)
library(rgdal)
library(fasterize)
library(labelled)
library(gdalUtils)


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

theme1 <- theme <- "env" # input your main theme e.g. "eco" (choose from list in "theme_folder")

# select subtheme / component

folders %>%
  filter(theme_folder == theme) %>% 
  dplyr::select(subtheme_folder, subtheme_name) %>% 
  print(n = Inf)

subtheme1 <- subtheme <- "topo" # input your subtheme based on the most appropriate option from the list. If your subtheme is missing, please add new subtheme to folders.txt first (and make sure to sync with shared location for WIOSym if applicable)

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
source1_dir <- source_dir <-  "./data_raw/reg/env/topo/gebco/gk2105281356/"
source1_id <- source_id <- "gk2105281356" # cut and past your source id here
source1_metasym <- read_tsv(paste(source_dir, source_id, "_metasym.txt", sep="")) # reads the associated metasym file
source1_metasym # check that metasym file is populated, and fix if not
# check what data is in data_raw directory and read in all files to use from the directory:
dir(source_dir, recursive=T) 
#File1: bathymetry from gebco
dir(source_dir) # check what data is in data_raw directory
source_file <- "GEBCO_2020_20_Aug_2020_771a55d8e1ac/gebco_2020_n18.0_s-42.0_w6.0_e80.0.tif"  # your data file(s)
source1_path1 <- source_path <- paste(source_dir, source_file, sep="")
source_path
gebco_tid <- inraster1 <- raster(source_path)

#File2: bathymetry uncertainty from gebco
dir(source_dir) # check what data is in data_raw directory
source_file <- "GEBCO_2020_20_Aug_2020_771a55d8e1ac/gebco_2020_tid_n18.0_s-42.0_w6.0_e80.0.tif"  # your data file(s)
source1_path2 <- source_path <- paste(source_dir, source_file, sep="")
source_path
gebco <- inraster2 <- raster(source_path)

  
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


# data5: bounding box shepfile
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
# from here and on you are free to write your own process, annotate well and use the style guide ref above when applicable 

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


# fit to old script ------------------------------------------------------------------

path_raw <- source1_dir 
path_raw
path_work <- work1_path
path_work

inraster <- gebco
inraster

inraster_path <- source1_path1
inraster_path
# warp depth to 250m grid --------------------------------------------------------------------

# write empty grid for gdalutil work
writeRaster(grid_1km_na, paste(path_work, "grid_1km_na_mean_depth_gebco2020.tif", sep= ""),  overwrite=T, COMPRESS=LZW)
writeRaster(grid_250m_na, paste(path_work, "grid_250m_na_mean_depth_gebco2020.tif", sep= ""),  overwrite=T, COMPRESS=LZW)

# write to file
dstfile <- paste(path_work, "grid_250m_na_mean_depth_gebco2020.tif", sep="")

crs(grid_1km)
# 250m grid transformation from high resolution water mask

sr <- crs(inraster)
tr <- crs(grid_1km)
#  "+proj=cea +lat_ts=-12 +lon_0=12 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"    

bat_250m <- gdalUtils::gdalwarp(srcfile = inraster_path,
                                dstfile = dstfile,
                                s_srs = sr,
                                t_srs = tr,
                                tr = c(250, 250),
                                tap = TRUE,
                                output_Raster = FALSE,
                                overwrite = TRUE,
                                r = "average",
                                multi = TRUE,
)

bat_250m <- raster(dstfile)

writeRaster(bat_250m , paste(path_work, "grid_250m_na_mean_depth_gebco2020_lzw.tif", sep=""), overwrite=F, COMPRESS=LZW)


dstfile <- paste(path_work, "grid_1km_na_mean_depth_gebco2020.tif", sep="")

bat_1km <- gdalUtils::gdalwarp(srcfile = inraster_path,
                               dstfile = dstfile,
                               s_srs = sr,
                               t_srs = tr,
                               tr = c(1000, 1000),
                               tap = TRUE,
                               output_Raster = TRUE,
                               overwrite = TRUE,
                               r = "average",
                               multi = TRUE,
                               srcnodata = "NA"
)

summary(bat_1km)

# NOT RUN high resolution depth (Allen Atlas) to 250m grid --------------------------------------------------------------------

# this data and code is not yet complete... needs to be reworked, gdalutils needs to be modified somehow and raster::aggregate is to slow. Also, source file have problems

# write empty grid for gdalutil work
#writeRaster(grid_1km_na, paste(path_work, "grid_1km_na_mean_depth_gebco2020.tif", sep= ""),  overwrite=T, COMPRESS=LZW)
writeRaster(grid_250m_na, paste(path_work, "grid_250m_na_mean_depth_allen2020.tif", sep= ""),  overwrite=T, COMPRESS=LZW)


# write to file

dstfile <- paste(path_work, "grid_250m_na_mean_depth_allen2020.tif", sep="")

crs(grid_1km)
# 250m grid transformation from high resolution water mask
inraster_path <- paste("./data_raw/env/depth/reg/allen/20200902/depth_mapping_west_indian_ocean_west_indian_ocean_islands_high_tide_normalized_sr_jan2018_jan2020_v2_shift_mosaic_20200624T151707Z_depth_composite_depth.tif", sep="")
inraster <- raster(inraster_path)

inraster[inraster == 0] <- NA

sr <- crs(inraster)
tr <- crs(grid_250m)
#  "+proj=cea +lat_ts=-12 +lon_0=12 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"    

bat_250m <- gdalUtils::gdalwarp(srcfile = inraster_path,
                                dstfile = dstfile,
                                s_srs = sr,
                                t_srs = tr,
                                tr = c(250, 250),
                                tap = TRUE,
                                output_Raster = FALSE,
                                overwrite = TRUE,
                                r = "average",
                                multi = TRUE,
)

bat_250m <- raster(dstfile)

writeRaster(bat_250m , paste(path_work, "grid_250m_na_mean_depth_allen2020_lzw.tif", sep=""), overwrite=T, COMPRESS=LZW)

inraster
inraster_agg10 <- raster::aggregate(inraster, fact=10, fun=mean, na.rm=T)
writeRaster(inraster_agg10 , paste(path_work, "depth_allen2020_sourcefile_agg10.tif", sep=""), overwrite=T, COMPRESS=LZW)

# mask depth to analysis grid -------------------------------------------------------------------------------------------------------
# read indata
bat_250m <- raster(paste(path_work, "grid_250m_na_mean_depth_gebco2020_lzw.tif", sep=""))
grid_250m


# crop/mask file to wiosym grid (three steps + 1 compression step, first crop extent and stack to include the mask band, the write stack to file, then use gdal utils to mask (for performance reasons), finally rewrite as compressed raster

# NOT USE masking raster to grid alteriantive 1 ------------------------------------------------------------------------------------
bat_250m_mask <- crop(bat_250m, extent(grid_250m))
bat_250m_stack <- stack(bat_250m_mask, grid_250m)

writeRaster(bat_250m_stack , paste(path_work, "grid_250m_na_mean_depth_gebco2020_maskband.tif", sep=""), overwrite=T, COMPRESS=LZW)

inraster_path <- paste(path_work, "grid_250m_na_mean_depth_gebco2020_maskband.tif", sep="")
dstfile <- paste(path_work, "env_depth_mean_grid_250m_v00.tif", sep="")

#write raster for gdal utils output
writeRaster(grid_250m_na, paste(path_work, "env_depth_mean_grid_250m_v00.tif", sep= ""),  overwrite=T, COMPRESS=LZW)

raster(dstfile)

sr <- crs(grid_250m)
tr <- crs(grid_250m)


bat_250m_mask2 <- gdal_translate(src_dataset = inraster_path,
                                 dst_dataset = dstfile,
                                 s_srs = sr,
                                 t_srs = tr,
                                 tr = c(250, 250),
                                 b = 1,
                                 mask = 2,
                                 output_Raster = FALSE,
                                 overwrite = TRUE,
                                 r = "nearest",
                                 co = "COMPRESS=LZW",  # info about compressions https://digital-geography.com/geotiff-compression-comparison/
                                 stats = TRUE
)


# masking raster to grid alteriantive 2 ------------------------------------------------------------------------------------
bat_250m_mask <- crop(bat_250m, extent(grid_250m))
bat_250m_mask <- mask(bat_250m_mask, grid_250m)

# reclassify values larger then 0 in bathymetry
bat_250m_mask[bat_250m_mask>1] <- 1


writeRaster(bat_250m_mask, paste(path_work, "env_depth_mean_grid_250m_", version, ".tif", sep=""), overwrite=T, COMPRESS=LZW, datatype = 'INT4S')

# aggregate to 1 km cells using the 250m grid

bat_1km <- raster::aggregate(bat_250m_mask, fact=4, fun=max, na.rm=T)
bat_1km
plot(bat_1km)
# write to file

writeRaster(bat_1km, paste(path_work, "env_depth_mean_grid_1km_", version, ".tif", sep=""), overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')



# Divide depth into individual rasters with depth zones ------------------------------------------------------------------------------------

# read data from disk again and change to shorter name

depth_250m <- raster(paste(path_out, "env_depth_mean_grid_250m_", version, ".tif", sep=""))
depth_1km <- raster(paste(path_out, "env_depth_mean_grid_1km_", version, ".tif", sep=""))

# explore data
hist(depth_1km,
     main = "Distribution of raster cell values in the WIO Gebco depth data",
     #breaks = c(0, -40, -200, -1000, -3000, -7000),
     xlab = "Depth (m)", ylab = "Area",
     col = "springgreen")

# create classification matrix for depth data using depth zones
reclass_df <- c(-Inf, -3000, 5,
                -3000, -1000, 4,
                -1000, -200, 3,
                -200, -40, 2,
                -40, Inf, 1)

reclass_df
#Reorder into matrix
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)

reclass_m

depth_250m_classified <- reclassify(depth_250m,
                                    reclass_m, include.lowest=TRUE)


depth_1km_classified <- reclassify(depth_1km,
                                   reclass_m, include.lowest=TRUE)

writeRaster(depth_250m_classified, paste(path_work, "env_depth_zone_grid_250m_", version, ".tif", sep=""), overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')
writeRaster(depth_1km_classified, paste(path_work, "env_depth_zone_grid_1km_", version, ".tif", sep=""), overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')


# export individual zones -----------------------------------------------------------------------
# create classification matrix for depth data using depth zones
reclass_df_zone <- c(1, 2, 1,
                     2, 5, NA)

reclass_df_zone
#Reorder into matrix
reclass_m_zone <- matrix(reclass_df_zone,
                         ncol = 3,
                         byrow = TRUE)

depth_250m_classified_zone <- reclassify(depth_250m_classified,
                                         reclass_m_zone, include.lowest=TRUE)


writeRaster(depth_250m_classified_zone, paste(path_work, "env_depth_0-40m_zone_grid_250m_", version, ".tif", sep=""), overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')


# zone 40-200m
reclass_df_zone <- c(1, 3, 1,
                     3, 5, NA)

reclass_df_zone
#Reorder into matrix
reclass_m_zone <- matrix(reclass_df_zone,
                         ncol = 3,
                         byrow = TRUE)

depth_250m_classified_zone <- reclassify(depth_250m_classified,
                                         reclass_m_zone, include.lowest=TRUE)

writeRaster(depth_250m_classified_zone, paste(path_work, "env_depth_0-200m_zone_grid_250m_", version, ".tif", sep=""), overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')



# zone 200-infm
reclass_df_zone <- c(1, 3, NA,
                     3, 5, 1)

reclass_df_zone
#Reorder into matrix
reclass_m_zone <- matrix(reclass_df_zone,
                         ncol = 3,
                         byrow = TRUE)

depth_250m_classified_zone <- reclassify(depth_250m_classified,
                                         reclass_m_zone, include.lowest=TRUE)

writeRaster(depth_250m_classified_zone, paste(path_work, "env_depth_200-Inf_zone_grid_250m_", version, ".tif", sep=""), overwrite=TRUE, COMPRESS=LZW, datatype = 'INT4S')

# calculate raster statistics and ratify -----------------------------------------------------
r_mosaic <- depth_250m_classified


df <- as.data.frame(r_mosaic)
df <- as_tibble(df)   %>%  
  mutate_all(as.integer) %>% 
  drop_na()
resolution = 250
stat <- df %>% 
  group_by(env_depth_mean_grid_250m_v00) %>% 
  summarize(n = n()) %>%
  mutate(km2 = round(n*resolution*resolution/1000000, 2)) %>% 
  mutate(percent = round(100 * n/sum(n), 2)) %>% 
  print()

r_mosaic <- ratify(r_mosaic)
r_mosaic_rat <- levels(r_mosaic)[[1]]

Value <- c(1, 2, 3, 4, 5)
class_name <- c("0-40m", "40-200m", "200-1000m", "1000m-3000m", ">3000m")
count <- as.factor(c(stat[["n"]][1],	stat[["n"]][2],	stat[["n"]][3],	stat[["n"]][4],	stat[["n"]][5]))
class_km2 <- as.factor(c(stat[["km2"]][1],	stat[["km2"]][2],	stat[["km2"]][3],	stat[["km2"]][4],	stat[["km2"]][5]))

stat <- stat %>%
  rename(Value = env_depth_mean_grid_250m_v00)

rat_table <- data.frame(Value, class_name) %>%
  as_tibble() %>% 
  mutate(Value = as.integer(Value)) %>% 
  left_join(stat, by="Value")

rat_table <- rat_table %>% 
  filter(n>=0) %>% 
  print()



r_mosaic_rat$Value <- rat_table$Value
r_mosaic_rat$class_name <- rat_table$class_name
r_mosaic_rat$count <- rat_table$n
r_mosaic_rat$class_km2 <- rat_table$km2
r_mosaic_rat$class_perc <- rat_table$percent
levels(r_mosaic) <- r_mosaic_rat
r_mosaic_rat

write.dbf(r_mosaic_rat, file = paste(path_work, "env_depth_zone_grid_250m_", version, ".m.tif.vat.dbf", sep=''), factor2char = TRUE, max_nchar = 254)



# OUTPUTS -----------------------------------------------------------------------------------------
# OBS only write final files to the root directory using the exampels from here and on, any files generated in the 
# "PROCESS Section will lake proper metadata (i.e. ..._sourcesym.txt) should be written only to ../proc/

# Check outputs  --------------------------------------------------------------------------------
# If needed, check your outputs in /proc/ before writing to file. In this case a manual check in ArcGIS was done.


# write to product file: product 1  -------------------------------------------------------------------------------------- 
# if multiple products are produced in a loop with same sources they can be written to file in one chunk (i.e. this section), 
# just make sure to save one sourcesym file for each individual produc. 

product1_path <- product_path <- paste(dest1_path, theme1, "_", subtheme1, "_additional_description", "_", version, d_version, ".format" )
product_path

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

product2_path <- product_path <- paste(dest1_path, theme1, "_", subtheme1, "_additional_description", "_", version, d_version, ".format" )
product_path

# create sourcesym file, not needed if its the same as product1
file_sourcesym # check so still sourcesym still, applies (e.g. identical sources as product 1), if not copy and adjust code from file 1 to create new file. 

# write to file 
write_tsv(file_sourcesym, paste(product_path, "_sourcesym.txt", sep=""))  # write unique sourcesym files for all products

product <- nameofobject # input and check your object to be written to file

writeRaster(product, product_path, COMPRESS=LZW, datatype = 'INT4S', overwrite=TRUE)

# write to product file: product 3...  -------------------------------------------------------------------------------------- 
# keep adding products according to above example

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





