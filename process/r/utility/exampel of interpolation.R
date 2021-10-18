#### Metrics from high resolution depth and backscatter ####
install.packages("fields")



library(tidyverse)   # installs sweeps of packages
library(sf)          # working with sf objects in R
library(raster)      # working with rasters in R
library(rgdal)       # links to the GDAL library
library(stars)	   # new package for large rasters and datacubes
library(fields)


#steps to do-----------------------------------------
#1. rastermask 5m, minimum geometry
#2. Fill/interpolate within mask extent for depth and bs high resolution
#3. resample depth and bs to rescaling (1m, 2,5, 5m), same pixels extents as 5m mask
#4. develop predictors for <=1m rasters, aggregate to 5m (mean, and some other variables)
#5. develop predictors for 5m and 2.5m


#setwd("C:/Users/guskag/Desktop/R-course/hob/")	#Set wd, remove when in R project

#Read in source data-----------------------------------------------
depth_raw <- raster("./data_raw/depth_50cm.tif")

bs_raw <- raster("./data_raw/bs_50cm.tif")

mask_raw <- read_sf("./data_raw/mask.shp") %>%
       st_set_crs("+proj=utm +zone=33 +ellps=GRS80") 

#aggregate if to change resolution
depth_raw <- aggregate(depth_raw, 2)
bs_raw <- aggregate(bs_raw, 2)

mask_res <- 5 #set to resolution of modelling rasters

raw_res <- sum(res(depth_raw))/2  #resolution of bathy/bs raw
raw_res

# Julian Code to create rastermask from shp ------------------------

extent(mask_raw) #modify extent belowe to be whole numbers (multiples of resolution)

target <- raster(xmn = 710345 - 10 * mask_res, 
                 xmx = 712210 + 10 * mask_res,
                 ymn = 6271480 - 10 * mask_res,
                 ymx = 6273820 + 10 * mask_res,
                 crs = "+proj=utm +zone=33 +ellps=GRS80",
                 res = mask_res)
extent(target)
mask <- rasterize(mask_raw, target) #creates raster mask out of polygon with modelling target resolution, buffer with 10 x cellsize
mask_buff <- buffer(mask, width = 10, doEdge=FALSE)
plot(mask)
plot(mask_buff)

#Create mask for high resolution depth/bs data

mask_raw_nodes <- mask_res/raw_res
mask_raw_nodes
mask_raw <- disaggregate(mask, fact=mask_raw_nodes)


#### Create buffer around NA values in depth/ bs data, and extract those

#adjust extent and mask raw raster

holes <- ra <- depth_raw %>%
			resample(mask_raw) %>%
			mask(mask_raw)

#create inverse masked where holes exist

holes[!is.na(holes[])] <- 0 
holes[is.na(holes[])] <- 1 
holes[holes[] == 0] <- NA 
plot(holes)


#grow raster 
buffer <- buffer(holes, width=10, doEdge=TRUE)
plot(buffer)

#extract values around buffer zones
ra_buffer <- mask(ra, buffer, inverse=FALSE)
plot(ra_buffer)

##interpoleringstest


a <- which(ra_buffer[]<=Inf) #extract cell numbers for all raster cells with values

xy <- ra_buffer %>% 
	xyFromCell(a)%>%
	as_tibble()

v <- ra_buffer[a] #extract values from cells 

xyv <- mutate(xy,v) #combine xy and value




xyv <- Tps(xy, v)

glimpse(which(ra_buffer[]!=1))

ra_interp <- interpolate(ra_buffer, xyv, ext = ra_buffer) %>%
  	mask(holes)
plot(ra_interp)

# mosaic to combined interpolated raster
ra_mosaik <- merge(ra, ra_interp)

plot(ra_mosaik)




# TEST ON LOW RES RASTER------------------------------------------------------------------

depth_100m <- aggregate(depth_raw, 100)

# not used # bs50cm.fill <- rfillspgaps(bs_50cm, maskPol=NULL, nmax=50) 
# not used # st_buffer #buffer line feature


#Create test dataset (not needed in normal scenario)----------------------
# Create an empty raster with the same geometry
my.mask <- depth_100m
my.mask[] <- 1 # Put some non-NA value

set.seed(100)
cells <- sample(1:ncell(my.mask), 20)
my.mask[cells] <- NA

rst_masked <- mask(depth_100m, my.mask)
#--------------------------------------------------------------------------









##code from francis


##interpoleringstest
ra <- depth_100m
xy <- tibble(data.frame(xyFromCell(ra, 1:ncell(ra))))
v <- tibble(getValues(ra))
tps <- Tps(xy, v)
depth_100m_interp <- interpolate(ra, tps)

ras<-as(ra, 'SpatialPixels')
rasd<-as(ra, 'SpatialGridDataFrame')





