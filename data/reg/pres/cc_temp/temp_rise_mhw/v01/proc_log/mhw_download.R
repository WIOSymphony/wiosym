library(tidyr)
library(dplyr)
library(httr)
library(rvest)

version = "v01"
d_version = ".0"
location = "reg"   #check
theme1 <- "pres"   #check
subtheme1 <- "cc_temp"   #check
name1 <- "temp_rise_mhw"   #check
name2 <- "" #check
(dest1_path <- paste("./data", location, theme1, subtheme1, name1, version, sep="/"))
(work1_path <- paste(dest1_path, "proc", sep="/"))
(proc1_path <- paste(dest1_path, "proc_log", sep="/"))
script_path <- "./process/r"
(script_path_windows <- paste(getwd(), "/process/r", sep=""))
(script_name <- paste(theme1, "_", subtheme1,"_", location, "_", name1,"_s01_",version, ".Rmd", sep=""))
(dest1_path_file1 <-paste0(dest1_path, "/", subtheme1, "_", name1, "_", version, d_version, ".tif"))

# Repeat for each year - data available from 1985
for(i in 1985:2023){

  yr <- i
  dir.create(paste0(work1_path, "/mhw_noaa/", yr))
  url <- paste0("https://www.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/marine_heatwave/v1.0.1/category/nc/", yr, "/")
  page <- GET(url)
  html_content <- content(page, "text")
  
  links <- read_html(html_content) %>%
    html_nodes("a") %>%
    html_attr("href") %>%
    na.omit()  # Remove NAs
  
  nc_links <- links[grepl("\\.nc$", links)]
  nc_links <- paste0(url, nc_links)
  
  for(link in nc_links){
    fn <- basename(link)
    fp <- paste0(work1_path, "/mhw_noaa/", yr, "/", fn)
    download.file(url = link, destfile = fp, mode = 'wb')
}
}