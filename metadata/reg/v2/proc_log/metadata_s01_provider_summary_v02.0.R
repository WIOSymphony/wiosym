#Script by gk & me Swedish WIO Symphony team
#Purpose: to make summary of all data sources used in tool (using sourcesym files updated from gitHUB) and write citations etc. to google sheet
#R v 4.1.3
#updated 20230209 by gk


install.packages("googlesheets4")

library(dplyr)
library(tidyverse)
library(rio)
library(googlesheets4)

# reading directly from google sheets, ask for permission to enter from Wiosym admin

layer_list <-read_sheet("https://docs.google.com/spreadsheets/d/1JIV92XisuBctBoVT53yqj0orPNRVyRQVEDtMq7R0fKw/edit#gid=0", sheet = "tool_metadata")

# extract all sourcesym names used for final products

sourcesym_list <- layer_list %>% dplyr::select(sourcesym) %>% dplyr::filter(!is.na(sourcesym)) %>% unique()

# list all metasym files and sourcesym files in wiosym cataloge

files_meta <- list.files(path = "./data_raw/", pattern = "_metasym.txt$", recursive=T)

files_source <- list.files(path = "./data/", pattern = "_sourcesym.txt$", recursive=T)

source_list <- c()


# extract filepath for all sourcesym files used in final products

for(i in sourcesym_list$sourcesym){
  file <- import(paste0("./data/", files_source[str_detect(files_source, pattern = i, negate = FALSE)]))
  #  file <- import(paste0("./metadata/wio_provider/wio_summary/Sourcesym/", i))
  source_list <- c(source_list, file$id)
  print(i)
}

#paste0("./data/", files_source[str_detect(files_source, pattern = i, negate = FALSE)])


# to check all sourcesym files in data, not needed if only to extract the used metasym files
#for(i in files_source){
#  file <- import(paste0("./data/", files_source[str_detect(files_source, pattern = i, negate = FALSE)]))
#  file <- import(paste0("./metadata/wio_provider/wio_summary/Sourcesym/", i))
#  source_list <- c(source_list, file$id)
#}


# combine information into a table by combining sourcesym IDs and metasym files, extract more information about providers in excel file (updated through wiosym github)

table(source_list)

meta_list <- tibble(id = c(), downloadDate = c(), provider = c(), provider_long = c(), citation = c(), citation_edit = c(), status = c(),theme = c(), subtheme = c(), tags = c(), freq = c())
fel <- c()
providers <- readxl::read_excel("./shiny_data_upload/modify_txt_files_v01.2.1.xlsx", sheet = "providers_master")

source_list <- unique(source_list) %>% discard(is.na)

for(j in unique(source_list)){  # run unique(source_list) to find missing file...
  tryCatch({
  file <- read_tsv(paste0("./data_raw/", files_meta[str_detect(files_meta, pattern = paste0(j, "_metasym.txt$"), negate = FALSE)]))
  dl_date  <- file[[3, 2]]
  provider  <- ifelse(file[[5, 2]] == "", file[[6, 2]], file[[5, 2]])
  provider_long <- providers$provider[which(providers$provider_val == provider)]
  if(is_empty(provider_long)){
    provider_long <- ""
  }
  citation <- file[[8, 2]]
  theme  <- file[[12, 2]]
  subtheme  <- file[[14, 2]]
  tags  <- file[[16, 2]]
  freq1 <- sum(j == source_list)
  tot <- tibble(id = j, downloadDate = dl_date, provider = provider, provider_long = provider_long, citation = citation, citation_edit = "", status = "", theme = theme, subtheme = subtheme, tags = tags, freq = freq1)
  meta_list <- rbind(meta_list, tot)
   })
}

meta_list <- arrange(meta_list, provider, downloadDate)
View(meta_list)

# write updated contributers list to google sheet, please make sure to check citations so looking ok, if update needed do it directly at the metasym level and push to github, then rerun code

write_sheet(meta_list, "https://docs.google.com/spreadsheets/d/1JIV92XisuBctBoVT53yqj0orPNRVyRQVEDtMq7R0fKw/edit#gid=0", sheet = "contributers")




