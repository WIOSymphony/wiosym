


#Script by gk & me Swedish WIO Symphony team
#Purpose: to update tool_metadata with data sources for each map
#R v 4.1.3
#updated 20230209 by gk

library(dplyr)
library(tidyverse)
library(stringi)
library(rio)
library(googlesheets4)


# reading directly from google sheets, ask for permission to enter from Wiosym admin

layer_list <-read_sheet("https://docs.google.com/spreadsheets/d/1JIV92XisuBctBoVT53yqj0orPNRVyRQVEDtMq7R0fKw/edit#gid=0", sheet = "tool_metadata")

# extract all sourcesym names used for final products

sourcesym_list <- layer_list %>% dplyr::select(title, theme, subtheme, sourcesym) %>% dplyr::filter(!is.na(sourcesym)) %>% unique()


# destinations
outfolder <- "./metadata/reg/v2/"



# list all metasym files and sourcesym files in wiosym cataloge

files_meta <- list.files(path = "./data_raw/", pattern = "_metasym.txt$", recursive=T)
files_source <- list.files(path = "./data/", pattern = "_sourcesym.txt$", recursive=T)
source_list <- tibble(title = c(), theme = c(), subtheme = c(), sourcesym= c(), data_sources = c())

# extract filepath for all sourcesym files used in final products

for(i in sourcesym_list$title){
  component <- sourcesym_list %>%  filter(title==i)
  file <- import(paste0("./data/", files_source[str_detect(files_source, pattern = component$sourcesym, negate = FALSE)]))
  metasym_id <- file$id
  tot <- tibble(title = component$title, theme = component$theme, subtheme = component$subtheme, sourcesym = component$sourcesym , data_sources = metasym_id)
  source_list <- rbind(source_list, tot)
}

source_list


# connect metasym id in source_list with citation and provider

providers <- readxl::read_excel("./shiny_data_upload/modify_txt_files_v01.2.1.xlsx", sheet = "providers_master")

citation_list <- tibble(data_sources = c(), provider_long = c(), citation = c())

for(j in source_list$data_sources){
  tryCatch({
    file <- read_tsv(paste0("./data_raw/", files_meta[str_detect(files_meta, pattern = paste0(j, "_metasym.txt$"), negate = FALSE)]))
    provider  <- ifelse(file[[5, 2]] == "", file[[6, 2]], file[[5, 2]])
    provider_long <- providers$provider[which(providers$provider_val == provider)]
    if(is_empty(provider_long)){
      provider_long <- ""
    }
    citation <- file[[8, 2]]
    
    tot <- tibble(data_sources = j, provider_long = provider_long, citation = citation)
    citation_list <- rbind(citation_list, tot)
  })
}

citation_list       # <- arrange(citation_list, data_sources, provider_long)

combined_list <- left_join(unique(citation_list), source_list, by = c("data_sources" = "data_sources"))

# collapse into unique rows for each component with citations and more seperated by ;

collapsed_list <- combined_list %>% group_by(title, theme, subtheme, sourcesym) %>% 
  dplyr::summarise(metasym_id=paste(data_sources, collapse="; "), providers_long=paste(provider_long, collapse="; "), `citation`=paste(citation, collapse="; "))

collapsed_list


# combine with original layer list in google sheet to ensure correct sorting, then write to sheet

updated_layer_list <- layer_list %>% left_join(collapsed_list) %>% 
  mutate(metasym=metasym_id) %>% 
  mutate(providers=providers_long) %>% 
  mutate(`data sources`=citation) #%>% 
  #select(!c(metasym_id, providers_long, citation))

View(updated_layer_list)

updated_layer_list_select <- updated_layer_list %>% select(c(providers, `data sources`, metasym))


# before write below, check all is safe and in order!

#range_write(data= updated_layer_list_select, "https://docs.google.com/spreadsheets/d/1JIV92XisuBctBoVT53yqj0orPNRVyRQVEDtMq7R0fKw/edit#gid=0", sheet = "tool_metadata", range = "O2", col_names=F)


#write_rds(combined_list, paste(outfolder, "combined_list.rds", sep=""))

