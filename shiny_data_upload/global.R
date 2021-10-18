# Import input files
folder_input <- read.table("folders.txt", sep = "\t", header = TRUE)
location_input <- read.table("locations.txt", sep = "\t", header = TRUE)
provider_input <- read.table("providers.txt", sep = "\t", header = TRUE)
tag_list <- read.table("tags.txt", sep = "\t", header = TRUE)


# Country list
locations <- as.vector(location_input$location_val)
names(locations) <- as.vector(location_input$location_name)

# Provider list
providers <- as.vector(provider_input$provider_val)
names(providers) <- as.vector(provider_input$provider_name)

# Copyright choices
copyright <- c("Unknown",
               "Open - product and raw data can be downloaded",
               "Partially open - product can be downloaded (resolution not finer than 250 m), raw data can be viewed, but not downloaded",
               "Restricted - product can be downloaded (resolution not finer than 250 m), raw data cannot be viewed or downloaded",
               "Restricted with special terms - product may be able to be downloaded under certain conditions (e.g. reduced resolution)"
               )

# Specify tag vectors
tag_act_names <- as.vector(tag_list$tag_act_name)[as.vector(tag_list$tag_act_name) != ""]
tag_act_vals <- as.vector(tag_list$tag_act_val)[as.vector(tag_list$tag_act_val) != ""]

tag_landact_names <- as.vector(tag_list$tag_landact_name)[as.vector(tag_list$tag_landact_name) != ""]
tag_landact_vals <- as.vector(tag_list$tag_landact_val)[as.vector(tag_list$tag_landact_val) != ""]

tag_adm_names <- as.vector(tag_list$tag_adm_name)[as.vector(tag_list$tag_adm_name) != ""]
tag_adm_vals <- as.vector(tag_list$tag_adm_val)[as.vector(tag_list$tag_adm_val) != ""]

tag_spec_names <- as.vector(tag_list$tag_spec_name)[as.vector(tag_list$tag_spec_name) != ""]
tag_spec_vals <- as.vector(tag_list$tag_spec_val)[as.vector(tag_list$tag_spec_val) != ""]

tag_hab_names <- as.vector(tag_list$tag_hab_name)[as.vector(tag_list$tag_hab_name) != ""]
tag_hab_vals <- as.vector(tag_list$tag_hab_val)[as.vector(tag_list$tag_hab_val) != ""]

tag_env_names <- as.vector(tag_list$tag_env_name)[as.vector(tag_list$tag_env_name) != ""]
tag_env_vals <- as.vector(tag_list$tag_env_val)[as.vector(tag_list$tag_env_val) != ""]

tag_pres_names <- as.vector(tag_list$tag_pres_name)[as.vector(tag_list$tag_pres_name) != ""]
tag_pres_vals <- as.vector(tag_list$tag_pres_val)[as.vector(tag_list$tag_pres_val) != ""]

tag_func_names <- as.vector(tag_list$tag_func_name)[as.vector(tag_list$tag_func_name) != ""]
tag_func_vals <- as.vector(tag_list$tag_func_val)[as.vector(tag_list$tag_func_val) != ""]


# Specifying folders
# Theme folders

folder_theme <- as.vector(unique(folder_input$theme_folder))
names(folder_theme) <- as.vector(unique(folder_input$theme_name))

# Sub-theme folders

folder_act <- as.vector(folder_input$subtheme_folder[which(folder_input$theme_folder == "act")])
names(folder_act) <- as.vector(folder_input$subtheme_name[which(folder_input$theme_folder == "act")])
folder_act

folder_adm <- as.vector(folder_input$subtheme_folder[which(folder_input$theme_folder == "adm")])
names(folder_adm) <- as.vector(folder_input$subtheme_name[which(folder_input$theme_folder == "adm")])
folder_adm

folder_eco <- as.vector(folder_input$subtheme_folder[which(folder_input$theme_folder == "eco")])
names(folder_eco) <- as.vector(folder_input$subtheme_name[which(folder_input$theme_folder == "eco")])
folder_eco

folder_env <- as.vector(folder_input$subtheme_folder[which(folder_input$theme_folder == "env")])
names(folder_env) <- as.vector(folder_input$subtheme_name[which(folder_input$theme_folder == "env")])
folder_env

folder_pres <- as.vector(folder_input$subtheme_folder[which(folder_input$theme_folder == "pres")])
names(folder_pres) <- as.vector(folder_input$subtheme_name[which(folder_input$theme_folder == "pres")])
folder_pres

folder_unknown <- as.vector(folder_input$subtheme_folder[which(folder_input$theme_folder == "unknown")])
names(folder_unknown) <- as.vector(folder_input$subtheme_name[which(folder_input$theme_folder == "unknown")])
folder_unknown

theme1_list <- data.frame("names" = names(folder_theme), "values" = folder_theme)

theme2_list_temp <- c(folder_act, folder_adm, folder_eco, folder_env, folder_pres, folder_unknown)

theme2_list <- data.frame("names" = names(theme2_list_temp), "values" = theme2_list_temp)


