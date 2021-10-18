# ========================================================================================================= #
#### About ####
# App that creates a unique metadata file (idyyyymmddhhmm_metasym.txt) for any data used in the WIO Symphony project, 
# which will make sure the incomming data is tracked in the processing of any products
# the App also creates the folder organisation for incomming data

#### Instructions ####
# Run the line below code without "#" if you have not yet installed the shiny package (e.g. install.packages("shiny")), 
# then click "> Run App" (top right corner of this window if you opened this file in R studio). 
# Follow the instructions provided in the app, once completed you need to move your data to the created folder. 
# Consult the instructional video and the Swedish WIO Symphony team if you have questions  
# ========================================================================================================= #

# Shiny App for metadata and folder structure for incomming data into WIO SYmphony --------------------------

#install.packages("shiny")
library(shiny)
source("global.R")

ui <- fluidPage(
  
  # Customization of UI
  
  tags$style(HTML("#head1{font-family:Arial; font-size: 22px; font-weight: bold; padding-bottom: 0px; margin-bottom: 3px; margin-top: 3px}
                   #head2{font-family:Arial; font-size: 20px; font-style: italic; padding-bottom: 0px; margin-bottom: 20px; margin-top: 3px}
                   #head3{font-family:Arial; font-size: 18px; font-weight: bold; padding-bottom: 0px; margin-bottom: 3px; margin-top: 3px}
                   #head4{font-family:Arial; font-size: 14px; font-style: italic; padding-bottom: 0px; margin-bottom: 3px; margin-top: 3px}
                  ")),
  
  # Custom CSS for 4 column tag layout
  
  tags$head(
    tags$style(HTML("
                     .multicol .shiny-options-group{
                     -webkit-column-count: 4; /* Chrome, Safari, Opera */
                     -moz-column-count: 4;    /* Firefox */
                     column-count: 4;
                     -moz-column-fill: balanced;
                     -column-fill: balanced;
                     }
                     .checkbox{
                     margin-top: 0px !important;
                     -webkit-margin-after: 0px !important; 
                     }
                     "
              )
    )
  ),
  
  # Submit button customization
  
  tags$head(
    tags$style(HTML("#save{background-color:PowderBlue; font-weight: bold; padding-bottom: 10px; margin-top: 12px; margin-bottom: 12px}
                    "))
  ),
  
  # Select input drop-down list customization
  
  tags$style(type='text/css',
             ".selectize-dropdown-content{max-height: 500px;}"),
  
  # Check-box customization
  
  tags$style("#tag_act {font-size:14px;height:300px;}
              #tag_landact {font-size:14px;height:300px;}
              #tag_adm {font-size:14px;height:300px;}
              #tag_spec {font-size:14px;height:300px;}
              #tag_hab {font-size:14px;height:300px;}
              #tag_env {font-size:14px;height:300px;}
              #tag_pres {font-size:14px;height:300px;}
              #tag_func {font-size:14px;height:300px;}
             "),
  
  titlePanel(img(src='logo.jpg', align = "right")
  ),
  
  fluidRow(
    column(width = 7, align = "left",
           
      h1(id = "head1", "WIOSym data upload"),     
      h1(id = "head2", "Please fill in the details of the dataset below. Once you have done that click 'Submit'. 
         Then you can add your data to the folder specified below."),
      h1(id="head2", "* required input"),
      
      # Load existing metadata
      
      h1(id = "head3" ,"Import existing metadata"),
      h1(id = "head4" ,"Here you can import existing metadata. Import the id_metasym.txt file associated with the data to load the metadata 
         into the shiny app."),
      fileInput(inputId = "import",
                label = NULL,
                multiple = FALSE),
      
      # Display existing metadata
      tableOutput("metatable"),
      
      
      # User
      
      h1(id = "head3" ,"Uploaded by user *"),
      h1(id = "head4" ,"Please provide your initials (e.g. 'es') - do not write your full name"),
      
      textInput('user', 
                label = NULL, 
                value = "",  width = "100%"),
      
      # Provider
      
      h1(id = "head3" ,"Data provider *"),
      h1(id = "head4" ,"If the data provider is not on the list, you may input the provider name manually.
                        If you use manual input for the provider, please use an abbreviated version of the provider name. 
                        For example, if the provider is Global Fishing Watch, you could put 'gfw' as the provider. You may add further information
                        about the provider in the 'Additional comments' section below."),

      
      tabsetPanel(
        
        tabPanel(title = "Select from list",
                 selectInput(inputId = 'provider2', 
                             label = NULL,
                             choices = providers,
                             selected = "No provider selected",
                             width = "100%")),
      
        tabPanel(title = "Manual input",
                 textInput(inputId = 'provider1', 
                           label = NULL, 
                           value = "", width = "100%"))
                 
      ), 

      # Source
      
      h1(id = "head3" ,"Data source (e.g. URL or DOI) *"),
      h1(id = "head4" ," If the dataset was obtained via personal communication, 
         please write the email of the person you obtained it from"),
      h1(id = "head4" , "Examples: https://globalfishingwatch.org/, https://doi.org/10.1111/csp2.221, johnsmith@gmail.com"),
      
      textInput('source', 
                label = NULL, 
                value = "",  width = "100%"),
      
      # Copyright
      
      h1(id = "head3" ,"Copyright *"),
      h1(id = "head4" ,"Please choose the level of copyright for the dataset. If you are unsure, or need to confirm the restrictions at a later date,
                        please select 'Unknown'"),
      
      selectInput(inputId = "copyright", 
                  label = NULL,
                  choices = copyright,
                  selected = "Unknown",
                  width = "100%"),
      
      # Copyright details
      
      h1(id = "head3" ,"Copyright details"),
      h1(id = "head4" ,"Please add any relevant comments concerning copyrights (e.g. details of special copyright terms)"),
      
      textAreaInput(inputId = "copyrightdet",
                    label = NULL,
                    width = "800px",
                    height = "80px",
                    resize = "none"),
      # Citation
      
      h1(id = "head3" ,"Citation"),
      h1(id = "head4" ,"If required, please provide customized citations for this dataset"),
      h1(id = "head4" , "Example: UNEP-WCMC & IUCN (2019). The World Database on Protected Areas (WDPA). www.protectedplanet.net."),
      
      textInput("citation", 
                label = NULL, 
                value = "",  width = "100%"),
      
      # Location
      
      h1(id = "head3" ,"Location *"),
      h1(id = "head4" ,"Please choose a source location for the dataset - Choose 'Regional' for international/regional datasets 
         in the WIO, country codes for national data, and 'Global' for datasets with a full global extent"),
      
      selectInput(inputId = "location", 
                  label = NULL,
                  choices = locations,
                  selected = "Regional",
                  width = "100%"),
      
      # Theme 1
      
      h1(id = "head3" ,"Theme *"),
      h1(id = "head4" ,"Choose the most relevant theme and sub-theme for the data.
                        If the dataset contains data relevant to multiple themes, choose one that is the most appropriate,
                        and specify other relevant themes for the dataset in the tags section below"),
      
      selectInput(inputId = "theme1", 
                  label = NULL,
                  choices = folder_theme, 
                  selected = "unknown",
                  width = "100%"),
      
      # Theme 2
      
      h1(id = "head3" ,"Sub-theme *"),
      
      selectInput(inputId = "theme2", 
                  label = NULL,
                  choices = folder_unknown,
                  width = "100%"),
      
      h1(id = "head3" ,"Tags"),
      h1(id = "head4" ,"Please choose any tags relevant to this dataset.
                        You may choose any tags regardless of the theme. E.g. if a dataset falls under the 'Activity' theme,
                        you can still choose species or habitat tags that may be relevant to the dataset."),
      h1(id = "head4" ,"Example: A dataset of tuna fisheries catch could be tagged with 'commercial fishing' and 'pelagic fish'"),
      
      tabsetPanel(
        tabPanel(title = "Administrative", 
                 tags$div(class = "multicol",
                          checkboxGroupInput(inputId = "tag_adm", label = "",
                                             choices = c("No choices yet"),
                                             inline = TRUE,
                                             width = "100%"))),
        tabPanel(title = "Environment", 
                 tags$div(class = "multicol",
                          checkboxGroupInput(inputId = "tag_env", label = "",
                                             choices = c("No choices yet"),
                                             inline = TRUE,
                                             width = "100%"))),
        tabPanel(title = "Activities", 
                 tags$div(class = "multicol",
                          checkboxGroupInput(inputId = "tag_act", label = "",
                                             inline = TRUE,
                                             width = "100%"))),
        tabPanel(title = "Land activities", 
                 tags$div(class = "multicol",
                          checkboxGroupInput(inputId = "tag_landact", label = "",
                                             inline = TRUE,
                                             width = "100%"))),
        tabPanel(title = "Pressures", 
                 tags$div(class = "multicol",
                          checkboxGroupInput(inputId = "tag_pres", label = "",
                                             choices = c("No choices yet"),
                                             inline = TRUE,
                                             width = "100%"))),
        tabPanel(title = "Habitats", 
                 tags$div(class = "multicol",
                          checkboxGroupInput(inputId = "tag_hab", label = "",
                                             choices = c("No choices yet"),
                                             inline = TRUE,
                                             width = "100%"))),
        tabPanel(title = "Species", 
                 tags$div(class = "multicol",
                          checkboxGroupInput(inputId = "tag_spec", label = "",
                                             choices = c("No choices yet"),
                                             inline = TRUE,
                                             width = "100%"))),
        tabPanel(title = "Functions", 
                 tags$div(class = "multicol",
                          checkboxGroupInput(inputId = "tag_func", label = "",
                                             choices = c("No choices yet"),
                                             inline = TRUE,
                                             width = "100%"))),
        selected = "Habitats"
      ),
      
      # Comments
      
      h1(id = "head3" ,"Additional comments"),
      
      textAreaInput(inputId = "comment",
                    label = NULL,
                    width = "800px",
                    height = "100px",
                    resize = "none"),
      
      
      htmlOutput("directory"),
      htmlOutput("localpath"),
      
      verbatimTextOutput("submit"),
      
      actionButton("save", "Submit", width = "50%")
      
    )
  )
)

server <- function(input , output, session){
    options(shiny.maxRequestSize=10000*1024^2)
  
    # Add reactive element for flexible variables
    flex <- reactiveValues()
    flex$theme_tag <- "NA"
    flex$subtheme_tag <- "NA"
    
    # Convert current time to number string, removing first two year digits, and removing last two second digits
    # Format is yymmddhhmm
    
    syst <- Sys.time()
    flex$downdate <- substr(syst, start = 1, stop = 10)
    t <- gsub("\\D", "", syst)
    t <- substring(t, 3, nchar(t)-2)
    
    
    # Determining folder name and id, and updating inputs if metadata is imported
    
    observe({
      
      user <- input$user
      user <- gsub(pattern = "[[:punct:]]", replacement = "", user)
      user <- gsub(pattern = "[[:blank:]]", replacement = "", user)
      user <- gsub(pattern = "[[:digit:]]", replacement = "", user)
      user <- gsub(pattern = "([[:upper:]])", perl = TRUE, replacement = "\\L\\1", user)
      
      flex$user <- user
      
      
      # If no metadata is uploaded
      if(is.null(input$import)) {
        
        flex$id <- paste0(user, t)
        
        updateCheckboxGroupInput(
          session,
          "tag_act",
          choiceNames = tag_act_names,
          choiceValues = tag_act_vals,
          selected = NULL
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_landact",
          choiceNames = tag_landact_names,
          choiceValues = tag_landact_vals,
          selected = NULL
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_adm",
          choiceNames = tag_adm_names,
          choiceValues = tag_adm_vals,
          selected = NULL
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_spec",
          choiceNames = tag_spec_names,
          choiceValues = tag_spec_vals,
          selected = NULL
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_hab",
          choiceNames = tag_hab_names,
          choiceValues = tag_hab_vals,
          selected = NULL
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_env",
          choiceNames = tag_env_names,
          choiceValues = tag_env_vals,
          selected = NULL
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_pres",
          choiceNames = tag_pres_names,
          choiceValues = tag_pres_vals,
          selected = NULL
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_func",
          choiceNames = tag_func_names,
          choiceValues = tag_func_vals,
          selected = NULL
        )
        
        updateSelectInput(
          session,
          "theme2",
          choices = if (input$theme1 == "act") {folder_act}
          else if (input$theme1 == "adm") {folder_adm}
          else if (input$theme1 == "eco") {folder_eco}
          else if (input$theme1 == "env") {folder_env}
          else if (input$theme1 == "pres") {folder_pres}
          else if (input$theme1 == "unknown") {folder_unknown}
        )
        
        }
      
      # If metadata is uploaded
      else if(!is.null(input$import)){
        
        # Storing the metadata variables in a data frame
        imp.file <- input$import$datapath
        imp.temp <- read.table(imp.file, sep = "\t")
        imp <- data.frame(matrix(ncol = 18, nrow = 1)) 
        colnames(imp) <- imp.temp[,1]
        imp[1,] <- imp.temp[,2]
        colnames(imp.temp) <- NULL
        output$metatable <- renderTable(imp.temp)
        
        # Updating the inputs
        flex$downdate <- imp$download_date
        flex$theme_tag <- imp$theme_tag
        flex$subtheme_tag <- imp$subtheme_tag
        flex$id <- imp$id                                                       # ID update
        updateTextInput(session, "user", value = imp$user)                      # User update
        updateSelectInput(session, "provider2", selected = imp$provider_list)   # Provider from list update
        updateSelectInput(session, "provider1", selected = imp$provider_manual) # Provider manual update
        updateTextInput(session, "source", value = imp$source)                  # Source update
        updateSelectInput(session, "copyright", selected = imp$copyright)       # Copyright update
        updateTextInput(session, "copyrightdet", value = imp$copyright_details) # Copyright details update
        updateTextInput(session, "citation", value = imp$citation)              # Citation update
        updateSelectInput(session, "location", selected = imp$location)         # Location update
        updateSelectInput(session, "theme1", selected = imp$theme_tag)          # Theme update (Sub-theme is updated below)
        updateTextInput(session, "comment", value = imp$comments)               # Comments update
        
        # Tag update
        tag_update <- strsplit(imp$tags, ";")
        tag_update <- tag_update[[1]]
        
        updateCheckboxGroupInput(
          session,
          "tag_act",
          choiceNames = tag_act_names,
          choiceValues = tag_act_vals,
          selected = tag_update
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_landact",
          choiceNames = tag_landact_names,
          choiceValues = tag_landact_vals,
          selected = tag_update
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_adm",
          choiceNames = tag_adm_names,
          choiceValues = tag_adm_vals,
          selected = tag_update
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_spec",
          choiceNames = tag_spec_names,
          choiceValues = tag_spec_vals,
          selected = tag_update
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_hab",
          choiceNames = tag_hab_names,
          choiceValues = tag_hab_vals,
          selected = tag_update
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_env",
          choiceNames = tag_env_names,
          choiceValues = tag_env_vals,
          selected = tag_update
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_pres",
          choiceNames = tag_pres_names,
          choiceValues = tag_pres_vals,
          selected = tag_update
        )
        
        updateCheckboxGroupInput(
          session,
          "tag_func",
          choiceNames = tag_func_names,
          choiceValues = tag_func_vals,
          selected = tag_update
        )
        
        updateSelectInput(
          session,
          "theme2",
          choices = if (imp$theme_tag == "act") {folder_act}
          else if (imp$theme_tag == "adm") {folder_adm}
          else if (imp$theme_tag == "eco") {folder_eco}
          else if (imp$theme_tag == "env") {folder_env}
          else if (imp$theme_tag == "pres") {folder_pres}
          else if (imp$theme_tag == "unknown") {folder_unknown},
          selected = imp$subtheme_tag
        )
        
      }
      
    })
  
    
    observeEvent(input$theme1, {

      if(!is.null(input$import) & input$theme1 == flex$theme_tag){
        
         updateSelectInput(
          session,
          "theme2",
          choices = if (input$theme1 == "act") {folder_act}
          else if (input$theme1 == "adm") {folder_adm}
          else if (input$theme1 == "eco") {folder_eco}
          else if (input$theme1 == "env") {folder_env}
          else if (input$theme1 == "pres") {folder_pres}
          else if (input$theme1 == "unknown") {folder_unknown},
          selected = flex$subtheme_tag
        )
        
      }
      
      else if (!is.null(input$import) & input$theme1 != flex$theme_tag){
        
        updateSelectInput(
          session,
          "theme2",
          choices = if (input$theme1 == "act") {folder_act}
          else if (input$theme1 == "adm") {folder_adm}
          else if (input$theme1 == "eco") {folder_eco}
          else if (input$theme1 == "env") {folder_env}
          else if (input$theme1 == "pres") {folder_pres}
          else if (input$theme1 == "unknown") {folder_unknown}
        )
        
      }
      
    }
    )
    
    
    # Determining folder location
  
    dest_dir <- reactiveValues(
      dir = "data_raw/quarantine/" 
    )
    
    
    observe({
      
      prov1 <- input$provider1
      prov1 <- gsub(pattern = "[[:punct:]]", replacement = "", prov1)
      prov1 <- gsub(pattern = "[[:blank:]]", replacement = "", prov1)
      prov1 <- gsub(pattern = "[[:digit:]]", replacement = "", prov1)
      prov1 <- gsub(pattern = "([[:upper:]])", perl = TRUE, replacement = "\\L\\1", prov1)
      
      prov2 <- input$provider2
      prov2 <- gsub(pattern = "[[:punct:]]", replacement = "", prov2)
      prov2 <- gsub(pattern = "[[:blank:]]", replacement = "", prov2)
      prov2 <- gsub(pattern = "[[:digit:]]", replacement = "", prov2)
      prov2 <- gsub(pattern = "([[:upper:]])", perl = TRUE, replacement = "\\L\\1", prov2)
      
      if(input$provider1 == "" & input$provider2 == "none") {flex$prov <- "unknown"}
      else if(input$provider1 != "") {flex$prov <- prov1}
      else if(input$provider1 == "" & input$provider2 != "none") {flex$prov <- prov2}
      
      if(input$location == "reg") {
        
        if(input$theme1 == "unknown") {
          
          dest_dir$dir <- paste("data_raw/reg/quarantine", flex$prov, flex$id, sep = "/")
          
        }
        
        else {
          
          dest_dir$dir <- paste("data_raw/reg", input$theme1, input$theme2, flex$prov, flex$id, sep = "/")
          
        }
        
      } 
      
      else if(input$location == "glo") {
        
        if(input$theme1 == "unknown") {
          
          dest_dir$dir <- paste("data_raw/glo/quarantine", flex$prov, flex$id, sep = "/")
          
        }
        
        else {
          
          dest_dir$dir <- paste("data_raw/glo", input$theme1, input$theme2, flex$prov, flex$id, sep = "/")
          
        }
        
      } 
      
      else if(input$location != "reg" & input$location != "glo") {
        
        if(input$theme1 == "unknown") {
          
          dest_dir$dir <- paste("data_raw/nat", input$location, "quarantine", flex$prov, flex$id, sep = "/")
          
        }
        
        else {

          dest_dir$dir <- paste("data_raw/nat", input$location, input$theme1, input$theme2, flex$prov, flex$id, sep = "/")
          
        }
        
      } 
      
        # output$directory <- renderText({
        #   HTML(paste0("Adding metadata to: <b>", dest_dir$dir,"</b><br/>"))
        # })
      
        output$localpath <- renderText({
          p <- gsub("/", "\\\\", paste0(dirname(getwd()),"/", dest_dir$dir))
          HTML(paste0("Adding metadata to: <b>", p, "</b><br/>"))
        })
    })
    

    # Writing the metadata file and saving tags
    
    metad <- reactive({
      
      sub_tag <- as.character(folder_input$subtheme_tag[folder_input$subtheme_folder == input$theme2])
      
      theme_tag <- ifelse(input$theme1 == "unknown", "", paste(input$theme1, sub_tag, sep = ";"))
      
      tag_result <- paste(";",
                          theme_tag,
                          paste(input$tag_adm, collapse = ";"),
                          paste(input$tag_env, collapse = ";"),
                          paste(input$tag_act, collapse = ";"),
                          paste(input$tag_landact, collapse = ";"),
                          paste(input$tag_pres, collapse = ";"),
                          paste(input$tag_hab, collapse = ";"),
                          paste(input$tag_spec, collapse = ";"),
                          paste(input$tag_func, collapse = ";"),
                          ";",
                          sep = ";")
      
      tag_result <- gsub('(\\;)\\1+', '\\1', tag_result)
      
      n <- c("id", 
             "file_path", 
             "local_path",
             "download_date",
             "user", 
             "provider_list",
             "provider_manual",
             "source", 
             "citation", 
             "copyright",
             "copyright_details",
             "location",
             "theme",
             "theme_tag",
             "subtheme",
             "subtheme_tag",
             "tags", 
             "comments")
      v <- c(paste0(flex$id), 
             gsub("/", "\\\\", paste0(dest_dir$dir)), # short file path
             gsub("/", "\\\\", paste0(dirname(getwd()))), # local file path
             flex$downdate,
             flex$user,
             input$provider2,
             input$provider1,
             input$source, 
             input$citation, 
             input$copyright,
             gsub('"', "'", input$copyrightdet, fixed=TRUE),
             input$location,
             paste(theme1_list$names[theme1_list$values == input$theme1]),
             input$theme1,
             paste(theme2_list$names[theme2_list$values == input$theme2]),
             input$theme2,
             tag_result, 
             gsub('"', "'", input$comment, fixed=TRUE))
      
      data.frame("Col1" = n, "Col2" = v)
      
    })
    
      
    # Submission button output
    
    observeEvent(input$save, {
      
      data_dir <- paste("..", dest_dir$dir, sep = "/")
      
      dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
      
      # Use this code to create a new updated metadata file if existing metadata is uploaded
      #if(!is.null(input$import)){file_name <- paste0(flex$id, "_metasym_update.txt")}
      #else {file_name <- paste0(flex$id, "_metasym.txt")}
      #meta.file <- paste0(data_dir, "/", file_name)
      
      # Use this code to overwrite old metadata file if existing metadata is uploaded
      file_name <- paste0(flex$id, "_metasym.txt")
      meta.file <- paste0(data_dir, "/", file_name)

      write.table(metad(), file = meta.file, sep = "\t", col.names = FALSE, row.names = FALSE)
      
      output$submit <- renderPrint({
        
        cat("Success! The metadata has been added to the directory specified above, in the file '", file_name, "'\n")
        cat("Now, please add the raw data to the '", flex$id, "' folder")
        
      })
      
      
    })   
    
}

shinyApp(ui = ui , server = server)
