SHINY APP READ ME

This folder contains the shiny app for uploading raw data to the WIOSym database.
To execute the app you need to have R and RStudio installed. To run the app open the "shiny_data_upload.R" file in RStudio. 
At the top right of the script window there will be a button that says "Run App". Press this to run the app.

Before running the app, you need to install the R shiny package with the code: 

	install.packages("shiny")

The app and folder are designed so that the "shiny_data_upload" folder containing everything should be in the same location
as the "data_raw" folder, i.e. they should be on the same level (NOT inside the data_raw folder!!). If you have not created a data_raw folder, the shiny app wil create one automatically.

OTHER FILES:

The other files in the folder are inputs to the shiny app that specify the locations/countries, themes, and providers that are options in the shiny app.

The "www" folder is just a folder containing the WIOSym logo that is displayed in the shiny app.

The "global.R" file is a script that shiny executes before running, so any variables/parameters can be specified in this script without
clogging up the main shiny script.
