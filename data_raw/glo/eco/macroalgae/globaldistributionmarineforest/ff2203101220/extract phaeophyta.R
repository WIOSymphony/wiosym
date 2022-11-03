
library(readr)
library(dplyr)
databaseAllPruned <- read.csv("C:/Users/frafre/Downloads/databaseAllPruned.csv", header = FALSE, sep = ";", nrows = 1)
View(databaseAllPruned)
library(data.table)
x <- fread("C:/Users/frafre/Downloads/databaseAllPruned.csv", select = c("decimalLatitude", "decimalLongitude", "acceptedName", 
                                                                         "kingdom", "phylum", "class", "order"))
z <- filter(x, phylum == "Ochrophyta")