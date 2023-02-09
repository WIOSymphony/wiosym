library(dplyr)
library(stringi)

sourcesym <- list.files("./././data", recursive = TRUE, pattern = "sourcesym.txt$", full.names = T)
metasym <- list.files("./././", recursive = TRUE, pattern = "metasym.txt$", full.names = T)

metasymID <- list()
citation <- data.frame(name  = c(), cit = c())

for(i in (1:length(sourcesym))){
  metasymID[[i]] <- read.table(sourcesym[i], header = T)
  names(metasymID[[i]]) <- sourcesym[i]
  cit <- ""
  for(j in 1:nrow(metasymID[[i]])){
    mID <- which(grepl(metasymID[[i]][j, 1], metasym, fixed=TRUE))
    a <- read.table(metasym[mID[1]], header = T, sep = "\t")
    b <- ifelse(a[8,2]=="", "EMPTY SOURCE", a[8,2])
    cit <- paste0(cit, b, "; ")
  }
  citation[i,1] <- sourcesym[i]
  citation[i,2] <- cit
}

citation |> show_in_excel()



spl <- strsplit(as.character(sourcesym), "/")
sourcesym_short <- sapply(lapply(spl, tail, 1), paste, collapse="_")

# sSymDoc <- read.table("clipboard", sep = "\t")
proc <- list.files("./././data", recursive = TRUE, pattern = ".Rmd$", full.names = T)

sID <- sourcesym[which(sourcesym_short %in% sSymDoc$V1)]
spl <- strsplit(as.character(sID), "/")
sourcesym_short <- sapply(lapply(spl, tail, 1), paste, collapse="_")

metasymID <- list()
citation <- data.frame(name  = c(), cit = c())
for(i in (1:length(sID))){
  metasymID[[i]] <- read.table(sID[i], header = T)
  names(metasymID[[i]]) <- sID[i]
  cit <- ""
  meta <- ""
  for(j in 1:nrow(metasymID[[i]])){
    mID <- which(grepl(metasymID[[i]][j, 1], metasym, fixed=TRUE))
    a <- read.table(metasym[mID[1]], header = T, sep = "\t")
    b <- ifelse(a[8,2]=="", "EMPTY SOURCE", a[8,2])
    cit <- paste0(b, ";", cit)
    meta <- paste0(metasymID[[i]][j, 1], ";", meta)
  }
  citation[i,1] <- sourcesym_short[i]
  citation[i,2] <- cit
  citation[i,3] <- substring(meta, 1, nchar(meta)-1)
}

p <-paste0(sapply(lapply(spl, head, -1), paste, collapse="/"), "/proc_log/")

method <- c()
for(i in 1:length(p)){
  list.files(p[i])
  method[i] <- paste0(list.files(p[i]), collapse = ";")
}





p1 <- "C:/Users/mner0007/OneDrive - Sveriges lantbruksuniversitet/wiosym_slu_internal/shiny_review/raster1km"

l <- list.files(p1, full.names = T)

a <- data.frame(name = c(), min = c(), max = c())
for(i in 1:length(l)){
  r <- raster(l[i])
  a <- rbind(a, data.frame(l1[i], values(r) |> na.omit() |> min(), values(r) |> na.omit() |> max()))
}