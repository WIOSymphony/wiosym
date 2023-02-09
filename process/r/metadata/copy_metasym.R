a <- list.files(path = "M:/marin/swoc/work/wiosym/data_raw/glo/", all.files = TRUE, full.names = TRUE, pattern = "_metasym.txt", recursive = TRUE)
file.copy(from = a, to = "M:/marin/swoc/work/wiosym/metadata/metasym_copy/", overwrite = TRUE)
