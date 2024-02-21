library(xml2)

ver <- "2.1"
dir <- paste0("../products/v", ver, "/output_geotiffs")
l <- list.files(path = dir, pattern = ".json", recursive = T, full.names = T)
l.short <- list.files(path = dir, pattern = ".json", recursive = T, full.names = F)

root <- xml_new_root("urlset", ns = "https://www.sitemaps.org/schemas/sitemap/0.9/")

for(i in 1:length(l)){
  turl <- paste0("https://raw.githubusercontent.com/WIOSymphony/wiosym/main/products/v2.1/output_geotiffs/", l.short[i])
  t <- file.info(l[i])$mtime
  t <- as.Date(t)
  t <- paste0(t)
  
  url_node <- xml_add_child(root, "url")
  xml_add_child(url_node, "loc", turl)
  xml_add_child(url_node, "lastmod", t)
}

doc <- as_xml_document(root)
#xml_declaration <- "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
#xml_doc_content <- paste0(xml_declaration, as.character(doc))

writeLines(xml_doc_content, "sitemap.xml")



