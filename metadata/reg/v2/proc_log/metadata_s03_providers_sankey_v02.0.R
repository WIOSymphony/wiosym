
#Script by me & gk Swedish WIO Symphony team
#Purpose: to visualize providers contribution to map layers in sankey diagram
#R v 4.1.3
#Updated: 20230209 by gk


# run metadata_s01 & metadata_s02 script prior to this code
#install.packages("devtools")

library(devtools)

#devtools::install_github("davidsjoberg/ggsankey")
#install.packages("tidyverse")

library(tidyverse)
library(ggsankey)
library(dplyr)
library(ggplot2)
library(rio)
library(googlesheets4)

# problems with purrr version... solved like this
#remove.packages("purrr")
#install.packages("purrr")
library(purrr)

#remove.packages("glue")
#install.packages("glue")
library(glue)

#install.packages("randomcoloR")

library(randomcoloR)


# destinations
outfolder <- "./metadata/reg/v2/"




## WIO SYmphony Sankey v2 -----------------------------------------------------------------

# read list of components and sourcedata from metadata_s02 script

combined_list <- read_rds(paste(outfolder, "combined_list.rds", sep=""))


# reading directly from google sheets, ask for permission to enter from Wiosym admin

#layer_list_metadata <-read_sheet("https://docs.google.com/spreadsheets/d/1JIV92XisuBctBoVT53yqj0orPNRVyRQVEDtMq7R0fKw/edit#gid=0", sheet = "tool_metadata")

#layer_list_providers <-read_sheet("https://docs.google.com/spreadsheets/d/1JIV92XisuBctBoVT53yqj0orPNRVyRQVEDtMq7R0fKw/edit#gid=0", sheet = "contributers")




## generate colors  group by ...
(sources <- combined_list %>% group_by(provider_long) %>% summarise())
(nr_sources <- sources %>% nrow())

(citations <- combined_list %>% group_by(citation) %>% summarise())
(nr_citations <- citations %>% nrow())



(eco_subtheme <- combined_list %>% group_by(theme, subtheme) %>% summarise() %>% filter(theme=="Ecosystem"))
(nr_eco_subtheme <- eco_subtheme %>% nrow())

(pres_subtheme <- combined_list %>% group_by(theme, subtheme) %>% summarise() %>% filter(theme=="Pressure"))
(nr_pres_subtheme <- pres_subtheme %>% nrow())

## generate colors

col_sources <- as_tibble(randomColor(count = nr_sources, hue = c("random")))
(col_sources <- bind_cols(sources, col_sources) %>% rename(hex_source = "value"))

col_eco_subtheme <- as_tibble(randomColor(count = nr_eco_subtheme, hue = c("green")))
(col_eco_subtheme <- bind_cols(eco_subtheme, col_eco_subtheme) %>% rename(hex = "value"))

col_pres_subtheme <- as_tibble(randomColor(count = nr_pres_subtheme, hue = c("red")))
col_pres_subtheme <- bind_cols(pres_subtheme, col_pres_subtheme) %>% rename(hex = "value")


col_subtheme <- bind_rows(col_eco_subtheme, col_pres_subtheme)


(combined_list_col <- combined_list %>% full_join(col_subtheme) %>% full_join(col_sources)) 

combined_list_col_order <- combined_list_col %>% group_by(theme, subtheme) %>% arrange(desc(title), .by_group = TRUE) %>% ungroup()


# set colors

named_col_source <- combined_list_col_order %>% select(provider_long, hex_source) %>% group_by(provider_long, hex_source) %>% summarise() 

(named_col_source <- named_col_source %>%  pull(hex_source, provider_long))


named_col<- combined_list_col_order %>% select(title, hex) %>% group_by(title, hex) %>% summarise() 

(named_col <- named_col %>%  pull(hex, title))


# attempt to order columns in sankey...
order_comp <- combined_list_col_order %>% group_by(theme, subtheme, title) %>% arrange(desc(theme), .by_group = TRUE) %>% summarise() %>% select(title) %>% pull()


combined_list_source_order <- combined_list_col %>% group_by(provider_long) %>% arrange(desc(citation), .by_group = TRUE) %>% summarise() %>% ungroup()

order_source <- combined_list_source_order %>% select(provider_long) %>% group_by(provider_long) %>% summarise()%>%  pull()



df <- combined_list_col_order |> distinct() |> make_long(provider_long, title) |> mutate(lab = ifelse(x=="provider_long", 1, -0.5)) 


order_source
order_comp

#df$node <- factor(df$node)
df$node <- factor(df$node,levels = c(order_comp, order_source))
#df$next_node <- factor(df$next_node)
df$next_node <- factor(df$next_node,levels = c(order_comp))

df$node
df$next_node


# https://ggplot2.tidyverse.org/reference/geom_text.html

ggplot(df,
       aes(x = x,
           next_x = next_x,
           node = node,
           next_node = next_node,
           fill = factor(node),
           # color = factor(node),
           label = node)
) +
  geom_alluvial(flow.alpha = .8, node.color = "transparent", width = 0.08) +
  geom_alluvial_text(size = 1.5, color = "black", hjust = "center", vjust = 0.2) +
  #geom_alluvial_text(size = 2.3, color = "black", hjust = 0.5) +
  #geom_alluvial_text(aes(label = ifelse(x == "provider_long", as.character(node), NA)), size = 1.8, color = "black", fill = "white", hjust = 2) +
  #geom_alluvial_text(aes(label = ifelse(x == "title", as.numeric(node), NA)), size = 1.3, color = "black", fill = "white", hjust = 0.5) +
  # scale_fill_viridis_d() +
  labs(x = NULL) +
  theme_void() +
  scale_fill_manual(values = c(named_col, named_col_source
  )) +
  theme(legend.position = "none")
  #geom_sankey()


ggsave(
  filename = paste(outfolder, "provider_sankey_v2.0.pdf", sep=""),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1.6,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

ggplot(df,
       aes(x = x,
           next_x = next_x,
           node = node,
           next_node = next_node,
           fill = factor(node),
           # color = factor(node),
           label = node)
) +
  geom_alluvial(flow.alpha = .8, node.color = "transparent", width = 0.08) +
  #geom_alluvial_text(size = 1.5, color = "black", hjust = "center", vjust = 0.2) +
  #geom_alluvial_text(size = 2.3, color = "black", hjust = 0.5) +
  geom_alluvial_text(aes(label = ifelse(x == "provider_long", as.character(node), NA)), size = 1.8, color = "black", fill = "white", hjust = 2) +
  #geom_alluvial_text(aes(label = ifelse(x == "title", as.numeric(node), NA)), size = 1.3, color = "black", fill = "white", hjust = 0.5) +
  # scale_fill_viridis_d() +
  labs(x = NULL) +
  theme_void() +
  scale_fill_manual(values = c(named_col, named_col_source
  )) +
  theme(legend.position = "none")
#geom_sankey()




ggsave(
  filename = paste(outfolder, "provider_sankey_v2.0_left_text.pdf", sep=""),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1.6,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



# alternativ för att bara spara text..!
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, label = node)) +
  geom_alluvial_text(size = 1.8, color = "black", hjust = "outward", vjust = 0.1) +
  theme_void()






