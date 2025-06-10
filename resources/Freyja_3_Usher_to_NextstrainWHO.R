rm(list=ls(all=FALSE))
ls()
character(0)

#library(data.table)
require(tidyverse)
#library(jsonlite)



clade_display_names<-read.table("clade_display_names.tsv", header = TRUE, sep = "\t")
alias_key<-read.table( "alias_key.tsv", header = TRUE, sep = "\t")
database_lineages_table<-read.table( "outbreakinfo_db_table_melted.tsv", header = TRUE, sep = "\t")
manual_lineages_table<-read.table( "manual_lineages_table.tsv", header = TRUE, sep = "\t")

# Prepare clade_display_names
clade_display_names<-clade_display_names %>%
  dplyr::arrange(-desc(Nextstrain)) %>%
  dplyr::mutate(legend.position = 1:nrow(clade_display_names))

# Prepare database_lineages_table2 = database_lineages_table + manual_lineages_table
database_lineages_table2<-manual_lineages_table %>%
  dplyr::add_row(database_lineages_table)

# Merge clade_display_names with database_lineages_table2
clade_display_names_WHO <- dplyr::left_join(clade_display_names, database_lineages_table2, by="Pango") %>%
  replace(is.na(.), "") %>%
  dplyr::filter(WHO!="") %>% # Get rid of row which have no WHO name
  dplyr::mutate(Pango = ifelse(Pango_replace != "", Pango_replace, Pango), #allow to change Pango (e.g., XDV.1 > XDV)
                Pango = ifelse(Pango == WHO, "", Pango), # Empty Pango if WHO = Pango column (indicate that Pango has a WHO name)
                Legend = ifelse(Pango!="", paste0(WHO, " (", Pango, ")"), WHO)) # Prepare legend for dashboard


write.table(clade_display_names_WHO, "ListVariantToDisplay.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")





