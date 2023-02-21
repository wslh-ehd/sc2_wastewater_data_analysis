rm(list=ls(all=FALSE))
ls()
set.seed(123)
library(outbreakinfo) #devtools::install_github("outbreak-info/R-outbreak-info")
library(dplyr)


##########################################################################################################
#### Extract info about lineages #########################################################################
##########################################################################################################

# Extract variants
curated = getCuratedLineages()

### Retrieve Outbreak.info data (all variants)
all.variant = curated %>% 
  select(who_name, char_muts_parent_query, variantType)
all.variant<-all.variant %>% filter(!grepl("[*]", who_name))

### Reformat all.variant
all.variant$char_muts_parent_query<-gsub(" OR ", ", ", all.variant$char_muts_parent_query)
names(all.variant)<-c("WHO", "Pango sublineages", "Status")


##########################################################################################################
#### Export ##############################################################################################
##########################################################################################################
write.table(all.variant, "database_lineages_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


