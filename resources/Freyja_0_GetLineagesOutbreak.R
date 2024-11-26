rm(list=ls(all=FALSE))
ls()
set.seed(123)
library(outbreakinfo) #devtools::install_github("outbreak-info/R-outbreak-info")
require(dplyr)
require(stringr)
require(reshape2)

##########################################################################################################
#### Extract info about lineages #########################################################################
##########################################################################################################

# Extract variants
curated = getCuratedLineages()

### Retrieve Outbreak.info data (all variants)
all.variant = curated %>% 
  select(who_name, char_muts_parent_query, variantType)
all.variant<-all.variant %>% filter(!grepl("[*]", who_name))


### Reformat all.variant 1/2
all.variant1<-all.variant
all.variant1$char_muts_parent_query<-gsub(" OR ", ", ", all.variant1$char_muts_parent_query)
names(all.variant1)<-c("WHO", "Pango sublineages", "Status")


### Reformat all.variant 2/2
all.variant2<-all.variant1 %>%
  dplyr::select(-Status)

all.variant2<-as.data.frame(cbind(all.variant2[,1], stringr::str_split_fixed(all.variant2$`Pango sublineages`, ',', max(stringr::str_count(as.character(all.variant2$`Pango sublineages`), ",")+1)))) # Split all Pango lineage is distinct columns

all.variant2<-reshape2::melt(all.variant2, id=c("V1")) %>%
  dplyr::mutate(WHO = V1,
                Pango = str_trim(value)) %>%  # Remove spaces before Pango
  dplyr::filter(Pango != "") %>% # Remove rows for which Pango = NULL
  dplyr::select(WHO, Pango)
  


##########################################################################################################
#### Export ##############################################################################################
##########################################################################################################
write.table(all.variant1, "outbreakinfo_db_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(all.variant2, "outbreakinfo_db_table_melted.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



