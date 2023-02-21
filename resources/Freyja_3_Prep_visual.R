rm(list=ls(all=FALSE))
ls()
set.seed(123)


#################################################################################
#setwd("/scratch/projects/SARS-CoV-2/Results/2022-10-20_removeRecombinant/freyja/")

min.proportion = 1 # 0-100

# WHO lineages with sublineages to display on the dashboard 
list.variants.with.sublineages<-c("Omicron", "Add lineage here")

# List sublineages
replace.lineages<-data.frame(Pango=c("B.1.1.529", "Add lineage here"),
                             WHO=c("Omicron", "Add lineage here"))

# Complete sublineages that are not listed Nextstrain
#complete.taxo<-data.frame(WHO=c("Omicron", "TESTOmicronTEST", "Add lineage here"),
#                          Pango=c("BA.3", "TESTBA.4.6TEST", "Add lineage here"))

# When simplify alias_key.json ignore the sublineages below
ignore.name.lineages<-data.frame(Pango.1=c("XBB", "BexampleA", "Add lineage here"),
                                 Pango.2=c("XBB", "Bexample.1.1.529", "Add lineage here"))


#################################################################################



library(tidyverse)
library(strex)
library(data.table)
library(stringi)
library(jsonlite)
'%!in%' <- function(x,y)!('%in%'(x,y)) 
'%!like%' <- function(x,y)!('%like%'(x,y)) 




#################################################################################
#### Import data
#################################################################################
freyja<-read.table("freyja_lineage.tsv", fill = TRUE, sep = "\t", h=T)
lineage <- read_json("alias_key.json", simplifyVector = TRUE)
taxo<-read.table("database_lineages_table.tsv", h=T, sep = "\t")
taxo.1<-read.table("NextstrainWHO.txt", header = TRUE, sep = "\t")



#################################################################################
#### Extract special lineages/recombinant info from cov_lineages.json
#################################################################################
# Convert json to df
lineage <- cbind(names(lineage), as.data.frame(matrix(c(lineage, rep(NA, length(lineage) %% length(lineage))), length(lineage))))
lineage[,1]<-dplyr::na_if(lineage[,1], "")
names(lineage)<-c("original", "final")
lineage$n<-str_count(as.character(lineage$final), ",")+1
lineage$temp<-gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", lineage$final, perl=T)
lineage<-as.data.frame(lapply(lineage, gsub, pattern='"', replacement=""))
lineage<-as.data.frame(lapply(lineage, gsub, pattern="[*]", replacement=""))
lineage<-as.data.frame(lapply(lineage, gsub, pattern=" ", replacement=""))
lineage$final<-ifelse(lineage$temp == "", lineage$final, lineage$temp)
# Remove A and B from lineages
lineage<-lineage[which(lineage$original %!in% c("A", "B")),]

lineage<-lineage[, 1:3]
lineage.backup<-lineage


# Identify VoC within lineages and rename them
lineage.sub.voc<-unique(gsub("\\..*","",taxo.1$Pango))
lineage.sub.voc<-lineage.backup %>% filter(original %in% lineage.sub.voc)
lineage$Pango<-lineage$final




#################################################################################
#### Taxo
#################################################################################

# List all WHO lineages not in "list.variants.with.sublineages" 
taxo<-taxo %>% filter(WHO %!in% list.variants.with.sublineages) %>%
  select("WHO", "Pango.sublineages") %>% 
  drop_na() 

taxo<-as.data.frame(cbind(taxo[,1], str_split_fixed(taxo$Pango.sublineages, ',', max(str_count(as.character(taxo$Pango.sublineages), ",")+1))))
taxo<-reshape2::melt(taxo, id=c("V1"))
taxo<-taxo[, c(1, 3)]; names(taxo)<-c("WHO", "Pango")
taxo$Pango<-str_trim(taxo$Pango) # Remove spaces before Pango

# Develop lineages in "list.variants.with.sublineages" 
taxo.1<-taxo.1 %>% filter(WHO %in% list.variants.with.sublineages) %>%
  filter(WHO %in% list.variants.with.sublineages) %>%
  select("WHO", "Pango")
# Add missing lineages
taxo.1<-rbind(taxo.1, replace.lineages) #taxo.1<-rbind(taxo.1, complete.taxo)
# Create the lineages displayed on the dashboard
taxo.1$Lineage<-ifelse(taxo.1$Pango %!in% replace.lineages$Pango, paste0(taxo.1$WHO, " (", taxo.1$Pango, ")"), taxo.1$WHO)

# Add full Pango name
taxo.1$Pango.full<-taxo.1$Pango
for(i in 1:nrow(lineage)){
  taxo.1$Pango.full<-gsub(paste0("^", lineage$original[i], "[.]"), paste0(lineage$Pango[i], "."), taxo.1$Pango.full, ignore.case = FALSE)
}
taxo.1$Pango.full<-ifelse(taxo.1$Pango %in% replace.lineages$Pango, taxo.1$Pango, taxo.1$Pango.full)





#################################################################################
#### Curate Freyja's output
#################################################################################
freyja.backup<-freyja
freyja.backup->freyja
freyja$OriginalLineages<-freyja$Lineage

# Convert long sub-lineages Pango into short Pango
## Convert all pango into LONG Pango names
freyja$Pango<-freyja$OriginalLineages
for(i in which(lineage$original %!in% ignore.name.lineages$Pango.1)){
  freyja$Pango<-gsub(paste0("^", lineage$original[i], "[.]"), paste0( lineage$Pango[i], "."), freyja$Pango, ignore.case = FALSE)
}
#freyja[which(freyja$Lineage =="DK.1"),]

# Assigned lineages that are not listed in "list.variants.with.sublineages" 
freyja$Lineage<-"Other"
freyja.Lineage<-freyja$Lineage
freyja.OriginalLineages<-freyja$OriginalLineages
freyja.Pango<-freyja$Pango

for(i in 1:nrow(taxo)){
  other.row<-which(freyja.Lineage == "Other")
  #lineage.row<-which(startsWith(freyja.OriginalLineages, taxo$Pango[i]))
  lineage.row<-which(freyja.OriginalLineages == taxo$Pango[i])
  lineage.row<-lineage.row[which(lineage.row %in% other.row)]
  freyja.Lineage[lineage.row]<-taxo$WHO[i]
}

# Assigned lineages that are listed in "list.variants.with.sublineages" 
## 1
for(i in 1:nrow(taxo.1)){
  other.row<-which(freyja.Lineage == "Other")
  lineage.row<-which(startsWith(freyja.Pango, taxo.1$Pango.full[i]))
  lineage.row<-lineage.row[which(lineage.row %in% other.row)]
  freyja.Lineage[lineage.row]<-taxo.1$Lineage[i]
}
freyja$Lineage<-freyja.Lineage

## 2

for(i in 1:nrow(taxo.1)){
  other.row<-which(freyja.Lineage == "Other")
  lineage.row<-which(startsWith(freyja.Pango, taxo.1$Pango[i]))
  lineage.row<-lineage.row[which(lineage.row %in% other.row)]
  freyja.Lineage[lineage.row]<-taxo.1$Lineage[i]
}

#freyja.Pango<-freyja$Pango
# for(i in 1:nrow(taxo.1)){
#   freyja.Lineage[which(paste0("^", freyja.Pango) %like% paste0("^", taxo.1$Pango[i]) & freyja.Lineage == "Other")]<-taxo.1$Lineage[i]
# }
freyja$Lineage<-freyja.Lineage



  
#################################################################################
# Filter proportions <X%
freyja<-freyja %>% filter(as.numeric(as.character(proportion))>=min.proportion)



#################################################################################               
## Adjust proportion to reach 100% per sample
sublineages<-freyja %>%
  group_by(samples, Lineage) %>%
  summarize(proportion = round(sum(as.numeric(as.character(proportion)), na.rm=TRUE), 2), 
            .groups = 'drop')

adj<-sublineages %>%
  group_by(samples) %>%
  summarize(SUM = sum(as.numeric(as.character(proportion)), na.rm=TRUE), 
            .groups = 'drop')


sublineages<-left_join(sublineages, adj, by="samples")
sublineages$proportion<-((1/sublineages$SUM)*sublineages$proportion)
sublineages<-sublineages %>% mutate_if(is.numeric, ~replace_na(., 0))



#################################################################################               
### Export data
write.table(sublineages[,c("samples", "Lineage", "proportion")], "./freyja_plot.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(freyja, "./freyja_lineages_conversion_backup.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")



