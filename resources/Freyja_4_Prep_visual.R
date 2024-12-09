rm(list=ls(all=FALSE))
ls()
set.seed(123)

require(tidyverse)
'%!in%' <- function(x,y)!('%in%'(x,y)) 




#################################################################################

# Variant discarded if porportion below this threshold
min.proportion = 1 # 0-100

# # List sublineages
replace.lineages<-data.frame(Pango=c("B.1.1.529", "Add lineage here"),
                             WHO=c("Omicron", "Add lineage here"))

#################################################################################





#################################################################################
#### Import data
#################################################################################
freyja<-read.table("./freyja_lineage.tsv", fill = TRUE, sep = "\t", h=T)
alias.key <- read.table("./alias_key.tsv", header = TRUE, sep = "\t")
list.clade.display <- read.table("./ListVariantToDisplay.tsv", header = TRUE, sep = "\t")
list.recombinant <- read.table("./alias_key_recombinant.tsv", header = TRUE, sep = "\t")
outbreakinfo.lineages<-read.table("./outbreakinfo_db_table_melted.tsv", header = TRUE, sep = "\t")





#################################################################################
#### Prepare files: Nextstrain clades (list.clade.display)
#################################################################################

list.clade.display$Pango.full<-list.clade.display$Pango # Will be useful for the loop
alias.key.norecombinant<-alias.key %>% dplyr::filter(n == 1) # Remove recombinant (n>2) from alias.key

# Add full Pango name (eg, B.1.1.1.1.16.1), except for recombinants
for(i in 1:nrow(alias.key.norecombinant)){
  list.clade.display$Pango.full<-gsub(paste0("^", alias.key.norecombinant$original[i], "[.]"), paste0(alias.key.norecombinant$final[i], "."), list.clade.display$Pango.full, ignore.case = FALSE)
}

# Some clades don't have a full Pango name (e.g., Omicron that is B.1.1.529). Fix it using the table replace.lineages
for(i in 1:nrow(replace.lineages)){
  list.clade.display<-list.clade.display %>%
    dplyr::mutate(Pango.full = ifelse(Legend %in% replace.lineages$WHO[i], replace.lineages$Pango[i], Pango.full))
}

# Add info if clade is a recombinant or not
list.recombinant1 <- list.recombinant %>% 
  dplyr::rename(Pango_simplified = original) %>%
  dplyr::select(Pango_simplified, n)

# List the clades for which there is Pango info in list.clade.display
list.clade.display.Pango<-list.clade.display %>%
  dplyr::filter(Pango.full != "")

#################################################################################
#### Prepare files: (outbreakinfo.lineages)
#################################################################################

# In outbreakinfo.lineages, discard the clades that are listed in list.clade.display.Pango
outbreakinfo.lineages.filter <- outbreakinfo.lineages %>%
  dplyr::mutate(WHO = ifelse(is.na(WHO), "Other", WHO)) %>%
  dplyr::filter(WHO %!in% list.clade.display.Pango$WHO)





#################################################################################
#### Freyja predictions (freyja)
#################################################################################


freyja<-read.table("freyja_lineage.tsv", fill = TRUE, sep = "\t", h=T)
freyja$OriginalLineages<-freyja$Lineage #Store Lineage into OriginalLineages

#################################################################################
##### Assign WHO/clade names

## Convert all Pango into LONG Pango names (Pango.full)
freyja$Pango.full<-freyja$Lineage #initiate variable for loop
for(i in 1:nrow(alias.key.norecombinant)){
  freyja$Pango.full<-gsub(paste0("^", alias.key.norecombinant$original[i], "[.]"), paste0(alias.key.norecombinant$final[i], "."), freyja$Pango.full, ignore.case = FALSE)
}
#freyja[which(freyja$Lineage =="DK.1"),]
#freyja[which(freyja$OriginalLineages =="KP.3.1.10"),] #Should be 1.1.529.2.86.1.1.11.1.3.1.10 (Pango.full)

## Work into variables

#freyja.backup<-freyja
#freyja<-freyja.backup
freyja$Lineage<-"Other" #Was called "Lineage" in the 1st version of the script, so will keep using "Lineage" instead of "Clade"
freyja.lineage<-freyja$Lineage
freyja.OriginalLineages<-freyja$OriginalLineages
freyja.Pango.full<-freyja$Pango.full
freyja.Pango.full2<-paste0(freyja$Pango.full,".") # Add dots to make sure KP.3.1.10 is not recognized as KP.3.1.1 :-)

list.clade.display.Pango<-dplyr::arrange(list.clade.display.Pango, desc(Pango.full)) # Sort lineages (to ensure they are all well detected in the db) - Otherwise, can misassign a clade
list.clade.display.Pango.full1<-list.clade.display.Pango$Pango.full
list.clade.display.Pango.full2<-paste0(list.clade.display.Pango$Pango.full,".") # Add dots to make sure KP.3.1.10 is not recognized as KP.3.1.1 :-)



## Assign the WHO/clades for the group of variants described in list.clade.display.Pango (on 2024-11-15, it is only Omicron's variants)
for(i in 1:nrow(list.clade.display.Pango)){
  other.row<-which(freyja.lineage == "Other")
  lineage.row<-c(which(startsWith(freyja.Pango.full2, list.clade.display.Pango.full2[i])), # Assign if start Pango.full match with the clade name
                 which(freyja.Pango.full %in% list.clade.display.Pango.full1[i])) # Assign if perfect match
  lineage.row<-lineage.row[which(lineage.row %in% other.row)]
  freyja.lineage[lineage.row]<-list.clade.display.Pango$Legend[i]
}
freyja$Lineage<-freyja.lineage




## Assign the WHO/clades for the group of variants described in outbreakinfo.lineages (on 2024-11-15, it is everything but Omicron's variants)
list.clade.display<-dplyr::arrange(list.clade.display, desc(Pango.full)) # Sort lineages (to ensure they are all well detected in the db) - Otherwise, can misassign a clade
for(i in 1:nrow(outbreakinfo.lineages.filter)){
  other.row<-which(freyja.lineage == "Other")
  lineage.row<-which(freyja.OriginalLineages == outbreakinfo.lineages.filter$Pango[i])
  lineage.row<-lineage.row[which(lineage.row %in% other.row)]
  freyja.lineage[lineage.row]<-outbreakinfo.lineages.filter$WHO[i]
}
freyja$Lineage<-freyja.lineage




#################################################################################
# Filter proportions <X%
freyja<-freyja %>% dplyr::filter(as.numeric(as.character(proportion))>=min.proportion)



#################################################################################
## Final changes: Rename header + add new column

# Rename Pango.full to Pango
freyja<-freyja %>% dplyr::rename(Pango = Pango.full) 

# Add Legend position
Legend.position.df<-list.clade.display %>% 
  dplyr::filter(Nextstrain %!in% c("21I", "21J")) %>% #Delta appears 3x time, just keep the 1st one (21A) to avoid errors below
  dplyr::mutate(legend.position = ifelse(Nextstrain == "21M", 22, legend.position)) %>% # 21M-Omicron appears after 21K-Omicron-BA.1 and 21L-Omicron-BA.2, which can be confusing for people. So, AJR artificially put Omicron before BA.1 and BA.2
  dplyr::rename(Lineage = Legend) %>% 
  dplyr::select(legend.position, Lineage) %>% 
  dplyr::add_row(legend.position = 1, Lineage = "Other") %>%
  unique() 


freyja<-dplyr::left_join(freyja, Legend.position.df, by="Lineage") %>%
  dplyr::mutate(legend.position = ifelse(Lineage == "Other", 1, legend.position))



#################################################################################               
## Adjust proportion to reach 100% per sample
sublineages<-freyja %>%
  dplyr::group_by(samples, Lineage, legend.position) %>%
  dplyr::summarize(proportion = round(sum(as.numeric(as.character(proportion)), na.rm=TRUE), 2), 
            .groups = 'drop')

adj<-sublineages %>%
  dplyr::group_by(samples) %>%
  dplyr::summarize(SUM = sum(as.numeric(as.character(proportion)), na.rm=TRUE), 
            .groups = 'drop')


sublineages<-dplyr::left_join(sublineages, adj, by="samples")
sublineages$proportion<-((1/sublineages$SUM)*sublineages$proportion)
sublineages<-sublineages %>% dplyr::mutate_if(is.numeric, ~replace_na(., 0))




#################################################################################               
### Export data
write.table(sublineages[,c("samples", "Lineage", "proportion", "legend.position")], "./freyja_plot.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(freyja, "./freyja_lineages_conversion_backup.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")



