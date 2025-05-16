rm(list=ls(all=FALSE))
#ls()

library(dplyr)

# Import databases
CovariantsMutations<-read.table("covariant_mutations.tsv", h=T, sep="\t", quote = "\\", stringsAsFactors=FALSE)
OutbreakLineages<-read.table("outbreak_lineages_table.tsv", h=T, sep = "\t")
NextstrainClades<-read.table("NextstrainWHO.txt", header = TRUE, sep = "\t")


# Format CovariantsMutations database
names(CovariantsMutations)[1]<-"lineages"
CovariantsMutations<-CovariantsMutations %>% dplyr::rename("mutation.nt" = "nuc_change")
CovariantsMutations<-CovariantsMutations %>% dplyr::rename("mutation.aa" = "aa_change")
CovariantsMutations<- CovariantsMutations %>% dplyr::separate(lineages, c('Lineages', 'WHO') )
CovariantsMutations<-CovariantsMutations %>% dplyr::filter(reversion != "y" & mutation.nt != "nuc_change")


# Select VoC (according to Outbreak)
CovariantsMutations<-CovariantsMutations %>% dplyr::filter(WHO %in% OutbreakLineages[which(OutbreakLineages$Status == "Variant of Concern"), c("WHO")])


# Prepare Nextstrain db
NextstrainClades<-NextstrainClades %>% dplyr::select(-WHO)
NextstrainClades <- NextstrainClades %>% dplyr::rename(Lineages = Nextstrain)


# Merge CovariantsMutations and Nextstrain db
CovariantsMutations<-dplyr::left_join(CovariantsMutations, NextstrainClades, by=("Lineages"))
CovariantsMutations$Lineages<-paste0(CovariantsMutations$Lineages, " (", CovariantsMutations$WHO, ")")



# Export
write.table(CovariantsMutations %>% dplyr::select(Lineages, mutation.nt, mutation.aa), "voc_taxo_mut_curated.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


