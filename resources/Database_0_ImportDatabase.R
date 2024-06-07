#rm(list=ls(all=FALSE))
#ls()
set.seed(123)


library(outbreakinfo) #devtools::install_github("outbreak-info/R-outbreak-info")
library(jsonlite)
library(dplyr)
library(data.table)
library(tidyr)



######################################################################
######################################################################
#### Functions #######################################################
'%!like%' <- function(x,y)!('%like%'(x,y)) 
'%!in%' <- function(x,y)!('%in%'(x,y)) 


extract_mutations_info <-function(x, output) {
  
  # Extract lineage
  count_loop<<-count_loop+1
  lineage = x[1]
  print(paste0(lineage, " (", count_loop, "/", number.lineage, ")"))
  
  # Extract nt mutations info
  #data.nt = read.table(paste0(query.nt, lineage), h=T)
  data.nt <- read_json(paste0(query.nt, lineage), simplifyVector = TRUE)
  if(is.data.frame(data.nt$data)){
    data.nt<-as.data.frame(data.nt$data)
    data.nt <- rename(data.nt, mutation.nt = mutation)
    data.nt$POS<-as.numeric(as.character(gsub(".([0-9]+).*$", "\\1", data.nt$mutation.nt)))
    data.nt<-as.data.frame(setDT(genome)[data.nt, on="POS"])
    data.nt<-data.nt %>% mutate(gene = replace_na(gene, "intergenic"))
    
    # Extract aa mutations info
    #data.aa<-read.table(paste0(query.aa, lineage), h=T)
    data.aa <- read_json(paste0(query.aa, lineage), simplifyVector = TRUE)
    if(is.data.frame(data.aa$data)){
      data.aa <- as.data.frame(data.aa$data)
      data.aa <- rename(data.aa, mutation.aa = mutation)
      data.aa$genePOS<-paste0(sub('\\:.*', '', data.aa$mutation.aa), ":", 
                              as.numeric(as.character(gsub(".([0-9]+).*$", "\\1", sub('.*\\:', '', data.aa$mutation.aa)))))
      
      # Combine aa and nt mutations info 
      data.mut<-as.data.frame(setDT(data.aa[, c("mutation.aa", "genePOS")])[data.nt, on="genePOS"])
      data.mut$lineage<-lineage
      
      } else {
        data.nt$mutation.aa<-NA
        data.mut<-data.nt
        data.mut$lineage<-lineage
      }
    
    # Return info
    return(data.mut)
    
    } else {
    
    data.nt <- read_json(paste0(query.nt.notUSA, lineage), simplifyVector = TRUE)
    if(is.data.frame(data.nt$data)){
      data.nt<-as.data.frame(data.nt$data)
      data.nt <- rename(data.nt, mutation.nt = mutation)
      data.nt$POS<-as.numeric(as.character(gsub(".([0-9]+).*$", "\\1", data.nt$mutation.nt)))
      data.nt<-as.data.frame(setDT(genome)[data.nt, on="POS"])
      data.nt<-data.nt %>% mutate(gene = replace_na(gene, "intergenic"))
      
      # Extract aa mutations info
      data.aa <- read_json(paste0(query.aa.notUSA, lineage), simplifyVector = TRUE)
      if(is.data.frame(data.aa$data)){
        data.aa <- as.data.frame(data.aa$data)
        data.aa <- rename(data.aa, mutation.aa = mutation)
        data.aa$genePOS<-paste0(sub('\\:.*', '', data.aa$mutation.aa), ":", 
                                as.numeric(as.character(gsub(".([0-9]+).*$", "\\1", sub('.*\\:', '', data.aa$mutation.aa)))))
        
        # Combine aa and nt mutations info 
        data.mut<-as.data.frame(setDT(data.aa[, c("mutation.aa", "genePOS")])[data.nt, on="genePOS"])
        data.mut$lineage<-lineage
        
        } else {
          
        data.nt$mutation.aa<-NA
        data.mut<-data.nt
        data.mut$lineage<-lineage
        }
      
      # Return info
      return(data.mut)
      
      }
  }
 
}
  




extract_lineages_info <-function(x, output) {
  variant.lineages<-t(data.frame(do.call('rbind', strsplit(as.character(x[2]),' OR ',fixed=TRUE))))
  variant.lineages<-as.data.frame(variant.lineages %>% unique())
  variant.lineages$variant<-x[1]
  variant.lineages$status<-x[3]
  return(variant.lineages)
}




######################################################################
######################################################################
#### Config ##########################################################

# Authenticate yourself
# outbreakinfo::authenticateUser()


# Set frequency of the mutation in the population
mut.freq=0.90

# Minimum number of sequences observed in a lineage to consider this lineage as representative
min.num.seq.in.lineage = 5



######################################################################
######################################################################
#### Import (from physical files) ####################################


# Database listing specific variants (Nextstrain names) we want to look at - Text file with variant names separated by commas
VoC.variant.force=toupper(as.character(read.table("Database_Outbreak_VoC_to_explore.txt", h=F, sep = ",")))

# SARS-CoV-2 gene location
gene<-read.table(url("https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/data/PositionGenes.txt"), h=T, sep = "\t")









#################################################################
# SARS-CoV-2 genome description            
#################################################################

for(i in 1:nrow(gene)){
  if(i==1){
    genome<-0
    genome<-as.data.frame(seq(gene[i, 2], gene[i, 3], by=1)); names(genome)<-"POS"
    genome$gene<-gene[i, 1]
    genome$aa.pos<-floor((genome$POS-min(genome$POS)+3)/3)
  }else{
    temp<-0
    temp<-as.data.frame(seq(gene[i, 2], gene[i, 3], by=1)); names(temp)<-"POS"
    temp$gene<-gene[i, 1]
    temp$aa.pos<-floor((temp$POS-min(temp$POS)+3)/3)
    genome<-rbind(genome, temp)
  }
}
genome<-as.data.frame(genome)
genome$POS<-as.numeric(as.character(genome$POS))
genome$genePOS<-paste0(genome$gene, ":", genome$aa.pos)





##########################################################################################################
#### Extract Pango lineages - using Outbreak.info R package ##############################################
##########################################################################################################

# Query the API
response <- fromJSON("https://lapis.cov-spectrum.org/open/v2/sample/aggregated?fields=pangoLineage")

# Check for errors
errors <- response$errors
if (length(errors) > 0) {
  stop("Errors")
}
# Check for deprecation
deprecationDate <- response$info$deprecationDate
if (!is.null(deprecationDate)) {
  warning(paste0("This version of the API will be deprecated on ", deprecationDate,
                 ". Message: ", response$info$deprecationInfo))
}
# The data is good to be used!
data.lineages <- response$data

# Selection (no NA and only look for lineages with a minimum of representative genomes)
data.lineages <- data.lineages %>% na.omit()
data.lineages <- data.lineages %>% dplyr::filter(count >= min.num.seq.in.lineage)
data.lineages <- data.lineages %>% dplyr::select(pangoLineage, count)



##########################################################################################################
#### Extract nt/aa mutations for each lineage - using covSPECTRUM.com website ############################
##########################################################################################################

## Pangolin
query.nt<-paste0("https://lapis.cov-spectrum.org/open/v2/sample/nucleotideMutations?countryExposure=USA&minProportion=", mut.freq, "&pangoLineage=")
query.aa<-paste0("https://lapis.cov-spectrum.org/open/v2/sample/aminoAcidMutations?countryExposure=USA&minProportion=", mut.freq, "&pangoLineage=")
query.nt.notUSA<-paste0("https://lapis.cov-spectrum.org/open/v2/sample/nucleotideMutations?minProportion=", mut.freq, "&pangoLineage=")
query.aa.notUSA<-paste0("https://lapis.cov-spectrum.org/open/v2/sample/aminoAcidMutations?minProportion=", mut.freq, "&pangoLineage=")


number.lineage<-nrow(data.lineages); count_loop = 0

data.mutations<-apply(data.lineages, 1, extract_mutations_info)
data.mutations.backup<-data.mutations
#load("data.mutations.backup.RData"); data.mutations<-data.mutations.backup
data.mutations<-data.mutations[lapply(data.mutations,length)>0]
data.mutations <- do.call("rbind", data.mutations)

# Add information
data.mutations$type<-ifelse(data.mutations$mutation.nt %like% "[-]", "deletion", 
                      ifelse(data.mutations$mutation.nt %like% "[+]", "deletion", "substitution"))

# For which lineages no mutation was found?
lineages.not.found<-levels(as.factor(data.lineages$pangoLineage))[which(levels(as.factor(data.lineages$pangoLineage)) %!in% levels(as.factor(data.mutations$lineage)))]
write.table(lineages.not.found, "database_MutationsNotFoundLineages.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




## Nextstrain
query.nt<-paste0("https://lapis.cov-spectrum.org/open/v2/sample/nucleotideMutations?&minProportion=", mut.freq, "&nextcladePangoLineage=")
query.aa<-paste0("https://lapis.cov-spectrum.org/open/v2/sample/aminoAcidMutations?&minProportion=", mut.freq, "&nextcladePangoLineage=")

data.lineage.ns<-data.frame(pangoLineage = VoC.variant.force,
                            count = 1)
number.lineage<-nrow(data.lineage.ns); count_loop = 0

data.mutations.ns<-apply(data.lineage.ns, 1, extract_mutations_info)
data.mutations.ns.backup<-data.mutations.ns
#load("data.mutations.backup.RData"); data.mutations<-data.mutations.backup
data.mutations.ns<-data.mutations.ns[lapply(data.mutations.ns,length)>0]
data.mutations.ns <- do.call("rbind", data.mutations.ns)

# Add information
data.mutations.ns$type<-ifelse(data.mutations.ns$mutation.nt %like% "[-]", "deletion", 
                            ifelse(data.mutations.ns$mutation.nt %like% "[+]", "deletion", "substitution"))




##########################################################################################################
#### Extract info about lineages #########################################################################
##########################################################################################################

# Extract variants
curated = getCuratedLineages()


### Retrieve Outbreak.info data (all variants)
all.variant = curated %>% 
  select(who_name, char_muts_parent_query, variantType)
all.variant<-all.variant %>% filter(!grepl("[*]", who_name))


# Overwrite outbreak.info data (to keep classify Omicron as "Variant of Concern")
row.overwrite<-which(all.variant$who_name %in% "Omicron")
all.variant[row.overwrite, c("variantType")]<-"Variant of Concern"


# List all lineages
data.lineages<-apply(all.variant, 1, extract_lineages_info)
data.lineages <- do.call("rbind", data.lineages)
names(data.lineages)[1]<-"lineage"


# Replace "NA" WHO naming by "Other
data.lineages$variant <- data.lineages$variant%>% replace_na('Other')


### Reformat all.variant
all.variant$char_muts_parent_query<-gsub(" OR ", ", ", all.variant$char_muts_parent_query)
names(all.variant)<-c("WHO", "Pango sublineages", "Status")




##########################################################################################################
#### Combine mutation and lineage databases ##############################################################
##########################################################################################################

data<-as.data.frame(setDT(data.mutations)[data.lineages, on="lineage"])
lineages.not.listed<-data.mutations[which(data.mutations$lineage %!in% unique(data$lineage)),]
lineages.not.listed$variant<-"NotClassified"
lineages.not.listed$status<-"NotClassified"
data<-rbind(data, lineages.not.listed)


data.agg<- data %>%
  filter(status != "Variant of Concern") %>% 
  group_by(mutation.nt) %>% 
  summarise(POS = mean(as.numeric(as.character(POS))),
            lineage = paste(lineage, collapse=","), 
            mutation.aa = paste(unique(mutation.aa), collapse=","),
            type = paste(unique(type), collapse=","),
            variant = paste(unique(variant), collapse=","),
            status = paste(unique(status), collapse=","),
            SNP.lineage.count = NA)

data.VoC.agg<- data %>% 
  filter(status == "Variant of Concern") %>% 
  group_by(mutation.nt) %>% 
  summarise(POS = mean(as.numeric(as.character(POS))),
            lineage = paste(lineage, collapse=","), 
            mutation.aa = paste(unique(mutation.aa), collapse=","),
            type = paste(unique(type), collapse=","),
            variant = paste(unique(variant), collapse=","),
            status = paste(unique(status), collapse=","),
            SNP.lineage.count = n())


data<-rbind(data.agg, data.VoC.agg)

data.notredundant<-rbind(data.agg %>% filter(data.agg$mutation.nt %!in%  data.VoC.agg$mutation.nt), data.VoC.agg)





##########################################################################################################
#### Export ##############################################################################################
##########################################################################################################
write.table(all.variant, "outbreak_lineages_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(data.lineages, "database_details_lineages.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.mutations, "database_details_mutations.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(data.mutations.ns, "database_gisaid.voc.force.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data, "database_outbreak_lineages.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.notredundant, "database_NR_outbreak_lineages.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



