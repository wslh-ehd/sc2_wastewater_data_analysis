rm(list=ls(all=FALSE))
ls()
set.seed(123)

library(tidyverse)
library(openxlsx)
library(lubridate)
library(data.table)

'%!in%' <- function(x,y)!('%in%'(x,y)) 
check.integer <- function(N){!grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))}


#################################################################################               
#################################################################################               
# Config
interval.to.investigate = 2









#################################################################################               
## Import data
data<-read.delim("CallVariantALLCompiled.tsv", h=T, sep = "\t")
depth<-read.delim("CallDepthCompiled.tsv", h=F, sep = "\t")

samplesinfo<-read.xlsx("ListSamples.xlsx", sheet="Samples", startRow = 3)
runsinfo<-read.xlsx("ListSamples.xlsx", sheet="Runs", startRow = 10)

mutations.variants<-read.delim("database_NR_outbreak_lineages.tsv", h=T, sep = "\t")

mut.covariants<-read.delim("voc_taxo_mut_curated.tsv", h=T, sep="\t", quote = "\\", stringsAsFactors=FALSE)
mutations.voc.force<-read.delim("database_gisaid.voc.force.tsv", h=T, sep = "\t")

gene<-read.delim(url("https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/data/PositionGenes.txt"), h=T, sep = "\t")





#################################################################
# SARS-CoV-2 genome description            
#################################################################
genes_names <- c(
  '1' = "ORF1a",
  '2' = "ORF1b", 
  '3' = "S", 
  '4' = "ORF3a", 
  '5' = "E",
  '6' = "M",
  '7' = "ORF6",
  '8' = "ORF7a",
  '9' = "ORF7b",
  '10' = "ORF8",
  '11' = "N",
  '12' = "ORF10",
  '13' = "intergenic")

for(i in 1:nrow(gene)){
  if(i==1){
    genome<-0
    genome<-as.data.frame(seq(gene[i, 2], gene[i, 3], by=1)); names(genome)<-"POS"
    genome$gene<-gene[i, 1]
    genome$geneplot<-gene[i, 4]
    genome$aa.pos<-floor((genome$POS-min(genome$POS)+3)/3)
  }else{
    temp<-0
    temp<-as.data.frame(seq(gene[i, 2], gene[i, 3], by=1)); names(temp)<-"POS"
    temp$gene<-gene[i, 1]
    temp$geneplot<-gene[i, 4]
    temp$aa.pos<-floor((temp$POS-min(temp$POS)+3)/3)
    genome<-rbind(genome, temp)
  }
}
genome<-as.data.frame(genome)
genome$POS<-as.numeric(as.character(genome$POS))
genome$codon.pos<-rep(seq(1,3,1), times = nrow(genome)/3)
genome$R<-paste0(genome$gene,".", genome$aa.pos, ".", genome$codon.pos)






#################################################################
# Preparation: data 1/3               
#################################################################

# Fix formatting
data$ALT_FREQ<-as.numeric(as.character(data$ALT_FREQ))
data$TOTAL_DP<-as.numeric(as.character(data$TOTAL_DP))
data$POS<-as.numeric(as.character(data$POS))

# Rename columns
names(data)[1:2]<-c("run", "samples")

# Remove rows that contains unwanted info
data<-data[which(is.na(data$POS)==FALSE), ]

# Select/create the columns of interest
data<-data %>% dplyr::select("run", "samples", "POS", "ALT_FREQ", "TOTAL_DP", "REF", "ALT", "REF_AA")
data$run_sample<-paste0(data$run, "_", data$samples)




#################################################################
# Preparation: data 2/3 - Add deletions/insertions            
#################################################################
data$type<-ifelse(data$ALT %like% '[-]', "deletion", ifelse(data$ALT %like% '[+]', "insertion", "substitution"))
data$POS<-as.numeric(as.character(data$POS)); genome$POS<-as.numeric(as.character(genome$POS))
data<-as.data.frame(setDT(genome)[data, on="POS"])

# Check if insertion/deletion are a multiple of 3, if so, determine the mutation
data$No_aa_indel<-ifelse(data$type != "substitution", nchar(data$ALT)-1, NA)/3
data$No_aa_indel<-ifelse(round(data$No_aa_indel) == data$No_aa_indel, data$No_aa_indel, NA)

data$deletion<-ifelse(data$No_aa_indel==1, 
                      paste0(data$gene, ":del", data$aa.pos+1), 
                      paste0(data$gene, ":del", data$aa.pos+1, "/", data$aa.pos+data$No_aa_indel))



#################################################################
#  Preparation: data 3/3 - Create new variables         
#################################################################
data$ALT_FREQ<-round(data$ALT_FREQ*100, digits = 2)
data$TOTAL_DP<-ifelse(data$TOTAL_DP>100,100,data$TOTAL_DP)
data$run_sample_POS<-paste0(data$run, "_", data$samples, "_",data$POS)




#################################################################
# Preparation: runsinfo and samplesinfo               
#################################################################
runsinfo<-runsinfo[, c("Run", "QC_Run")]
samplesinfo <- left_join(samplesinfo, runsinfo, by=c("Run"))

samplesinfo <- samplesinfo %>% rename(samples = FilesNames)
samplesinfo <- samplesinfo %>% rename(sites = Location)

samplesinfo <-samplesinfo %>% dplyr::filter(!is.na(samples))
samplesinfo$Date<-as.Date(as.numeric(as.character(samplesinfo$Date)), origin = "1899-12-30")
samplesinfo$Date<-as.Date(samplesinfo$Date, "%m%d%y")
samplesinfo$Samples<-paste0(samplesinfo$sites, "_", samplesinfo$Date)
samplesinfo$Month<-lubridate::interval(samplesinfo$Date, lubridate::today()) %/% months(1)
samplesinfo$run_sample<-paste0(samplesinfo$Run, "_", samplesinfo$samples)

samplesinfo<-samplesinfo %>% dplyr::select("Samples", "sites", "Date", "Sample.type", "QC_Sample",  "QC_Run", "run_sample", "Month")



#################################################################
# Merge data & samplesinfo               
#################################################################
data <- left_join(data, samplesinfo, by=c("run_sample"))



#################################################################
# Select samples
#################################################################
data <-data %>% dplyr::filter(Sample.type == "Sample" & QC_Run == "qc_passed" & Month < interval.to.investigate)
data <-data %>% dplyr::filter(QC_Sample %in% c("hold", "qc_passed"))




#################################################################
# Preparation: Depth               
#################################################################

# Rename columns
names(depth)[1:2]<-c("run", "samples")
names(depth)[4]<-"POS"
names(depth)[6]<-"TOTAL_DP"

# Fix formatting
depth$TOTAL_DP<-as.numeric(as.character(depth$TOTAL_DP))
depth$POS<-as.numeric(as.character(depth$POS))

# Select/create the columns of interest
depth$run_sample<-paste0(depth$run, "_", depth$samples)
depth$run_sample_POS<-paste0(depth$run, "_", depth$samples, "_",depth$POS)
depth$TOTAL_DP<-ifelse(depth$TOTAL_DP>100,100,depth$TOTAL_DP)
depth$ALT_FREQ<-0

depth<-depth %>% 
  dplyr::select("TOTAL_DP", "POS", "ALT_FREQ", "run_sample", "run_sample_POS") %>%
  dplyr::filter(run_sample %in% as.character(unique(data$run_sample)))
  




#################################################################
# Add mutation info to data
#################################################################

# Replace NA describing intergenic positions by "inter"
data$temp<-paste0(data$gene, ":", data$REF_AA, data$aa.pos)
data$gene<-ifelse(is.na(data$gene), "inter", data$temp)
data<-data %>% dplyr::select(-temp)


# Create mutation.nt variable
data$mutation.nt<-ifelse(data$type=="substitution",
                             paste0(data$REF, data$POS, data$ALT),
                             paste0(data$deletion))
data$mutation.nt<-ifelse(data$mutation.nt=="NA",
                             paste0(data$REF, data$POS, "(", data$ALT, ")"),
                             data$mutation.nt)
data$mutation.nt<-ifelse(data$deletion %like% "NA:del",
                         paste0(data$gene, "_", data$POS, ":DEL", data$ALT),
                         data$mutation.nt)


data$mutation.nt<-toupper(data$mutation.nt)





#save.image(file="data.RData")



###############################################################################################
####### Visual preparation ####################################################################
###############################################################################################


## Visual 1: VoC overview (based on covariant.org)

# Preparation: NEXTSTRAIN
data.voc.nextstrain<- left_join(mut.covariants, data, by=c("mutation.nt"))
data.voc.nextstrain<-data.voc.nextstrain %>% dplyr::filter(!is.na(POS))
data.voc.nextstrain<-data.voc.nextstrain[, c("mutation.nt", "Lineages", "mutation.aa", "sites", "Date", "POS", "gene", "geneplot", "ALT_FREQ", "TOTAL_DP", "run", "type", "Samples", "run_sample", "run_sample_POS", "QC_Run", "QC_Sample", "Sample.type")]
order1_xaxis<-data.voc.nextstrain %>% dplyr::arrange(Date)
# Add missing depth: 1. Identify the mutations in data.voc.nextstrain that don't have depth info
data.voc<-reshape2::dcast(data.voc.nextstrain, run_sample~mutation.nt,
                               value.var="mutation.nt",fill=0,fun.aggregate=length)
data.voc<-reshape2::melt(data.voc, id.vars = "run_sample")
data.voc<-data.voc %>% dplyr::filter(run_sample != "NA" & value==0)
data.voc$run_sample_POS<-paste0(data.voc$run_sample, "_", parse_number(as.character(data.voc$variable)))
# Add missing depth: 2. Extract missing depths
depth.voc.nextstrain<-depth %>% dplyr::filter(run_sample_POS %in% unique(data.voc$run_sample_POS))
# Add missing depth: 3. Add samples info to depth
depth.voc.nextstrain<-as.data.frame(setDT(unique(data.voc.nextstrain[, c("sites", "Date",  "run", "Samples", "run_sample", "QC_Run", "QC_Sample", "Sample.type")]))[depth.voc.nextstrain, on="run_sample"])
depth.voc.nextstrain<-as.data.frame(setDT(unique(data.voc.nextstrain[, c("mutation.nt", "Lineages", "POS", "mutation.aa", "gene", "geneplot",  "type")]))[depth.voc.nextstrain, on="POS",allow.cartesian=TRUE])
depth.voc.nextstrain<-depth.voc.nextstrain %>% dplyr::filter(paste0(run_sample_POS) %!in% paste0(data.voc.nextstrain$run_sample_POS))
# Add missing depth: 4. Merge data and depth datasets
data.voc.nextstrain<-full_join(data.voc.nextstrain, depth.voc.nextstrain, by=intersect(colnames(data.voc.nextstrain), colnames(depth.voc.nextstrain)))
# Add tooltips info
data.voc.nextstrain$display<-paste0(data.voc.nextstrain$gene, "_", data.voc.nextstrain$mutation.nt)
# Clean data
rm(data.voc); rm(depth.voc.nextstrain)






## Visual 2: Overview of variants of interest (manual entry [forced])

# Preparation: CDC https://covid.cdc.gov/covid-data-tracker/#variant-proportions
data.voc.force<- left_join(mutations.voc.force %>% dplyr::select(-POS, -gene, -aa.pos, -type), data, by=c("mutation.nt"))
data.voc.force<-data.voc.force %>% dplyr::filter(!is.na(POS))
#data.voc.force <- data.voc.force %>% rename(Lineages = i.lineage)
data.voc.force<-data.voc.force[, c("mutation.nt", "lineage", "mutation.aa", "sites", "Date", "POS", "gene", "geneplot", "ALT_FREQ", "TOTAL_DP", "run", "type", "Samples", "run_sample", "run_sample_POS", "QC_Run", "QC_Sample", "Sample.type")]
order2_xaxis<-data.voc.force %>% dplyr::arrange(Date)
# Add missing depth: 1. Identify the mutations in data.voc.force that don't have depth info
data2.voc<-reshape2::dcast(data.voc.force, run_sample~mutation.nt,
                               value.var="mutation.nt",fill=0,fun.aggregate=length)
data2.voc<-reshape2::melt(data2.voc, id.vars = "run_sample")
data2.voc<-data2.voc %>% dplyr::filter(run_sample != "NA" & value==0)
data2.voc$run_sample_POS<-paste0(data2.voc$run_sample, "_", parse_number(as.character(data2.voc$variable)))
# Add missing depth: 2. Extract missing depths
depth.voc.force<-depth %>% dplyr::filter(run_sample_POS %in% unique(data2.voc$run_sample_POS))
# Add missing depth: 3. Add samples info to depth
depth.voc.force<-as.data.frame(setDT(unique(data.voc.force[, c("sites", "Date",  "run", "Samples", "run_sample", "QC_Run", "QC_Sample", "Sample.type")]))[depth.voc.force, on="run_sample"])
depth.voc.force<-as.data.frame(setDT(unique(data.voc.force[, c("mutation.nt", "lineage", "POS", "mutation.aa",  "gene", "geneplot",  "type")]))[depth.voc.force, on="POS",allow.cartesian=TRUE])
depth.voc.force<-depth.voc.force %>% dplyr::filter(paste0(run_sample_POS) %!in% paste0(data.voc.force$run_sample_POS))
# Add missing depth: 4. Merge data and depth datasets
data.voc.force<-full_join(data.voc.force, depth.voc.force, by=intersect(colnames(data.voc.force), colnames(depth.voc.force)))
# Add tooltips info
data.voc.force$display<-paste0(data.voc.force$gene, "_", data.voc.force$mutation.nt)
# Clean data
rm(data2.voc); rm(depth.voc.force)





## Visual 3: VoC-associated mutations (based on Outbreak database)

# Merge mutations & data 
data.snp.voc <- right_join(data, mutations.variants %>% dplyr::filter(status == "Variant of Concern") %>% dplyr::select(mutation.nt, mutation.aa, lineage, SNP.lineage.count), by=c("mutation.nt"))
data.snp.voc<-data.snp.voc %>% dplyr::filter(!is.na(POS))
data.snp.voc<-data.snp.voc[, c("mutation.nt", "lineage", "SNP.lineage.count", "mutation.aa", "sites", "Date", "POS", "gene", "geneplot", "ALT_FREQ", "TOTAL_DP", "run", "type", "Samples", "run_sample", "run_sample_POS", "QC_Run", "QC_Sample", "Sample.type")]
order3_xaxis<-data.snp.voc %>% dplyr::arrange(Date)
# Add missing depth: 1. Identify the mutations in data.snp.voc that don't have depth info
data.snp.voc.2<-reshape2::dcast(data.snp.voc, run_sample~mutation.nt,
                          value.var="mutation.nt",fill=0,fun.aggregate=length)
data.snp.voc.2<-reshape2::melt(data.snp.voc.2, id.vars = "run_sample")
data.snp.voc.2<-data.snp.voc.2 %>% dplyr::filter(run_sample != "NA" & value==0)
data.snp.voc.2$run_sample_POS<-paste0(data.snp.voc.2$run_sample, "_", parse_number(as.character(data.snp.voc.2$variable)))
# Add missing depth: 2. Extract missing depths
depth.snp.voc<-depth %>% dplyr::filter(run_sample_POS %in% unique(data.snp.voc.2$run_sample_POS))
# Add missing depth: 3. Add samples info to depth
depth.snp.voc<-as.data.frame(setDT(unique(data.snp.voc[, c("sites", "Date",  "run", "Samples", "run_sample", "QC_Run", "QC_Sample", "Sample.type")]))[depth.snp.voc, on="run_sample"])
depth.snp.voc<-as.data.frame(setDT(unique(data.snp.voc[, c("mutation.nt", "lineage", "SNP.lineage.count", "POS", "mutation.aa", "gene", "geneplot",  "type")]))[depth.snp.voc, on="POS",allow.cartesian=TRUE])
depth.snp.voc<-depth.snp.voc %>% dplyr::filter(paste0(run_sample_POS) %!in% paste0(data.snp.voc$run_sample_POS))
# Add missing depth: 4. Merge data and depth datasets
data.snp.voc<-full_join(data.snp.voc, depth.snp.voc, by=intersect(colnames(data.snp.voc), colnames(depth.snp.voc)))
# Add tooltips info
data.snp.voc$display<-paste0(data.snp.voc$gene, "_", data.snp.voc$mutation.nt)
# Clean data
rm(data.snp.voc.2); rm(depth.snp.voc)






## Visual 4: NOT VoC-associated mutations (based on Outbreak database)

# Merge mutations & data 
data.snp.notvoc<-data %>% dplyr::filter(mutation.nt %!in% unique(data.snp.voc$mutation.nt) & !is.na(POS))
data.snp.notvoc <- full_join(data.snp.notvoc, mutations.variants %>% dplyr::filter(status != "Variant of Concern") %>% dplyr::select(mutation.nt, mutation.aa, lineage), by=c("mutation.nt"))
data.snp.notvoc<-data.snp.notvoc[, c("mutation.nt", "lineage", "mutation.aa", "sites", "Date", "POS", "gene", "geneplot", "ALT_FREQ", "TOTAL_DP", "run", "type", "Samples", "run_sample", "run_sample_POS", "QC_Run", "QC_Sample", "Sample.type")]
order4_xaxis<-data.snp.notvoc %>% dplyr::arrange(Date)
# Add missing depth: 1. Identify the mutations in data.snp.notvoc that don't have depth info
data.snp.notvoc.2<-reshape2::dcast(data.snp.notvoc, run_sample~mutation.nt,
                            value.var="mutation.nt",fill=0,fun.aggregate=length)
data.snp.notvoc.2<-reshape2::melt(data.snp.notvoc.2, id.vars = "run_sample")
data.snp.notvoc.2<-data.snp.notvoc.2 %>% dplyr::filter(run_sample != "NA" & value==0)
data.snp.notvoc.2$run_sample_POS<-paste0(data.snp.notvoc.2$run_sample, "_", parse_number(as.character(data.snp.notvoc.2$variable)))
# Add missing depth: 2. Extract missing depths
depth.snp.notvoc<-depth %>% dplyr::filter(run_sample_POS %in% unique(data.snp.notvoc.2$run_sample_POS))
# Add missing depth: 3. Add samples info to depth
depth.snp.notvoc<-as.data.frame(setDT(unique(data.snp.notvoc[, c("sites", "Date",  "run", "Samples", "run_sample", "QC_Run", "QC_Sample", "Sample.type")]))[depth.snp.notvoc, on="run_sample"])
depth.snp.notvoc<-as.data.frame(setDT(unique(data.snp.notvoc[, c("mutation.nt", "lineage", "POS", "mutation.aa", "gene", "geneplot",  "type")]))[depth.snp.notvoc, on="POS",allow.cartesian=TRUE])
depth.snp.notvoc<-depth.snp.notvoc %>% dplyr::filter(paste0(run_sample_POS) %!in% paste0(data.snp.notvoc$run_sample_POS))
# Add missing depth: 4. Merge data and depth datasets
data.snp.notvoc<-full_join(data.snp.notvoc, depth.snp.notvoc, by=intersect(colnames(data.snp.notvoc), colnames(depth.snp.notvoc)))
# Add tooltips info
data.snp.notvoc$display<-paste0(data.snp.notvoc$gene, "_", data.snp.notvoc$mutation.nt)
# Clean data
rm(data.snp.notvoc.2); rm(depth.snp.notvoc)


###############################################################################################
####### Export data for plotting ##############################################################
###############################################################################################

# Prepare data for plot
data.voc.nextstrain<- data.voc.nextstrain %>% dplyr::select(Lineages, mutation.nt, mutation.aa, sites, Date, ALT_FREQ, TOTAL_DP, run, Samples, display, QC_Run, QC_Sample, Sample.type)
data.voc.force<-data.voc.force %>% dplyr::select(lineage, mutation.nt, mutation.aa, sites, Date, ALT_FREQ, TOTAL_DP, run, Samples, display, QC_Run, QC_Sample, Sample.type)
data.snp.voc<- data.snp.voc %>% dplyr::select(mutation.nt, mutation.aa, lineage, SNP.lineage.count, sites, Date, ALT_FREQ, TOTAL_DP, run, Samples, display, POS, QC_Run, QC_Sample, Sample.type)
data.snp.notvoc<- data.snp.notvoc %>% dplyr::select(mutation.nt, mutation.aa, lineage, sites, Date, ALT_FREQ, TOTAL_DP, run, Samples, display, POS, QC_Run, QC_Sample, Sample.type)

order1_xaxis<-order1_xaxis$Samples
order2_xaxis<-order2_xaxis$Samples
order3_xaxis<-order3_xaxis$Samples
order4_xaxis<-order4_xaxis$Samples

save(data.voc.nextstrain, data.voc.force, data.snp.voc, data.snp.notvoc, order1_xaxis, order2_xaxis, order3_xaxis, order4_xaxis, mutations.variants, mut.covariants, mutations.voc.force,  file = "plot.RData")


