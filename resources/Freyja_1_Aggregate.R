rm(list=ls(all=FALSE))
ls()
set.seed(123)


library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(tidyverse)

'%!in%' <- function(x,y)!('%in%'(x,y))
'%!like%' <- function(x,y)!('%like%'(x,y))


################################################################################
#### VoC
################################################################################


#### Import data #### 
file.list <- list.files("./bootstraps/", pattern='*_summarized.csv')
for (i in 1:length(file.list)){
  temp<-as.data.frame(t(read.table(paste0("./bootstraps/", file.list[i]), sep = ",", h=F)))
  temp$sample<-file.list[i]
  if(i==1){data<-temp[-1,]}else{data<-rbind(data, temp[-1,])}
} 

#### Fix factors ####
data.VoC<-data
names(data.VoC)<-c("Lineage", "q2.5", "q5", "q25", "q50", "q75", "q95", "q97.5", "samples")
data.VoC$samples<-gsub("_summarized.csv", "", data.VoC$samples)


####  Rename A and B #### 
data.VoC$Lineage<-ifelse(data.VoC$Lineage %in% c("A", "B"), "Other", data.VoC$Lineage)


## Adjust abundance to reach 100% per sample
VoC<-data.VoC %>%
  group_by(samples, Lineage, q50) %>%
  summarize(proportion = round(sum(as.numeric(as.character(q50)), na.rm=TRUE)*100,2), 
            .groups = 'drop')

adj.VoC<-VoC %>%
  group_by(samples) %>%
  summarize(SUM = sum(as.numeric(as.character(proportion)), na.rm=TRUE), 
            .groups = 'drop')

VoC<-left_join(VoC, adj.VoC, by="samples")
VoC$proportion<-((1/VoC$SUM)*VoC$proportion)*100


# Remove the assignments with a proportion of 0 after rounding
VoC <- VoC %>% filter(proportion > 0)


#################################################################################               
### Export data
#write.table(data.VoC, "./freyja_bootstraps_VoC.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(VoC[,c("samples", "Lineage", "proportion")], "./freyja_VoC.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")




























################################################################################
#### Lineages
################################################################################

#### Import data #### 
file.list <- list.files("./bootstraps/", pattern='*_lineages.csv')
for (i in 1:length(file.list)){
  temp<-as.data.frame(t(read.table(paste0("./bootstraps/", file.list[i]), sep = ",", h=F)))
  temp$sample<-file.list[i]
  if(i==1){data<-temp[-1,]}else{data<-rbind(data, temp[-1,])}
} 

#### Fix factors ####
data.lineages<-data
names(data.lineages)<-c("Lineage", "q2.5", "q5", "q25", "q50", "q75", "q95", "q97.5", "samples")
data.lineages$samples<-gsub("_lineages.csv", "", data.lineages$samples)


####  Rename A and B #### 
data.lineages$Lineage<-ifelse(data.lineages$Lineage %in% c("A", "B"), "Other", data.lineages$Lineage)


## Adjust abundance to reach 100% per sample
lineage<-data.lineages %>%
  group_by(samples, Lineage, q50) %>%
  summarize(proportion = round(sum(as.numeric(as.character(q50)), na.rm=TRUE)*100,2), 
            .groups = 'drop')

adj.lineage<-lineage %>%
  group_by(samples) %>%
  summarize(SUM = sum(as.numeric(as.character(proportion)), na.rm=TRUE), 
            .groups = 'drop')

lineage<-left_join(lineage, adj.lineage, by="samples")
lineage$proportion<-((1/lineage$SUM)*lineage$proportion)*100


# Remove the assignments with a proportion of 0 after rounding
lineage <- lineage %>% filter(proportion > 0)


#################################################################################               
### Export data
#write.table(data.lineages, "./freyja_bootstraps_lineages.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(lineage[,c("samples", "Lineage", "proportion")], "./freyja_lineage.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
