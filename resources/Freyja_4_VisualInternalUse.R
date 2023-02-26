rm(list=ls(all=FALSE))
ls()
set.seed(123)


library(openxlsx)
library(data.table)
library(lubridate)
library(tidyverse)





########### Generate data ########### 

# Import data
freyja <-read.table("freyja_plot.tsv", header = TRUE, sep = "\t")
freyja.raw <-read.table("freyja_lineages_conversion_backup.tsv", header = TRUE, sep = "\t")
samplesinfo<-read.xlsx("ListSamples.xlsx", sheet="Samples", startRow = 3)
runsinfo<-read.xlsx("ListSamples.xlsx", sheet="Runs", startRow = 10)

# Merge three datasets
runsinfo<-runsinfo[, c("Run", "QC_Run")]
samplesinfo <- left_join(samplesinfo, runsinfo, by=c("Run"))
samplesinfo$samples<-paste0(samplesinfo$Run, "@", samplesinfo$FilesNames)

# Set (new) variables
samplesinfo <- samplesinfo %>% rename(sites = Location)
samplesinfo$Date<-as.Date(as.numeric(as.character(samplesinfo$Date)), origin = "1899-12-30")
samplesinfo$Date[samplesinfo$Sample.type=="Control"] <- Sys.Date()  # Add "today's date" to all Controls


# Select "Samples" 
samplesinfo<-samplesinfo %>% filter(samples != "NA")
samplesinfo$Month<-format(as.Date(samplesinfo$Date), "%b %Y")
samplesinfo$month<-format(as.Date(samplesinfo$Date), "%m %Y")
samplesinfo$week<-as.numeric(as.character(strftime(samplesinfo$Date, format = "%V")))
samplesinfo$year<-format(as.Date(samplesinfo$Date), "%Y")
samplesinfo$weekID<-ifelse((samplesinfo$week %% 2) == 0, samplesinfo$week-1, samplesinfo$week)

# Prepare plots

# Summarized lineages
freyja <- left_join(freyja, samplesinfo, by=c("samples"))
freyja<-freyja %>% drop_na(Run)
freyja$Display<-paste0(as.Date(freyja$Date), "_", freyja$samples)
last.run<-levels(as.factor(freyja$Run))[length(levels(as.factor(freyja$Run)))]
# Raw output
freyja.raw <- left_join(freyja.raw, samplesinfo, by=c("samples"))
freyja.raw<-freyja.raw %>% drop_na(Run)
freyja.raw$Display<-paste0(as.Date(freyja.raw$Date), "_", freyja.raw$samples)



save(freyja.raw, freyja, last.run, file="InternalUseData.RData")
