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
runsinfo<-openxlsx::read.xlsx("ListSamples.xlsx", sheet="Runs", startRow = 10)

# Merge three datasets
runsinfo<-runsinfo[, c("Run", "QC_Run")]
samplesinfo <- dplyr::left_join(samplesinfo, runsinfo, by=c("Run")) %>%
  dplyr::mutate(samples = paste0(Run, "@", FilesNames))

# Set (new) variables
samplesinfo <- samplesinfo %>% 
  dplyr::rename(sites = Location) %>%
  dplyr::mutate(Date = as.Date(as.numeric(as.character(Date)), origin = "1899-12-30"))
samplesinfo$Date[samplesinfo$Sample.type=="Control"] <- Sys.Date()  # Add "today's date" to all Controls


# Select "Samples" 
samplesinfo<-samplesinfo %>% 
  dplyr::filter(samples != "NA") %>%
  dplyr::mutate(Month = format(as.Date(Date), "%b %Y"),
                month = format(as.Date(Date), "%m %Y"),
                week = as.numeric(as.character(strftime(Date, format = "%V"))),
                year = format(as.Date(Date), "%Y"),
                weekID = ifelse((week %% 2) == 0, week-1, week))

# Prepare plots

# Summarized lineages
freyja <- dplyr::left_join(freyja, samplesinfo, by=c("samples")) %>%
  tidyr::drop_na(Run) %>%
  dplyr::mutate(Display = paste0(as.Date(Date), "_", samples))
last.run<-levels(as.factor(freyja$Run))[length(levels(as.factor(freyja$Run)))]
# Raw output
freyja.raw <- dplyr::left_join(freyja.raw, samplesinfo, by=c("samples")) %>% 
  tidyr::drop_na(Run) %>%
  dplyr::mutate(Display = paste0(as.Date(Date), "_", samples))



save(freyja.raw, freyja, last.run, file="InternalUseData.RData")
