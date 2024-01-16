rm(list=ls(all=FALSE))
ls()
set.seed(123)



library(openxlsx)
library(data.table)
library(lubridate)
library(tidyverse)
library(hrbrthemes) #sudo apt -y install libfontconfig1-dev; sudo apt-get install libcairo2-dev
library(sf)
library(zipcodeR)


'%!in%' <- function(x,y)!('%in%'(x,y)) 
'%!like%' <- function(x,y)!('%like%'(x,y)) 
'geometric_mean' <- function(x,na.rm=TRUE) { exp(mean(log(x),na.rm=na.rm)) }
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))




################################################################################
######## Import data 
################################################################################

# Import sequencing related data
freyja <-read.table("freyja_plot.tsv", header = TRUE, sep = "\t")
samplesinfo<-read.xlsx("ListSamples.xlsx", sheet="Samples", startRow = 3)
runsinfo<-read.xlsx("ListSamples.xlsx", sheet="Runs", startRow = 10)
wwtp.info<-read.table(url("https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/data/wwtp_info.txt"), h=T, sep = "\t")

# Import SARS-CoV-2 loads related data
load("./SARSCoV2concentrationData.RData")

################################################################################
######## Preparation 
################################################################################

# Merge run and sample metadata0
runsinfo<-runsinfo[, c("Run", "QC_Run")]
samplesinfo <- left_join(samplesinfo, runsinfo, by=c("Run"))
samplesinfo$samples<-paste0(samplesinfo$Run, "@", samplesinfo$FilesNames)


# Set (new) variables
samplesinfo <- samplesinfo %>% rename(sites = Location)
samplesinfo$Date<-as.Date(as.numeric(as.character(samplesinfo$Date)), origin = "1899-12-30")
samplesinfo$Date<-as.Date(samplesinfo$Date, "%m%d%y")
samplesinfo$Samples<-paste0(samplesinfo$sites, "_", samplesinfo$date)


# Select "Samples" 
samplesinfo<-samplesinfo %>% filter(samples != "NA", Sample.type == "Sample", QC_Run == "qc_passed", QC_Sample == "qc_passed")
samplesinfo<-samplesinfo %>% filter(sites %!like% "MadisonP")

samplesinfo$Month<-format(as.Date(samplesinfo$Date), "%b %Y")
samplesinfo$month<-format(as.Date(samplesinfo$Date), "%m %Y")
samplesinfo$week<-as.numeric(as.character(strftime(samplesinfo$Date, format = "%V")))
samplesinfo$year<-format(as.Date(samplesinfo$Date), "%Y")
samplesinfo$WEEK<-paste0(samplesinfo$year, "_", strftime(samplesinfo$Date, format = "%V"))
samplesinfo$WEEKID<-ifelse((samplesinfo$week %% 2) == 0, 
                           ifelse(nchar(samplesinfo$week-1)==1, 
                                  paste0(samplesinfo$year, "_0", samplesinfo$week-1),
                                  paste0(samplesinfo$year, "_", samplesinfo$week-1)),
                           ifelse(nchar(samplesinfo$week)==1, 
                                  paste0(samplesinfo$year, "_0", samplesinfo$week),
                                  paste0(samplesinfo$year, "_", samplesinfo$week)))
samplesinfo$weekID<-ifelse((samplesinfo$week %% 2) == 0, samplesinfo$week-1, samplesinfo$week)

# Merge Freyja and SampleInfo
freyja <- left_join(freyja, samplesinfo, by=c("samples"))
freyja<-freyja %>% drop_na(Run)




################################################################################
######## Estimation of the SARS levels
################################################################################

# Preparation
data.SARS.level.MKE<-data.SARS.level.MKE |>
  dplyr::mutate(Date = as.Date(sample_collect_date, format="%m/%d/%Y")) |>
  dplyr::select(wwtp_name, sample_id, Date, pcr_gene_target, pcr_target_avg_conc, flow_rate)
  
data.SARS.level.WSLH<-data.SARS.level.WSLH |>
  dplyr::mutate(Date = as.Date(sample_collect_date_stop, format="%Y-%m-%d")) |>
  dplyr::filter(!grepl("not representative", tolower(wwtp_comments)), 
                pcr_target_avg_conc >= 0,
                sample_location == "wwtp") |>
  dplyr::select(wwtp_name, sample_id, Date, pcr_gene_target, pcr_target_avg_conc, flow_rate)

# Merge MKE + WSLH SARS-levels datasets
data.SARS.level <- rbind(data.SARS.level.MKE, data.SARS.level.WSLH) |>
  dplyr::filter(pcr_gene_target %in% c("n1", "n2")) %>%
  dplyr::mutate(week = as.numeric(as.character(strftime(Date, format = "%V"))),
                weekID = ifelse((week %% 2) == 0, week-1, week),
                year = format(as.Date(Date), "%Y"))

# Add population size
wwtp.info.levels <- wwtp.info.levels |>
  dplyr::select(wwtp_GenomicDashboard, wwtp_HorizonExtract, PopulationServed) |> 
  dplyr::rename(City = wwtp_GenomicDashboard,
                wwtp_name = wwtp_HorizonExtract) |>
  distinct()

data.SARS.level <- left_join(data.SARS.level, wwtp.info.levels, by=c("wwtp_name"))


# Aggregate the data bi-weekly
data.SARS.level.summarized <- data.SARS.level %>%
  dplyr::group_by(as.factor(sample_id), City, weekID, year) %>%
  dplyr::summarise(sars = (geometric_mean(pcr_target_avg_conc+1, na.rm = TRUE) * mean(flow_rate) * 3.7854*1e6 )/ mean(PopulationServed),
                   .groups = 'drop')
data.SARS.level.summarized.state <- data.SARS.level.summarized %>%
  dplyr::group_by(weekID, year) %>%
  dplyr::summarise(sars.per.week = geometric_mean(sars+1, na.rm = TRUE),
                   .groups = 'drop') %>%
  dplyr::mutate(City = "All cities combined")

data.SARS.level.summarized.state[is.nan(data.SARS.level.summarized.state)] <- 0
data.SARS.level.summarized.cities <- data.SARS.level.summarized %>%
  dplyr::group_by(weekID, year, City) %>%
  dplyr::summarise(sars.per.week = geometric_mean(sars+1, na.rm = TRUE),
                   .groups = 'drop')
data.SARS.level.summarized.cities[is.nan(data.SARS.level.summarized.cities)] <- 0
data.SARS.level.summarized<-rbind(data.SARS.level.summarized.state, data.SARS.level.summarized.cities)



################################################################################
######## STATEWIDE: Aggregate data per week/bi-weekly
################################################################################

### WEEK ###
freyja.all.week <-as.data.frame(freyja %>% 
  select(proportion, Lineage, WEEK, week, year) %>%
  group_by(week, WEEK, Lineage) %>%
  summarise(Sum = sum(proportion, na.rm = TRUE),
            n = n(),
            .groups = 'drop'))

freyja.all.week <- as.data.frame(setDT(
  as.data.frame(freyja %>% 
  select(proportion, year, week, WEEK) %>%
  group_by(week, WEEK, year) %>%
  summarise(Total = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))
  )[freyja.all.week, on="WEEK"])

freyja.all.week$proportion<-round(freyja.all.week$Sum/freyja.all.week$Total*100,2)
freyja.all.week$Date<-parse_date_time(paste(freyja.all.week$year, freyja.all.week$week, 1, sep="/"),'Y/W/w')
freyja.all.week$n<-freyja.all.week$Total/100


### BI-WEEKLY ###
freyja.all.biweekly <-as.data.frame(freyja %>% 
  select(proportion, Lineage, year, weekID, WEEKID) %>%
  group_by(WEEKID, weekID, Lineage) %>%
  summarise(Sum = sum(proportion, na.rm = TRUE),
            n = n(),
            .groups = 'drop'))

freyja.all.biweekly <- as.data.frame(setDT(
  as.data.frame(freyja %>% 
  select(proportion, year, weekID, WEEKID) %>%
  group_by(WEEKID, weekID, year) %>%
  summarise(Total = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))
  )[freyja.all.biweekly, on="WEEKID"])

freyja.all.biweekly$proportion<-round(freyja.all.biweekly$Sum/freyja.all.biweekly$Total*100,2)
freyja.all.biweekly$Date<-parse_date_time(paste(freyja.all.biweekly$year, freyja.all.biweekly$weekID, 1, sep="/"),'Y/W/w')

freyja.all.biweekly$City<-"All cities combined"
freyja.all.biweekly$hoover<-paste0("weeks ", freyja.all.biweekly$weekID, "-", freyja.all.biweekly$weekID+1, " (starting ", format(freyja.all.biweekly$Date, format="%b %d, %y"), ")<br>",
                                   freyja.all.biweekly$Lineage, "<br>",
                                   round(freyja.all.biweekly$proportion, 1), "% <br>Average of ",
                                   round(freyja.all.biweekly$Total, 0), " sample(s)")




################################################################################
######## SEWERSHED LEVEL: Aggregate data per week/bi-weekly
################################################################################
freyja.cities<-freyja

### WEEK ###
freyja.cities$weekcity<-paste0(freyja.cities$week, freyja.cities$sites, freyja.cities$year)

freyja.city.week <-as.data.frame(freyja.cities %>% 
  select(proportion, Lineage, week, weekcity, sites, year) %>%
  group_by(weekcity, Lineage) %>%
  summarise(Sum = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))

freyja.city.week <- as.data.frame(setDT(
  as.data.frame(freyja.cities %>% 
  select(proportion, year, sites, weekcity, week, sites) %>%
  group_by(weekcity, sites, year, week) %>%
  summarise(Total = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))
  )[freyja.city.week, on="weekcity"])

freyja.city.week$proportion<-round(freyja.city.week$Sum/freyja.city.week$Total*100,2)
freyja.city.week$Date<-parse_date_time(paste(freyja.city.week$year, freyja.city.week$week, 1, sep="/"),'Y/W/w')
freyja.city.week$n<-freyja.city.week$Total/100


### BI-WEEKLY ###
freyja.cities$weekIDcity<-paste0(freyja.cities$weekID, freyja.cities$sites, freyja.cities$year)

freyja.city.biweekly <-as.data.frame(freyja.cities %>% 
  select(proportion, Lineage, weekID, weekIDcity, sites, year) %>%
  group_by(weekIDcity, Lineage) %>%
  summarise(Sum = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))

freyja.city.biweekly <- as.data.frame(setDT(
  as.data.frame(freyja.cities %>% 
  select(proportion, year, sites, weekIDcity, weekID, sites) %>%
  group_by(weekIDcity, sites, year, weekID) %>%
  summarise(Total = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))
  )[freyja.city.biweekly, on="weekIDcity"])

freyja.city.biweekly$proportion<-round(freyja.city.biweekly$Sum/freyja.city.biweekly$Total*100,2)
freyja.city.biweekly$Date<-parse_date_time(paste(freyja.city.biweekly$year, freyja.city.biweekly$weekID, 1, sep="/"),'Y/W/w')
freyja.city.biweekly$n<-freyja.city.biweekly$Total/100

freyja.city.biweekly$City<-freyja.city.biweekly$sites
freyja.city.biweekly$hoover<-paste0("weeks ", freyja.city.biweekly$weekID, "-", freyja.city.biweekly$weekID+1, " (starting ", format(freyja.city.biweekly$Date, format="%b %d, %y"), ")<br>",
                                    freyja.city.biweekly$Lineage, "<br>",
                                    round(freyja.city.biweekly$proportion, 1), "% <br>Average of ",
                                    round(freyja.city.biweekly$Total, 0), " sample(s)")




################################################################################
######## VISUAL PREPARATION
################################################################################
colors<-c('#002f45','#dcbeff','#6bade3',"#872b8b",'#9A6324','#64b34d','#fb9a99','#ffc500','#1e90ff',"#e31a1c",'#f58231','#98fb98',"#6a3d9a",'#800000','#ccebc5','#75fff1','#ffff54','#cab2d6','#fdbf6f','#ff1493')
colors.grey<-"#9db8c7"



  

################################################################################
######## VISUAL: BARPLOT
################################################################################
freyja.barplot<-rbind(freyja.all.biweekly[, c("weekID", "year", "Lineage", "proportion", "Date", "City", "hoover")], 
                      freyja.city.biweekly[, c("weekID", "year", "Lineage", "proportion", "Date",  "City", "hoover")])

old.lvl<-levels(as.factor(freyja.barplot$Lineage))
freyja.barplot$Lineage<-factor(freyja.barplot$Lineage, levels=c("Other", sort(old.lvl[old.lvl!="Other"], decreasing=T)))
freyja.barplot<-dplyr::left_join(freyja.barplot, data.SARS.level.summarized, by=c("weekID", "year", "City"))

colors.plot<-c(colors.grey, rev(colors[1:length(old.lvl)-1]))
#scales::show_col(colors.plot)
colors.plot.barplot<-colors.plot


################################################################################
######## VISUAL: HEATMAP
################################################################################
freyja.heatmap<-freyja.city.week
freyja.heatmap$alpha<-freyja.heatmap$proportion/100

# Determine the predominant variant in a sample
freyja.heatmap$predominance<-freyja.heatmap$Lineage
temp<-freyja.heatmap %>%
  group_by(sites, week) %>%
  filter(proportion == max(proportion, na.rm=TRUE))
temp$predominance<-"Predominant variants"
temp$alpha<-1
temp$tooltip<-paste0(temp$Lineage, "<br>",
                     temp$sites, "<br>",
                     "week ", temp$week, " (", temp$Date, ")")


# Add 0% when lineage not observed in a sample 
level.lineages<-levels(as.factor(freyja.heatmap$Lineage))
temp.1<-as.data.frame(temp)
temp.1$Lineage = temp.1$predominance = temp.1$tooltip <-"NA"; temp.1$proportion = temp.1$Sum = temp.1$alpha <-0;
for (i in 1:length(level.lineages)){
  temp.2<-freyja.heatmap %>% filter(Lineage == level.lineages[i])
  temp.3<-temp.1[which(temp.1$weekcity %!in% temp.2$weekcity), c("weekcity", "sites", "year", "week", "alpha", "Total", "Lineage", "Sum", "proportion", "Date", "n", "predominance")]
  temp.3$Lineage = temp.3$predominance <-level.lineages[i]
  freyja.heatmap<-rbind(freyja.heatmap, temp.3)
}
freyja.heatmap$tooltip<-paste0(freyja.heatmap$Lineage, "<br>",
                     "Relative abundance: ", round(freyja.heatmap$proportion, 1), "%", "<br>",
                     freyja.heatmap$sites, "<br>",
                     "week ", freyja.heatmap$week, " (", freyja.heatmap$Date, ")")

# Merge  predominant data with all lineages (including 0%) data
freyja.heatmap<-rbind(freyja.heatmap, temp)

freyja.heatmap <- freyja.heatmap %>% filter(Lineage != "Other")
freyja.heatmap$`Relative abundance (%)`<-freyja.heatmap$proportion
freyja.heatmap$sites<-factor(freyja.heatmap$sites, levels = c(sort(unique(freyja.heatmap$site), decreasing=T)))
freyja.heatmap.alphabetical<-freyja.heatmap

# Define colors
old.lvl<-levels(as.factor(freyja.heatmap$Lineage))
freyja.heatmap$Lineage<-factor(freyja.heatmap$Lineage, levels=c(sort(old.lvl, decreasing=F)))
colors.plot<-(c(colors[1:length(old.lvl)]))
colors.plot.heatmap<-colors.plot
colors.heatmap<-as.data.frame(cbind(colors.plot.heatmap, old.lvl))
names(colors.heatmap)<-c("colors", "Lineage")
colors.heatmap <- colors.heatmap %>% filter(Lineage %in% unique(temp$Lineage))
#colors.plot.heatmap<-colors.heatmap$colors
for (i in 1:nrow(colors.heatmap)){
  variant=colors.heatmap[i, 2]
  hex=colors.heatmap[i, 1]
  pre<-assign(variant, hex)
  if(i == 1){colors.plot.heatmap <-pre} else {colors.plot.heatmap<-c(pre, colors.plot.heatmap)}
}

# Add population served
freyja.heatmap <- left_join(freyja.heatmap, wwtp.info, by=c("sites"))

#Setting up for dashboard
freyja.heatmap<- freyja.heatmap %>%
  select(sites, PopulationServed, Lineage, `Relative abundance (%)`, tooltip, predominance, Date)
freyja.heatmap$sites<-as.character(freyja.heatmap$sites)
colors.plot.heatmap<-rev(colors.plot.heatmap)




################################################################################
######## VISUAL: MAP
################################################################################
zipcode<- wwtp.info %>% select(sites, zipcode) %>% filter(zipcode != "NA")

freyja.map<-freyja 
freyja.map <- freyja.map %>% filter(Lineage != "Other")

geo<-as.data.frame(geocode_zip(zipcode$zipcode))
zipcode$zipcode<-as.character(zipcode$zipcode)
zipcode<-full_join(geo, zipcode, by="zipcode")
names(zipcode)<-c("zipcode", "lat", "long", "sites")

freyja.map<-left_join(freyja.map, zipcode, by="sites")
old.lvl<-levels(as.factor(freyja.map$Lineage))
freyja.map$Lineage<-factor(freyja.map$Lineage, levels=c(sort(old.lvl, decreasing=F)))
colors.plot<-(c(colors[1:length(old.lvl)]))
colors.plot.map<-colors.plot

## Extract first occurence of each variant
freyja.1st.occurence <- freyja.map %>%
  filter(proportion > 0) %>%
  select(Lineage, Date) %>%
  group_by(Lineage) %>%
  top_n(n=1, wt=desc(Date))
freyja.1st.occurence<-as.data.frame(freyja.1st.occurence)

### Dashboard
freyja.map$Week <- paste0(freyja.map$Month, " week ", freyja.map$week)
freyja.map<-freyja.map %>%
  select(sites, Date, Lineage, proportion, lat, long, Month, Week)

## Add "missing" samples = Attribute a value for each collection date (even if everything = 0)
freyja.map.empty<-freyja.map
freyja.map.empty$sites = freyja.map.empty$Lineage = freyja.map.empty$lat = freyja.map.empty$long = "NA"; freyja.map.empty$proportion = 0
freyja.map.empty<-unique(freyja.map.empty)
level.sites<-levels(as.factor(freyja.map$sites))

for (i in 1:length(level.sites)){
  temp.2<-freyja.map %>% filter(sites == level.sites[i])
  temp.3<-freyja.map.empty[which(freyja.map.empty$Date %!in% temp.2$Date), ]
  temp.3$sites <-level.sites[i]; temp.3$lat<-temp.2[1, c("lat")]; temp.3$long<-temp.2[1, c("long")]
  freyja.map<-rbind(freyja.map, temp.3)
}

freyja.map <- freyja.map %>%
  group_by(sites, Lineage, lat, long, Month) %>%
  summarise(abundance=round(sum(proportion*100), 2),
            n = n(),
            .groups = 'drop') 
freyja.map$proportion<-round(freyja.map$abundance/freyja.map$n, digits=2)

# Format df for dashbaord
freyja.map<-reshape2::dcast(freyja.map, sites+lat+long+Month~Lineage, sum, value.var = "proportion")
freyja.map<-freyja.map[ , -which(names(freyja.map) %in% c("NA"))]
freyja.map$Month<-as.Date(paste0("01 ", freyja.map$Month), format="%d %b %Y")




################################################################################
######## VISUAL: TIMESTAMP
################################################################################
TimeStamp<-as.character(Sys.Date())




################################################################################
######## EXPORT
################################################################################
save(freyja.barplot, colors.plot.barplot, freyja.heatmap, colors.plot.heatmap, freyja.map, colors.plot.map, TimeStamp, freyja.1st.occurence, file = "DashboardData.RData")

