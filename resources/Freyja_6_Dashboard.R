rm(list=ls(all=FALSE))
ls()
set.seed(123)



library(openxlsx)
library(data.table)
library(lubridate)
library(ISOweek)
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
wwtp.info.levels<-read.table("wwtp_info.tsv", h=T, sep = "\t")


################################################################################
######## Preparation 
################################################################################

# Merge run and sample metadata0
runsinfo<-runsinfo[, c("Run", "QC_Run")]
samplesinfo <- dplyr::left_join(samplesinfo, runsinfo, by=c("Run"))


# Set new variables
samplesinfo <- samplesinfo %>% 
  dplyr::rename(sites = Location) %>%
  dplyr::mutate(samples = paste0(Run, "@", FilesNames), 
                Samples = paste0(sites, "_", Date),
                
                Date = as.Date(as.numeric(as.character(Date)), origin = "1899-12-30"),
                Date = as.Date(Date, "%m%d%y"),
                week = as.numeric(as.character(strftime(Date, format = "%V"))),
                #week = ISOweek::ISOweek(Date),
                
                # Split dates by month - MAP
                Month = format(as.Date(Date), "%b %Y"),
                Month_num = format(as.Date(Date), "%m %Y"),
                
                # Split dates by week - HEATMAP
                Week = ISOweek::ISOweek(Date),
                temp.Date_Week = ifelse(!is.na(Week),(paste0(Week, "-1")), NA), # Somehow this line is important otherwise it crashes the line below
                Date_Week = ISOweek::ISOweek2date(temp.Date_Week),
                
                # Split dates bi-weekly - BARPLOT
                biWeek = ifelse((week %% 2) == 0, ISOweek::ISOweek(Date-7), ISOweek::ISOweek(Date)), # Bi-weekly week, i.e., week 1 = 1, week 2 = 1, week 3 = 3, etc.  
                temp.Date_biWeek = ifelse(!is.na(biWeek),(paste0(biWeek, "-1")), NA), # Somehow this line is important otherwise it crashes the line below
                Date_biWeek = ISOweek::ISOweek2date(temp.Date_biWeek)) %>%
  
  dplyr::select(-temp.Date_Week, -temp.Date_biWeek) # Remove temporary columns



# Filter "Samples" 
samplesinfo<-samplesinfo %>% dplyr::filter(samples != "NA", Sample.type == "Sample", QC_Run == "qc_passed", QC_Sample == "qc_passed")
samplesinfo<-samplesinfo %>% dplyr::filter(sites %!like% "MadisonP")

# Merge Freyja and SampleInfo
freyja <- dplyr::left_join(freyja, samplesinfo, by=c("samples")) %>%
  tidyr::drop_na(Run) # Remove rows containing no Run info




################################################################################
######## Estimation of the SARS levels
################################################################################

# Preparation
data.SARS.level.UWM<-data.SARS.level.UWM |>
  dplyr::mutate(Date = as.Date(sample_collect_date, format="%Y-%m-%d")) |>
  dplyr::select(wwtp_name, sample_id, Date, pcr_target, pcr_target_avg_conc, flow_rate)
  
data.SARS.level.MDH<-data.SARS.level.MDH |>
  dplyr::mutate(Date = as.Date(sample_collect_date, format="%m/%d/%Y")) |>
  dplyr::select(wwtp_name, sample_id, Date, pcr_target, pcr_target_avg_conc, flow_rate)

data.SARS.level.WSLH<-data.SARS.level.WSLH |>
  dplyr::mutate(Date = as.Date(sample_collect_date_stop, format="%Y-%m-%d")) |>
  dplyr::filter(!grepl("not representative", tolower(wwtp_comments)), 
                pcr_target_avg_conc >= 0,
                sample_location == "wwtp") |>
  dplyr::select(wwtp_name, sample_id, Date, pcr_target, pcr_target_avg_conc, flow_rate)

# Merge MKE + WSLH SARS-levels datasets
data.SARS.level <- rbind(data.SARS.level.UWM, data.SARS.level.MDH, data.SARS.level.WSLH) |>
  dplyr::filter(pcr_target %in% c("sars-cov-2"), 
                Date > "2022-01-01") %>%
  dplyr::mutate(week = as.numeric(as.character(strftime(Date, format = "%V"))),
                biWeek = ifelse((week %% 2) == 0, ISOweek::ISOweek(Date-7), ISOweek::ISOweek(Date))) # Bi-weekly week, i.e., week 1 = 1, week 2 = 1, week 3 = 3, etc.  

# Add population size
wwtp.info.levels <- wwtp.info.levels |>
  dplyr::select(wwtp_GenomicDashboard, wwtp_HorizonExtract, PopulationServed) |> 
  dplyr::rename(City = wwtp_GenomicDashboard,
                wwtp_name = wwtp_HorizonExtract) |>
  distinct()

data.SARS.level <- left_join(data.SARS.level, wwtp.info.levels, by=c("wwtp_name"))


# Aggregate the data bi-weekly for BARPLOT visual
data.SARS.level.summarized <- data.SARS.level %>%
  dplyr::mutate_at(c('flow_rate', 'pcr_target_avg_conc', 'PopulationServed'), as.numeric)  %>%
  dplyr::group_by(as.factor(sample_id), City, biWeek) %>%
  dplyr::summarise(sars = ((geometric_mean(pcr_target_avg_conc+1, na.rm = TRUE)-1) * mean(flow_rate) * 3.7854*1e6 )/ mean(PopulationServed),
                   .groups = 'drop')
data.SARS.level.summarized.state <- data.SARS.level.summarized %>%
  dplyr::group_by(biWeek) %>%
  dplyr::summarise(sars.per.week = geometric_mean(sars+1, na.rm = TRUE)-1,
                   .groups = 'drop') %>%
  dplyr::mutate(City = "All cities combined")

data.SARS.level.summarized.state[is.nan(data.SARS.level.summarized.state)] <- 0
data.SARS.level.summarized.cities <- data.SARS.level.summarized %>%
  dplyr::group_by(biWeek, City) %>%
  dplyr::summarise(sars.per.week = geometric_mean(sars+1, na.rm = TRUE)-1,
                   .groups = 'drop')
data.SARS.level.summarized.cities[is.nan(data.SARS.level.summarized.cities)] <- 0
data.SARS.level.summarized<-rbind(data.SARS.level.summarized.state, data.SARS.level.summarized.cities)



################################################################################
######## STATEWIDE: Aggregate data weekly and bi-weekly
################################################################################

### WEEK - HEATMAP ###
freyja.all.week <-as.data.frame(freyja %>% 
  select(proportion, Lineage, Week, Date_Week, legend.position) %>%
  group_by(Week, Date_Week, Lineage, legend.position) %>%
  summarise(Sum = sum(proportion, na.rm = TRUE),
            n = n(),
            .groups = 'drop'))

freyja.all.week <- as.data.frame(setDT(
  as.data.frame(freyja %>% 
  select(proportion, Week, Date_Week) %>%
  group_by(Week, Date_Week) %>%
  summarise(Total = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))
  )[freyja.all.week, on="Date_Week"])

freyja.all.week<-freyja.all.week %>%
  dplyr::mutate(proportion = round(Sum/Total*100,2),
                n<-Total/100)  %>%
  dplyr::rename(Date = Date_Week)


### BI-WEEKLY - BARPLOT ###
freyja.all.biweekly <-as.data.frame(freyja %>% 
  select(proportion, Lineage, biWeek, Date_biWeek, legend.position) %>%
  group_by(biWeek, Date_biWeek, Lineage, legend.position) %>%
  summarise(Sum = sum(proportion, na.rm = TRUE),
            n = n(),
            .groups = 'drop'))

freyja.all.biweekly <- as.data.frame(setDT(
  as.data.frame(freyja %>% 
  select(proportion, biWeek, Date_biWeek) %>%
  group_by(biWeek, Date_biWeek) %>%
  summarise(Total = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))
  )[freyja.all.biweekly, on="Date_biWeek"])

freyja.all.biweekly<-freyja.all.biweekly %>%
  dplyr::mutate(proportion=round(Sum/Total*100,2),
                City = "All cities combined",
                week.num = as.numeric(str_sub(biWeek, start= -2)),
                hoover = paste0("weeks ", week.num, "-", week.num+1, " (starting ", format(Date_biWeek, format="%b %d, %y"), ")<br>",
                                   Lineage, "<br>",
                                   round(proportion, 1), "% <br>Average of ",
                                   round(Total, 0), " sample(s)")) %>%
  dplyr::rename(Date = Date_biWeek)




################################################################################
######## SEWERSHED LEVEL: Aggregate data weekly and bi-weekly
################################################################################
freyja.cities<-freyja

### WEEK ###
freyja.cities$weekcity<-paste0(freyja.cities$Week, freyja.cities$sites)

freyja.city.week <-as.data.frame(freyja.cities %>% 
  select(proportion, Lineage, Week, Date_Week, weekcity, sites, legend.position) %>%
  group_by(weekcity, Lineage,legend.position) %>%
  summarise(Sum = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))

freyja.city.week <- as.data.frame(setDT(
  as.data.frame(freyja.cities %>% 
  select(proportion, Week, Date_Week, sites, weekcity, sites) %>%
  group_by(weekcity, sites, Week, Date_Week) %>%
  summarise(Total = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))
  )[freyja.city.week, on="weekcity"])

freyja.city.week<-freyja.city.week %>%
  dplyr::mutate(proportion = round(Sum/Total*100,2),
                n = Total/100)  %>%
  dplyr::rename(Date = Date_Week)


### BI-WEEKLY ###
freyja.cities$weekIDcity<-paste0(freyja.cities$biWeek, freyja.cities$sites)

freyja.city.biweekly <-as.data.frame(freyja.cities %>% 
  select(proportion, Lineage, biWeek, Date_biWeek, weekIDcity, sites, legend.position) %>%
  group_by(weekIDcity, Lineage, legend.position) %>%
  summarise(Sum = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))

freyja.city.biweekly <- as.data.frame(setDT(
  as.data.frame(freyja.cities %>% 
  select(proportion, biWeek, Date_biWeek, sites, weekIDcity, sites) %>%
  group_by(weekIDcity, sites, biWeek, Date_biWeek) %>%
  summarise(Total = sum(proportion, na.rm = TRUE),
            .groups = 'drop'))
  )[freyja.city.biweekly, on="weekIDcity"])

freyja.city.biweekly<-freyja.city.biweekly %>%
  dplyr::mutate(proportion=round(Sum/Total*100,2),
                n = Total/100,
                City = sites,
                week.num = as.numeric(str_sub(biWeek, start= -2)),
                hoover = paste0("weeks ", week.num, "-", week.num+1, " (starting ", format(Date_biWeek, format="%b %d, %y"), ")<br>",
                                Lineage, "<br>",
                                round(proportion, 1), "% <br>Average of ",
                                round(Total, 0), " sample(s)")) %>%
  dplyr::rename(Date = Date_biWeek)


################################################################################
######## VISUAL PREPARATION
################################################################################

colors<-c("#00425f", #Delta
          "#BDF271", #Omicron
          "#e53935", #BA.1
          "#0eeaff", #BA.2
          "#33691E", #BA.4
          "#26A69A", #BA.5
          "#f79100", #BA.2.12.1
          "#a2decc", #BA.2.75
          "#b9121b", #BQ.1
          "#f6e92b", #XBB
          "#04668C", #XBB.1.5
          "#04BFBF", #XBB.1.16
          "#553285", #CH.1
          "#ff1d23", #XBB.1.9
          "#900B0A", #XBB.2.3
          "#45BF55", #EG.5
          "#BEEB9F", #XBB.1.5.70
          "#020873", #HK.3
          "#9768d1", #BA.2.86
          "#ff822e", #JN.1
          "#ACF0F2", #JN.1.11.1
          "#ffd10f", #KP.3
          "#43A047", #XDV.1   
          "#0092B2", #KP3.1.1
          "#EF5350", #XEC
          "#35478C", #KP.2.3
          "#ffa600", #LF.7
          "#D50000", #MV.1
          "#96CA2D", #LP.8.1
          "#ff9b78",
          "#6b14a6",
          "#0288D1",
          "#ff2d00",
          "#0EEAFF",
          "#36175e",
          "#add5f7",
          "#1976d2",
          "#ffd34e",
          "#f77a52",
          "#ab47bc",
          "#3498db",
          "#f2b705",
          "#553285",
          "#168039",
          "#f4511e",
          "#00305a",
          "#ffd933",
          "#7abaf2",
          "#5c0002",
          "#f0c755",
          "#0d47a1")
colors.grey<-"#9db8c7"


lvl.lineage.color <- freyja %>%           
  dplyr::select(legend.position, Lineage) %>%
  unique() %>%
  dplyr::arrange(legend.position) %>% # sort dataframe according to legend.position order
  dplyr::mutate(Lineage = factor(Lineage, unique(Lineage)),  # reset your factor-column based on that order
                color = c(colors.grey, colors[1:length(Lineage)-1]))
  

################################################################################
######## VISUAL: BARPLOT
################################################################################
# Extract relevant data
freyja.barplot<-rbind(freyja.all.biweekly[, c("biWeek", "Lineage", "proportion", "Date", "City", "hoover", "legend.position")], 
                      freyja.city.biweekly[, c("biWeek", "Lineage", "proportion", "Date",  "City", "hoover", "legend.position")])

# Add SARS data
freyja.barplot<-dplyr::left_join(freyja.barplot, data.SARS.level.summarized, by=c("biWeek", "City")) %>%
  dplyr::arrange(legend.position) %>%               # sort dataframe according to legend.position order
  dplyr::mutate(Lineage = factor(Lineage, unique(Lineage))) # reset your factor-column based on that order

# Define colors
colors.plot.barplot<-lvl.lineage.color %>%
  dplyr::filter(Lineage %in% freyja.barplot$Lineage) %>%
  dplyr::pull(color) # Only keep column 'color' and convert to vector



################################################################################
######## VISUAL: HEATMAP
################################################################################
# Extract relevant data
freyja.heatmap<-freyja.city.week

# Determine the predominant variant in a sample
freyja.heatmap$predominance<-freyja.heatmap$Lineage
temp<-freyja.heatmap %>%
  dplyr::group_by(sites, Week) %>%
  dplyr::filter(proportion == max(proportion, na.rm=TRUE), 
                Date > Sys.Date()-2.1*365) %>% ############### IF THE GRADIENT GLITCH HAPPENS AGAIN, change 2.1 with lower value (e.g., 2, 1.9, 1.8)
  dplyr::mutate(predominance = "Predominant variants",
                week.num = as.numeric(str_sub(Week, start= -2)), 
                tooltip = paste0(Lineage, "<br>",
                     sites, "<br>",
                     "week ", week.num, " (", Date, ")")) %>%
  as.data.frame()


# Add 0% when lineage not observed in a sample 
level.lineages.heatmap<-levels(as.factor(freyja.heatmap$Lineage))
temp.1<-temp %>%
  dplyr::mutate(Lineage = NA, 
                predominance = NA, 
                tooltip = NA, 
                proportion = 0, 
                Sum = 0)
for (i in 1:length(level.lineages.heatmap)){
  temp.2<-freyja.heatmap %>% dplyr::filter(Lineage == level.lineages.heatmap[i])
  temp.3<-temp.1[which(temp.1$weekcity %!in% temp.2$weekcity), c("weekcity", "sites", "Week", "Total", "Lineage", "Sum", "proportion", "Date", "n", "predominance", "legend.position")]
  temp.3$Lineage = temp.3$predominance <-level.lineages.heatmap[i]
  freyja.heatmap<-rbind(freyja.heatmap, temp.3)
}

# # Merge predominant data with all lineages (including 0%) data: merge freyja.heatmap and temp 
freyja.heatmap <- freyja.heatmap %>%
  dplyr::mutate(week.num = as.numeric(str_sub(Week, start= -2)),
                tooltip = paste0(Lineage, "<br>",
                                 "Relative abundance: ", round(proportion, 1), "%", "<br>",
                                 sites, "<br>",
                                 "week ", week.num, " (", Date, ")"),
                Date = as.POSIXct(Date, format="%Y-%m-%d")) %>%
  dplyr::add_row(temp)



# Add population served
freyja.heatmap <- dplyr::left_join(freyja.heatmap, wwtp.info, by=c("sites"))


# Last edits to freyja.heatmap for dashboard
freyja.heatmap <- freyja.heatmap %>% 
  dplyr::filter(Lineage != "Other") %>%  #Remove from the rolling menu the "Other" Lineage
  dplyr::mutate(`Relative abundance (%)` = proportion) %>% # Rename 'proportion'
  dplyr::arrange(legend.position) %>%               # sort dataframe according to legend.position order
  dplyr::mutate(Lineage = factor(Lineage, unique(Lineage))) %>% # reset your factor-column based on that order
  dplyr::select(sites, PopulationServed, Lineage, `Relative abundance (%)`, tooltip, predominance, Date) # Select variables of interest


# Define colors
colors.plot.heatmap<-lvl.lineage.color %>%
  dplyr::arrange(desc(legend.position)) %>%
  dplyr::filter(Lineage %in% as.factor(temp$Lineage),
                Lineage != "Other") %>% # Remove thecolor associated to Other
  dplyr::pull(color) # Only keep column 'color' and convert to vector





################################################################################
######## VISUAL: MAP
################################################################################
zipcode<- wwtp.info %>% select(sites, zipcode) %>% filter(zipcode != "NA")


# # Define colors - THIS DOESN"T WORK!!
# colors.plot.map<-lvl.lineage.color %>%
#   dplyr::filter(Lineage %in% as.factor(names(freyja.map)),
#                 Lineage != "Other") %>% # Remove thecolor associated to Other
#   dplyr::arrange(-desc(Lineage)) %>% ## DOES NOT WORK!!!!!
#   dplyr::mutate(Lineage = factor(Lineage, unique(Lineage))) %>%
#   dplyr::pull(color)

# Define colors - somehow, arrange() was not working (above), so I had to find a way around to make it work :-|
colors.plot.map <- data.frame(Lineage = levels(as.factor(freyja$Lineage))) %>% # List Lineages in alphabetical order!!
  dplyr::filter(Lineage != "Other") # Remove the color associated to Other
colors.plot.map <- dplyr::left_join(colors.plot.map, lvl.lineage.color, by="Lineage") %>% # Add HEX to lineage
  dplyr::pull(color) # extract colors



# Extract relevant data
freyja.map<-freyja %>% 
  dplyr::filter(Lineage != "Other")

# Convert zipcode into long and lat
geo<-as.data.frame(geocode_zip(zipcode$zipcode))
zipcode$zipcode<-as.character(zipcode$zipcode)
zipcode<-dplyr::full_join(geo, zipcode, by="zipcode")
names(zipcode)<-c("zipcode", "lat", "long", "sites")
freyja.map<-left_join(freyja.map, zipcode, by="sites")

## Extract first occurence of each variant
freyja.1st.occurence <- freyja.map %>%
  dplyr::filter(proportion > 0) %>%
  dplyr::select(Lineage, Date) %>%
  dplyr::group_by(Lineage) %>%
  dplyr::top_n(n=1, wt=desc(Date)) %>%
  as.data.frame()

### Dashboard
#freyja.map$Week <- paste0(freyja.map$Month, " week ", freyja.map$week)
freyja.map<-freyja.map %>%
  dplyr::select(sites, Date, Lineage, proportion, lat, long, Month, Week, legend.position)

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
  dplyr::group_by(sites, Lineage, lat, long, Month) %>%
  dplyr::summarise(abundance=round(sum(proportion*100), 2),
            n = n(),
            .groups = 'drop') %>%
  dplyr::mutate(proportion = round(abundance/n, digits=2))

# Format df for dashboard
freyja.map<-freyja.map %>%
  dplyr::select(Lineage, sites, lat, long, Month, proportion) %>%
  tidyr::spread(key=Lineage, value=proportion, fill=0) %>%
  dplyr::select(-"NA") %>%              
  dplyr::mutate(Month = as.Date(paste0("01 ", Month), format="%d %b %Y")) 




################################################################################
######## VISUAL: TIMESTAMP
################################################################################
TimeStamp<-as.character(Sys.Date())




################################################################################
######## EXPORT
################################################################################
save(freyja.barplot, colors.plot.barplot, freyja.heatmap, colors.plot.heatmap, freyja.map, colors.plot.map, TimeStamp, freyja.1st.occurence, file = "DashboardData.RData")

