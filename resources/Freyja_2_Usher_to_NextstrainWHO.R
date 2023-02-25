rm(list=ls(all=FALSE))
ls()
character(0)

library(data.table)
library(stringr)
library(dplyr)
library(jsonlite)


# Convert json to a df
jdata <- read_json("nameTable.json", simplifyVector = TRUE)
jdata <- as.data.frame(jdata)
for(i in 1:nrow(jdata)){
  temp<-0
  temp<-cbind(jdata[i,1], ifelse(length(jdata[[2]][[i]]$name)!=0, jdata[[2]][[i]]$name, ""))
  if(i==1){
    NextstrainWHO<-temp
  } else {
    NextstrainWHO<-rbind(temp, NextstrainWHO)
  }
}

# Prepare NextstrainWHO
NextstrainWHO<-as.data.frame(NextstrainWHO)
names(NextstrainWHO)<-c("Nextstrain", "Pango")
NextstrainWHO$WHO<-gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", NextstrainWHO$Nextstrain, perl=T)
NextstrainWHO$Nextstrain<-sub(" .*", "", NextstrainWHO$Nextstrain)
write.table(NextstrainWHO, "NextstrainWHO.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


