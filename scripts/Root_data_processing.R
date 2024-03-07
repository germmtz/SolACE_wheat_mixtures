#### G. Montazeaud
#### 13.09.2023

### This script aim to process raw root trait data obtained from the 4PMi highthroughput phenotyping platform. 

## Package loading
library(lme4)

## Data loading
dates <- c("28_06_2019","07_07_2019","15_07_2019") # root traits were measured at three different dates 
lots <- c("Calib01","Calib02","Calib03","Calib04","Clim","Rec") # there are three types of samples in the data set: calibration samples (used to calibrate image analysis), climatic samples (used to assess the effect of the local climate on plant growth), and the harvest ("Rec" for "Récolte") data which corresponds to the main experimental data. 
har <-  read.table("data/raw_data/biomass_data/harvest_data.csv", header=T, sep=";") # harvest data is used to have the genotype IDs for each RT

## We first load a list of all rhizotubes for each type of samples (calibration samples, climatic samples, and harvest samples)
ref <- read.table("data/raw_data/List_rhyzotubes.csv", header=T, sep=",")
RT_Rec <- ref[which(ref$Lot=="Rec"),"Num"]

all_root <- data.frame()

for (d in dates) {
  root <- read.table(paste("data/raw_data/root_data/roots_",d, ".csv", sep = ""), header = T, sep = ";")
  
  ### keeping only rhizotubes from the "Rec" lot
  root <- root[which(root$plantid%in%RT_Rec),]
  
  ### removing doublons
  if(d=="28_06_2019"){
    root <- root[-which(root$plantid==30 &root$taskid==19285),]
    root <- root[-which(root$plantid==573 &root$taskid==19291),]
    root <- root[-which(root$plantid==472 &root$taskid==19291),]
    
  }
  
  ### removing empty columns
  root <- root[,colSums(is.na(root))<nrow(root)]
  
  ### removing useless columns
  root <- root[,-c(16:25,27:36)]
  
  ### adding a "date" column
  root$date <- rep(which(dates==d),nrow(root))
  
  ### merging with harvest data to have access to genotype id, etc
  colnames(root)[which(colnames(root)=="plantid")] <- "RT_ID"
  colnames(root)[which(colnames(root)=="ABCDEF")] <- "Plant_ID"
  
  root <- merge(root, har[,c("RT_ID", "Plant_ID", "Focal", "Neighbour","stand","Treatment","Block")], by=c("RT_ID","Plant_ID"))
  
  ### reordering columns
  root <- root[,c(17,1,2,18:22,3:16)]
  
  all_root <- rbind(all_root,root)
}


## Counting the number of data points per rhizotube per date
RT_count <- data.frame()

for (d in c(1:3)) {
  root_tmp <- all_root[which(all_root$date==d),]
  RTS <- unique(root_tmp$RT_ID)
  
  for (r in RTS) {
    
    root_RT <- root_tmp[which(root_tmp$RT_ID==r),]
    n <- nrow(unique(root_RT[,c(6:22)])) 
    df <-
      data.frame(
        Date = d,
        RT_ID = r,
        Nb_measures = n
      )
    RT_count <- rbind(RT_count, df)
  }
}

table(RT_count$Date, RT_count$Nb_measures)

## ---->>>> Most rhizotubes have at least two plants merged together, which means we can only analyze root data at the rhyzotube level

### Computing root trait data
all_root_trait <- data.frame()
for (d in unique(all_root$date)) {
  for(r in unique(all_root$RT_ID)) {
    
    root_tmp <- all_root[which(all_root$date==d & all_root$RT_ID==r),]
    root_trait <- root_tmp[1,c(1:2,4:8)]
    root_trait$BE_BoitEng_Hauteur_mm <- max(root_tmp$BE_BoitEng_Hauteur_mm)
    root_trait$SPUR_Squelette_Corrige_mm <- sum(unique(root_tmp$SPUR_Squelette_Corrige_mm))
    root_trait$SURF_Surface_Projetee_mm2 <- sum(unique(root_tmp$SURF_Surface_Projetee_mm2))
    
    all_root_trait <- rbind(all_root_trait,root_trait)
  }
}

## data formatting
all_root_trait$pair <- as.factor(paste(all_root_trait$Focal, all_root_trait$Neighbour, sep=";"))
all_root_trait$Treatment <- as.factor(all_root_trait$Treatment)
all_root_trait$Treatment  = factor(all_root_trait$Treatment ,levels(all_root_trait$Treatment )[c(2,1)])
all_root_trait$stand <- as.factor(all_root_trait$stand)
all_root_trait$stand = factor(all_root_trait$stand,levels(all_root_trait$stand)[c(2,1)])
all_root_trait$Block <- as.factor(all_root_trait$Block)
all_root_trait$date <- as.factor(all_root_trait$date)

### Removing some strange values
all_root_trait[which(all_root_trait$SPUR_Squelette_Corrige_mm>100000),"SPUR_Squelette_Corrige_mm"] <- NA

### Overall stats
mod <- lmer(BE_BoitEng_Hauteur_mm ~ date + Treatment + date:Treatment + Treatment/Block + stand + Treatment:stand + date:stand + date:Treatment:stand + (1|pair) , data=all_root_trait)
summary(mod)
#anova(mod, ddf = "Kenward-Roger") 

mod <- lmer(SPUR_Squelette_Corrige_mm ~ date + Treatment  + date:Treatment +  Treatment/Block + stand + Treatment:stand + date:stand + date:Treatment:stand + (1|pair), data=all_root_trait)
summary(mod)
#anova(mod, ddf = "Kenward-Roger") 

mod <- lmer(SURF_Surface_Projetee_mm2 ~ date + Treatment +  date:Treatment + Treatment/Block + stand + Treatment:stand + date:stand + (1|pair), data=all_root_trait)
summary(mod)
#anova(mod, ddf = "Kenward-Roger") 

## Subsetting data per date
all_root_trait_d1 <- all_root_trait[which(all_root_trait$date=="1"),]
all_root_trait_d2 <- all_root_trait[which(all_root_trait$date=="2"),]
all_root_trait_d3 <- all_root_trait[which(all_root_trait$date=="3"),]

## Outputting final files
write.csv(all_root_trait_d1,file="data/processed_data/Root_traits_D1.csv",row.names=F)
write.csv(all_root_trait_d2,file="data/processed_data/Root_traits_D2.csv",row.names=F)
write.csv(all_root_trait_d3,file="data/processed_data/Root_traits_D3.csv",row.names=F)
