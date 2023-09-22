#### G. Montazeaud
#### 12.09.2023

### This script aim to process raw spectrometric data that was measured on wheat leaves as part of the SOlACE experiment. Spectrum where used to measure leaf N content (%). 
### Spectrum where measured over 4 days (15, 16, 17, and 18 July 2019), we have one seperate file per day of measurement

### Loading spectrometric information
info_J1 <- read.csv("data/raw_data/spectrometric_data/NIRS_spectrum_to_RT_number_J1.csv", header=T, sep=",")
info_J2 <- read.csv("data/raw_data/spectrometric_data/NIRS_spectrum_to_RT_number_J2.csv", header=T, sep=",")
info_J3 <- read.csv("data/raw_data/spectrometric_data/NIRS_spectrum_to_RT_number_J3.csv", header=T, sep=",")
info_J4 <-read.csv("data/raw_data/spectrometric_data/NIRS_spectrum_to_RT_number_J4.csv", header=T, sep=",")

### Loading leafN data
data_J1 <- read.csv("data/raw_data/spectrometric_data/NIRS_Nvalues_J1.csv", header=T, sep=",")
data_J2 <- read.csv("data/raw_data/spectrometric_data/NIRS_Nvalues_J2.csv", header=T, sep=",")
data_J3 <- read.csv("data/raw_data/spectrometric_data/NIRS_Nvalues_J3.csv", header=T, sep=",")
data_J4 <- read.csv("data/raw_data/spectrometric_data/NIRS_Nvalues_J4.csv", header=T, sep=",")


##################
##### Day 1 ######
##################


### Merging spectrum information with leafN data
J1 <- merge(data_J1, info_J1, by= "spectre")

### Removing problematic rhyzotubes (some spectrum are problematic, or the number of measurements per rhyzotube does not match the expected number of plants (6))
J1 <- J1[-which(J1$comment=="PB"),]
PB <- names(table(J1$rhyzotube)[which(!table(J1$rhyzotube)%in%c(6,12))])
J1 <- J1[-which(J1$rhyzotube%in%PB),]
unique(table(J1$rhyzotube))

## Adding the plant_ID
rhyzo_12 <- names(table(J1$rhyzotube)[which(table(J1$rhyzotube)==12)])
J1_rhyzo_12 <- J1[which(J1$rhyzotube%in%rhyzo_12),]
J1_rhyzo_12 <- J1_rhyzo_12[order(J1_rhyzo_12$rhyzotube,J1_rhyzo_12$spectre),]
J1_rhyzo_12$Plant_ID <- rep(rep(c("A","B","C","D","E","F"), each=2),length(unique(J1_rhyzo_12$rhyzotube)))
J1_rhyzo_12 <- aggregate(azote~rhyzotube+Plant_ID, FUN="mean", data=J1_rhyzo_12)

rhyzo_6 <- names(table(J1$rhyzotube)[which(table(J1$rhyzotube)==6)])
J1_rhyzo_6 <- J1[which(J1$rhyzotube%in%rhyzo_6),]
J1_rhyzo_6 <- J1_rhyzo_6[order(J1_rhyzo_6$rhyzotube,J1_rhyzo_6$spectre),]
J1_rhyzo_6$Plant_ID <- rep(c("A","B","C","D","E","F"),length(unique(J1_rhyzo_6$rhyzotube)))
J1_rhyzo_6 <- J1_rhyzo_6[,c(3,5,2)]

J1 <- rbind(J1_rhyzo_12,J1_rhyzo_6)
J1 <- J1[order(J1$rhyzotube,J1$Plant_ID),]
J1$Day_of_leafNmeasurement <- 1

##################
##### Day 2 ######
##################

J2 <- merge(data_J2, info_J2, by= "spectre")

### Removing problematic rhyzotubes (some spectrum are problematic, or the number of measurements per rhyzotube does not match the expected number of plants (6))
J2 <- J2[-which(J2$comment=="PB"),]
PB <- names(table(J2$rhyzotube)[which(!table(J2$rhyzotube)%in%c(6,12))])
J2 <- J2[-which(J2$rhyzotube%in%PB),]
unique(table(J2$rhyzotube))

## Adding the plant_ID
rhyzo_6 <- names(table(J2$rhyzotube)[which(table(J2$rhyzotube)==6)])
J2 <- J2[which(J2$rhyzotube%in%rhyzo_6),]
J2 <- J2[order(J2$rhyzotube,J2$spectre),]
J2$Plant_ID <- rep(c("A","B","C","D","E","F"),length(unique(J2$rhyzotube)))
J2 <- J2[,c(3,5,2)]

J2 <- J2[order(J2$rhyzotube,J2$Plant_ID),]
J2$Day_of_leafNmeasurement <- 2



##################
##### Day 3 ######
##################

J3 <- merge(data_J3, info_J3, by= "spectre")

### Removing problematic rhyzotubes (some spectrum are problematic, or the number of measurements per rhyzotube does not match the expected number of plants (6))
J3 <- J3[-which(J3$comment=="PB"),]
PB <- names(table(J3$rhyzotube)[which(!table(J3$rhyzotube)%in%c(6,12))])
J3 <- J3[-which(J3$rhyzotube%in%PB),]
unique(table(J3$rhyzotube))

## Adding the plant_ID
rhyzo_6 <- names(table(J3$rhyzotube)[which(table(J3$rhyzotube)==6)])
J3 <- J3[which(J3$rhyzotube%in%rhyzo_6),]
J3 <- J3[order(J3$rhyzotube,J3$spectre),]
J3$Plant_ID <- rep(c("A","B","C","D","E","F"),length(unique(J3$rhyzotube)))
J3 <- J3[,c(3,5,2)]

J3 <- J3[order(J3$rhyzotube,J3$Plant_ID),]
J3$Day_of_leafNmeasurement <- 3

##################
##### Day 4 ######
##################

J4 <- merge(data_J4, info_J4, by= "spectre")

### Removing problematic rhyzotubes (some spectrum are problematic, or the number of measurements per rhyzotube does not match the expected number of plants (6))
J4 <- J4[-which(J4$comment=="PB"),]
PB <- names(table(J4$rhyzotube)[which(!table(J4$rhyzotube)%in%c(6,12))])
J4 <- J4[-which(J4$rhyzotube%in%PB),]
unique(table(J4$rhyzotube))

## Adding the plant_ID
rhyzo_6 <- names(table(J4$rhyzotube)[which(table(J4$rhyzotube)==6)])
J4 <- J4[which(J4$rhyzotube%in%rhyzo_6),]
J4 <- J4[order(J4$rhyzotube,J4$spectre),]
J4$Plant_ID <- rep(c("A","B","C","D","E","F"),length(unique(J4$rhyzotube)))
J4 <- J4[,c(3,5,2)]

J4 <- J4[order(J4$rhyzotube,J4$Plant_ID),]
J4$Day_of_leafNmeasurement <- 4

########################
##### Combining all data
########################

leafN <- do.call(rbind,list(J1,J2,J3,J4))
colnames(leafN) <- c("RT_ID","Plant_ID","leaf_N","Day_leafN_measurement")

##############
## Data quality control ##
##############

### Checking if some rhyzotubes have been measured several times
leafN$RT_ID <- as.factor(leafN$RT_ID)
count <- as.data.frame.matrix(table(leafN$RT_ID, leafN$Day_leafN_measurement))
names(which(rowSums(count)>6))
### Removing data for RT 313 on day 1 (was also measured on day 4)
leafN <- leafN[-which(leafN$RT_ID==313 & leafN$Day_leafN_measurement==1),]


### Checking leafN distribution
hist(leafN$leaf_N)

### Removing extremal outlier values
leafN <- leafN[-which(leafN$leaf_N<2 | leafN$leaf_N>5),]
hist(leafN$leaf_N)

### outputting leafN data
write.csv(leafN,"data/processed_data/leafN_processed.csv", row.names = F)

