### Script to analyze the water status and nutrient solution provisionning of the Rhizotubes
library(ggplot2)
library(ggpubr)
Sys.setlocale("LC_ALL", "English")


### I. Daily provision of nutrient solution
nutri <- read.csv("./data/raw_data/greenhouse_data/Daily_nutrient_solution_inputs.csv", header=T, sep=";", dec=",")
nutriC <- nutri[,c("date","C")]
colnames(nutriC)[2] <- "Dose"
nutriS <- nutri[c("date","S")]
colnames(nutriS)[2] <- "Dose"
nutri <- rbind(nutriC,nutriS)
nutri$Treatment <- c(rep("C",nrow(nutriC)), rep("S", nrow(nutriS)))
summary(nutri)
date_format <- "%d/%m/%Y %H:%M"
nutri$date <- as.POSIXct(nutri$date,format = date_format, tz = "GMT")
nutri$Treatment <- as.factor(nutri$Treatment)

p1 <- ggplot(data=nutri[which(nutri$date>=as.POSIXct("23/06/2019 9:00",format = date_format, tz = "GMT")),], aes(x=date, y=Dose, fill=Treatment))+
  geom_bar(position = position_dodge(), stat="identity", width=25000)+
  scale_fill_manual(values=c("blue","red"))+
  labs(y="Dose (mL)",x="")+
  annotate(geom="text",
           x=as.POSIXct("24/06/2019 7:00",format = date_format, tz = "GMT"),
           y=1100,
           label="Seedling\ntransfer",
           color="black",
           size=3,
           vjust=1)+
  annotate(geom="text",
           x=as.POSIXct("15/07/2019 5:00",format = date_format, tz = "GMT"),
           y=1100,
           label="Harvest",
           color="black",
           size=3)+
  theme_bw()


### II. Dynamic water status of the Rhizotubes (average over all Rhizotubes)
water <- read.table("./data/raw_data/greenhouse_data/Average_water_storage.csv", header=T, sep=";", dec=",")
waterC <- water[,c("date","C")]
colnames(waterC)[2] <- "water_status"
waterS <- water[c("date","S")]
colnames(waterS)[2] <- "water_status"
water <- rbind(waterC,waterS)
water$Treatment <- c(rep("C",nrow(waterC)), rep("S", nrow(waterS)))
summary(water)
date_format <- "%d/%m/%y %H:%M"
water$date <- as.POSIXct(water$date,format = date_format, tz = "GMT")

p2 <- ggplot(aes(x=date, y=water_status, color=Treatment), data=water[which(water$date>=as.POSIXct("23/6/19 9:00",format = date_format, tz = "GMT")),])+
  geom_line()+
  scale_color_manual(values=c("blue","red"))+
  labs(y="% water storage capacity",x="")+
  annotate(geom="text",
           x=as.POSIXct("24/6/19 7:00",format = date_format, tz = "GMT"),
           y=0.5,
           label="Seedling\ntransfer",
           color="black",
           size=3,
           vjust=-0.2)+
  annotate(geom="text",
           x=as.POSIXct("15/7/19 5:00",format = date_format, tz = "GMT"),
           y=0.5,
           label="Harvest",
           color="black",
           size=3,
           vjust=-0.2)+
  theme_bw()


ggarrange(plotlist=list(p1,p2),
          nrow=2,
          ncol=1)
ggsave("./manuscript/Figures_&_tables/FigS1_rhizotube_monitoring.png", dpi=300, height=5, width=6)
