### G. Montazeaud

### Analysis of wheat mixture experimental data. 36 genotypes were grown in monocultures and in 54 two-way mixtures with two treatments (C = non limiting water and nutrients, S = limiting water and nutrients) in controled conditions (Dijon, 4PMI platform)

## Mixing varieties reduces root competition between wheat seedlings grown under resource-limited conditions
######################################
### PACKAGES ###
######################################
library(lme4)
library(arm)
library(ggplot2)
library(cowplot)
library(agricolae)
library(lmerTest)
library(glmulti)
library(ggpubr)

####################################
### COLOURS ########################
cols_Treatment <- c("blue","red")
cols_stand <- c("grey","bisque2")

###########################################
#### I. ANALYZING DATA MEASURED AT THE INDIVIDUAL PLANT LEVEL
###########################################

### we here only consider harvest trait data, as root trait data is not available at the single plant level

### Loading harvest data
har <- read.table("data/raw_data/biomass_data/harvest_data.csv", header=T, sep=";")

### Merging harvest data with leaf Nitrogen data (measured with Near-Infrared Spectroscopy)
leafN <- read.csv("data/processed_data/leafN_processed.csv", header=T)
har <- merge(har,leafN, by=c("RT_ID","Plant_ID"),all.x=T)

## Re-ordering factor levels
har$Treatment <- as.factor(har$Treatment)
har$stand <- as.factor(har$stand)
har$stand <- factor(har$stand,levels(har$stand)[c(2,1)])
har$Block <- as.factor(har$Block)

## Adding +1 tiller to all plants so that the baseline tiller number is 1 and not 0
har$Tiller_nb = har$Tiller_nb+1

### Data distribution
missing <- format(round(sum(is.na(har$Leaf_nb_main_stem))/nrow(har)*100,2), nsmall=2)
plt_leaf <- ggplot(har, aes(x=Leaf_nb_main_stem)) +
  geom_histogram(binwidth=1, colour="black", fill="grey")+
  labs(x="Nb leaves on the main stem",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(har)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(2000,0.95*2000),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(har$Tiller_nb))/nrow(har)*100,2), nsmall=2)
plt_tiller <- ggplot(har, aes(x=Tiller_nb)) +
  geom_histogram(binwidth=1, colour="black", fill="grey")+
  labs(x="Nb of tillers",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(har)), paste("%NA = ", missing, sep="")),
           x=1,
           y=c(1900,0.95*1900),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(har$leaf_N))/nrow(har)*100,2), nsmall=2)
plt_leafN <- ggplot(har, aes(x=leaf_N)) +
  geom_histogram(binwidth=0.1, colour="black", fill="grey")+
  labs(x="Leaf nitrogen content (%)",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(har)), paste("%NA = ", missing, sep="")),
           x=2,
           y=c(500,0.95*500),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(har$Root_Shoot_Ratio))/nrow(har)*100,2), nsmall=2)
plt_RSratio <- ggplot(har, aes(x=Root_Shoot_Ratio)) +
  geom_histogram(binwidth=0.1, colour="black", fill="grey")+
  labs(x="Root:Shoot ratio",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(har)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(800,0.95*800),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(har$Shoot_DW))/nrow(har)*100,2), nsmall=2)
plt_Shoot_DW <- ggplot(har, aes(x=Shoot_DW)) +
  geom_histogram(binwidth=25, colour="black", fill="grey")+
  labs(x="Shoot biomass (mg)",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(har)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(500,0.95*500),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(har$Root_DW))/nrow(har)*100,2), nsmall=2)
plt_Root_DW <- ggplot(har, aes(x=Root_DW)) +
  geom_histogram(binwidth=25, colour="black", fill="grey")+
  labs(x="Root biomass (mg)",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(har)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(800,0.95*800),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(har$Total_DW))/nrow(har)*100,2), nsmall=2)
plt_Total_DW <- ggplot(har, aes(x=Total_DW)) +
  geom_histogram(binwidth=50, colour="black", fill="grey")+
  labs(x="Total biomass (mg)",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(har)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(700,0.95*700),
           hjust=0)+
  theme_bw()

plot_grid(plotlist = list(plt_leaf, plt_tiller, plt_leafN, plt_RSratio, plt_Shoot_DW, plt_Root_DW, plt_Total_DW),
          nrow=2,
          ncol = 4)
ggsave("outputs/plots/single_plant_trait_distribution.pdf", dpi=300, height=8, width=10)


## Checking 2 plants with very high R:S ratios
har[which(har$Root_Shoot_Ratio>2),]

## Checking data quality (comments on plant damages, missing plants, etc)
summary(har)
har$Comments <- as.factor(har$Comments)
levels(har$Comments)
har_pb <- har[grep("^P",har$Comments),]
har_no_pb <- har[!grepl("^P",har$Comments),]
har_pb_unref <- har_no_pb[!complete.cases(har_no_pb$Leaf_nb_main_stem, har_no_pb$Tiller_nb, har_no_pb$Root_DW, har_no_pb$Shoot_DW, har_no_pb$Total_DW, har_no_pb$Root_Shoot_Ratio),]
pb <- rbind(har_pb, har_pb_unref)
length(unique(pb$RT_ID))

## Checking pairwise correlations between variables

## Defining a graphical function to add correlation coefficients (Pearson'r) and correlation significance to classical pairwise scatter plots. 
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- cor(x, y, use="pairwise")
  txt <- format(round(r, 2), nsmall = 2)
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 1/(3*strwidth(txt)) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex) 
  text(.7, .8, Signif, cex=cex, col=2) 
}


pdf("outputs/plots/pairwise_correlations.pdf", pointsize = 2)
par(cex.axis=1.5, las=1, pch=19)
pairs(har[,c("Leaf_nb_main_stem","Tiller_nb","leaf_N", "Root_Shoot_Ratio","Shoot_DW","Root_DW","Total_DW")], upper.panel = panel.cor, lower.panel = panel.smooth, labels=c("Nb  leaves\non the main stem", "Nb  tillers", "Leaf N\ncontent (%)", "R:S Ratio","Shoot\nbiomass (mg)","Root\nbiomass (mg)","Total\nbiomass (mg)"), cex.labels=1.9)
dev.off()


###########################################
#### II. ANALYSIS AT THE RHIZOTUBE LEVEL
###########################################

################################
### INTEGRATION OF ROOT TRAIT DATA ###
################################

all_root_trait <- read.csv("data/processed_data/Root_traits_D3.csv", header=T)

########## MERGING ROOT TRAIT DATA AND ABOVEGROUND DATA

## Summing all aboveground trait (but leafN) per RT and per genotype
har_ag <- aggregate(.~RT_ID+stand+Treatment+Block+Sampling_ID, data=har[,-which(colnames(har)%in%c("Plant_ID","Root_Shoot_Ratio","Comments","Focal","Neighbour","leaf_N","Day_leafN_measurement"))], FUN=sum)
har_ag$Root_Shoot_Ratio <- har_ag$Root_DW/har_ag$Shoot_DW

## Averaging LeafN per rhyzotube
leafN_ag <- aggregate(leaf_N~RT_ID+stand+Treatment+Block+Day_leafN_measurement, data=har[,-which(colnames(har)%in%c("Plant_ID","Root_Shoot_Ratio","Comments","Focal","Neighbour"))], FUN=mean)

## Merging leafN data and other aboveground traits
har_ag <- merge(har_ag,leafN_ag,by=c("RT_ID", "stand", "Treatment", "Block"), all.x=T)

## merging above and belowground data
all_trait <- merge(har_ag, all_root_trait, by=c("RT_ID", "stand", "Treatment", "Block"))

## Reordering columns
all_trait <- all_trait[,c(20,1:13,17:19)]

## Plotting pot-level trait distribution
missing <- format(round(sum(is.na(all_trait$Leaf_nb_main_stem))/nrow(all_trait)*100,2), nsmall=2)
plt_leaf <- ggplot(all_trait, aes(x=Leaf_nb_main_stem)) +
  geom_histogram(binwidth=1, colour="black", fill="grey")+
  labs(x="Nb leaves on the main stem",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(all_trait)), paste("%NA = ", missing, sep="")),
           x=10,
           y=c(125,0.9*125),
           hjust=0)+
  theme_bw()


missing <- format(round(sum(is.na(all_trait$Tiller_nb))/nrow(all_trait)*100,2), nsmall=2)
plt_tiller <- ggplot(all_trait, aes(x=Tiller_nb)) +
  geom_histogram(binwidth=1, colour="black", fill="grey")+
  labs(x="Nb of tillers",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(all_trait)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(100,0.9*100),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(all_trait$leaf_N))/nrow(all_trait)*100,2), nsmall=2)
plt_leafN <- ggplot(all_trait, aes(x=leaf_N)) +
  geom_histogram(binwidth=0.1, colour="black", fill="grey")+
  labs(x="Leaf N (%)",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(all_trait)), paste("%NA = ", missing, sep="")),
           x=2,
           y=c(110,0.9*110),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(all_trait$Root_Shoot_Ratio))/nrow(all_trait)*100,2), nsmall=2)
plt_RSratio <- ggplot(all_trait, aes(x=Root_Shoot_Ratio)) +
  geom_histogram(binwidth=0.1, colour="black", fill="grey")+
  labs(x="Root:Shoot ratio",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(all_trait)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(190,0.9*190),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(all_trait$Shoot_DW))/nrow(all_trait)*100,2), nsmall=2)
plt_Shoot_DW <- ggplot(all_trait, aes(x=Shoot_DW)) +
  geom_histogram(binwidth=100, colour="black", fill="grey")+
  labs(x="Shoot biomass (mg)",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(all_trait)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(90,0.9*90),
           hjust=0)+
  theme_bw()


missing <- format(round(sum(is.na(all_trait$Root_DW))/nrow(all_trait)*100,2), nsmall=2)
plt_Root_DW <- ggplot(all_trait, aes(x=Root_DW)) +
  geom_histogram(binwidth=100, colour="black", fill="grey")+
  labs(x="Root biomass (mg)",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(all_trait)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(140,0.9*140),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(all_trait$Total_DW))/nrow(all_trait)*100,2), nsmall=2)
plt_Total_DW <- ggplot(all_trait, aes(x=Total_DW)) +
  geom_histogram(binwidth=200, colour="black", fill="grey")+
  labs(x="Total biomass (mg)",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(all_trait)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(110,0.9*110),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(all_trait$BE_BoitEng_Hauteur_mm))/nrow(all_trait)*100,2), nsmall=2)
plt_BE_BoitEng_Hauteur_mm <- ggplot(all_trait, aes(x=BE_BoitEng_Hauteur_mm)) +
  geom_histogram(binwidth=50, colour="black", fill="grey")+
  labs(x="Root depth (mm)",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(all_trait)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(500,0.9*500),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(all_trait$SPUR_Squelette_Corrige_mm))/nrow(all_trait)*100,2), nsmall=2)
plt_SPUR_Squelette_Corrige_mm <- ggplot(all_trait, aes(x=SPUR_Squelette_Corrige_mm)) +
  geom_histogram(binwidth=5000, colour="black", fill="grey")+
  labs(x="Root length (mm)",
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(all_trait)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(200,0.9*200),
           hjust=0)+
  theme_bw()

missing <- format(round(sum(is.na(all_trait$SURF_Surface_Projetee_mm2))/nrow(all_trait)*100,2), nsmall=2)
plt_SURF_Surface_Projetee_mm2 <- ggplot(all_trait, aes(x=SURF_Surface_Projetee_mm2)) +
  geom_histogram(binwidth=5000, colour="black", fill="grey")+
  labs(x=expression("Root Surface (mm"^2*")"),
       y="Nb obs")+
  annotate("text",
           label=c(paste("N obs = ",nrow(all_trait)), paste("%NA = ", missing, sep="")),
           x=0,
           y=c(200,0.9*200),
           hjust=0)+
  theme_bw()


plot_grid(plotlist = list(plt_leaf, plt_tiller, plt_leafN, plt_RSratio, plt_Shoot_DW, plt_Root_DW, plt_Total_DW, plt_BE_BoitEng_Hauteur_mm, plt_SPUR_Squelette_Corrige_mm, plt_SURF_Surface_Projetee_mm2),
          nrow=3,
          ncol = 4)
ggsave("outputs/plots/rhyzotube_trait_distribution.pdf", dpi=300, height=8, width=10)

## Statistical analysis
## We use mixed models to compare monocultures vs mixtures
# Response variable to be changed
mod_leaf_nb <- lmer(Leaf_nb_main_stem ~ Sampling_ID + Block + Treatment + stand + Treatment:stand + (1+Treatment|pair) , data=all_trait)
summary(mod_leaf_nb)
anova(mod_leaf_nb, ddf = "Kenward-Roger") 
write.csv(as.data.frame(anova(mod_leaf_nb, ddf = "Kenward-Roger")),"./outputs/data/leaf_nb_anova.csv", row.names=T)
## --> Significant effect of the treatment, the sampling date, and the experimental block 

mod_till_nb <- lmer(Tiller_nb ~ Sampling_ID + Block + Treatment + stand + Treatment:stand + (1+Treatment|pair) , data=all_trait)
summary(mod_till_nb)
anova(mod_till_nb, ddf = "Kenward-Roger")
write.csv(as.data.frame(anova(mod_till_nb, ddf = "Kenward-Roger")),"./outputs/data/till_nb_anova.csv", row.names=T)
## --> Significant effect of the treatment, the sampling date, and the experimental block 

mod_leafN <- lmer(leaf_N ~ Sampling_ID + Block + Treatment + stand + Treatment:stand + (1+Treatment|pair) , data=all_trait)
summary(mod_leafN)
anova(mod_leafN, ddf = "Kenward-Roger") 
write.csv(as.data.frame(anova(mod_leafN, ddf = "Kenward-Roger")),"./outputs/data/leafN_anova.csv", row.names=T)
## --> Significant effect of the treatment, and the experimental block 

mod_RSratio <- lmer(Root_Shoot_Ratio ~ Sampling_ID + Block + Treatment + stand + Treatment:stand + (1+Treatment|pair), data=all_trait)
summary(mod_RSratio)
anova(mod_RSratio, ddf = "Kenward-Roger") 
write.csv(as.data.frame(anova(mod_RSratio, ddf = "Kenward-Roger")),"./outputs/data/RS_ratio_anova.csv", row.names=T)
## --> Significant effect of the treatment, the sampling date, and the experimental block 

mod_Shoot_DW <- lmer(Shoot_DW ~ Sampling_ID + Block + Treatment + stand + Treatment:stand + (1+Treatment|pair), data=all_trait)
summary(mod_RSratio)
anova(mod_RSratio, ddf = "Kenward-Roger") 
write.csv(as.data.frame(anova(mod_Shoot_DW, ddf = "Kenward-Roger")),"./outputs/data/Shoot_DW_anova.csv", row.names=T)
## --> Significant effect of the treatment, the sampling date, and the experimental block 

mod_Root_DW <- lmer(Root_DW ~ Sampling_ID + Block + Treatment + stand + Treatment:stand + (1+Treatment|pair), data=all_trait)
summary(mod_Root_DW)
anova(mod_Root_DW, ddf = "Kenward-Roger")
write.csv(as.data.frame(anova(mod_Root_DW, ddf = "Kenward-Roger")),"./outputs/data/Root_DW_anova.csv", row.names=T)
## --> Significant effect of the treatment, the sampling date, and the experimental block 

mod_Total_DW <- lmer(Total_DW ~ Sampling_ID + Block + Treatment + stand + Treatment:stand + (1+Treatment|pair), data=all_trait)
summary(mod_Total_DW)
anova(mod_Total_DW, ddf = "Kenward-Roger") 
write.csv(as.data.frame(anova(mod_Total_DW, ddf = "Kenward-Roger")),"./outputs/data/Total_DW_anova.csv", row.names=T)
## --> Significant effect of the treatment, the sampling date, and the experimental block 

mod_BE_BoitEng_Hauteur_mm <- lmer(BE_BoitEng_Hauteur_mm ~ Block + Treatment + stand + Treatment:stand + (1+Treatment|pair) , data=all_trait)
summary(mod_BE_BoitEng_Hauteur_mm)
anova(mod_BE_BoitEng_Hauteur_mm, ddf = "Kenward-Roger") 
write.csv(as.data.frame(anova(mod_BE_BoitEng_Hauteur_mm, ddf = "Kenward-Roger")),"./outputs/data/prof_rac_anova.csv", row.names=T)
## --> Significant effect of the treatment only

mod_SPUR_Squelette_Corrige_mm <- lmer(SPUR_Squelette_Corrige_mm ~ Block + Treatment + stand + Treatment:stand + (1+Treatment|pair) , data=all_trait)
summary(mod_SPUR_Squelette_Corrige_mm)
anova(mod_SPUR_Squelette_Corrige_mm, ddf = "Kenward-Roger") 
write.csv(as.data.frame(anova(mod_SPUR_Squelette_Corrige_mm, ddf = "Kenward-Roger")),"./outputs/data/long_rac_anova.csv", row.names=T)
## --> Significant effect of the treatment, and the experimental block 

mod_SURF_Surface_Projetee_mm2 <- lmer(SURF_Surface_Projetee_mm2 ~ Block + Treatment + stand + Treatment:stand + (1+Treatment|pair), data=all_trait)
summary(mod_SURF_Surface_Projetee_mm2)
as.data.frame(anova(mod_SURF_Surface_Projetee_mm2, ddf = "Kenward-Roger") )
write.csv(as.data.frame(anova(mod_SURF_Surface_Projetee_mm2, ddf = "Kenward-Roger")),"./outputs/data/surf_rac_anova.csv", row.names=T)
## --> Significant effect of the treatment, and the experimental block 


## Plotting the effect of the treatment on above an belowground traits
plt_leaf <- ggplot(data = all_trait, aes(x=Treatment, y=Leaf_nb_main_stem, color=Treatment)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none") +
  annotate("text",
           x=c(1:2),
           y=13,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$Treatment=="C"),"Leaf_nb_main_stem"])), length(na.omit(all_trait[which(all_trait$Treatment=="S"),"Leaf_nb_main_stem"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=44,
           vjust=1,
           label="***",
           size=6)+
  ylim(12,45)+
  labs(y = "# leaves",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.1)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_tiller <-  ggplot(data = all_trait, aes(x=Treatment, y=Tiller_nb, color=Treatment)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none") +
  ylim(2,35)+
  annotate("text",
           x=c(1:2),
           y=3,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$Treatment=="C"),"Tiller_nb"])), length(na.omit(all_trait[which(all_trait$Treatment=="S"),"Tiller_nb"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=34,
           vjust=1,
           label="***",
           size=6)+
  labs(y = "# tillers",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.1)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_leafN <- ggplot(data = all_trait, aes(x=Treatment, y=leaf_N, color=Treatment)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none") +
  ylim(2.3,4)+
  annotate("text",
           x=c(1:2),
           y=2.35,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$Treatment=="C"),"leaf_N"])), length(na.omit(all_trait[which(all_trait$Treatment=="S"),"leaf_N"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=3.9,
           vjust=1,
           label="***",
           size=6)+
  labs(y = "Leaf N (%)",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.1)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_Shoot_DW <- ggplot(data = all_trait, aes(x=Treatment, y=Shoot_DW, color=Treatment)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none") +
  ylim(0,3300)+
  annotate("text",
           x=c(1:2),
           y=150,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$Treatment=="C"),"Shoot_DW"])), length(na.omit(all_trait[which(all_trait$Treatment=="S"),"Shoot_DW"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=3250,
           vjust=1,
           label="***",
           size=6)+
  labs(y = "Shoot biomass (mg)",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.1)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_Root_DW <- ggplot(data = all_trait, aes(x=Treatment, y=Root_DW, color=Treatment)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none") +
  ylim(100,1600)+
  annotate("text",
           x=c(1:2),
           y=150,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$Treatment=="C"),"Root_DW"])), length(na.omit(all_trait[which(all_trait$Treatment=="S"),"Root_DW"]))),
           size=3.5)+
  annotate("text",
             x=1.5,
             y=1550,
             vjust=1,
             label="***",
             size=6)+
  labs(y = "Root biomass (mg)",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.1)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))


plt_Total_DW <- ggplot(data = all_trait, aes(x=Treatment, y=Total_DW, color=Treatment)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none") +
  ylim(250,4800) +
  annotate("text",
           x=c(1:2),
           y=350,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$Treatment=="C"),"Total_DW"])), length(na.omit(all_trait[which(all_trait$Treatment=="S"),"Total_DW"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=4700,
           vjust=1,
           label="***",
           size=6)+
  labs(y = "Total biomass (mg)",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.1)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_RSratio <- ggplot(data = all_trait, aes(x=Treatment, y=Root_Shoot_Ratio, color=Treatment)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  ylim(0.20,1.25)+
  scale_color_manual(values=cols_Treatment, guide="none") +
  annotate("text",
           x=c(1:2),
           y=0.23,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$Treatment=="C"),"Root_Shoot_Ratio"])), length(na.omit(all_trait[which(all_trait$Treatment=="S"),"Root_Shoot_Ratio"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=1.22,
           vjust=1,
           label="***",
           size=6)+
  labs(y = "Root:Shoot ratio",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.1)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_SPUR_Squelette_Corrige_mm <- ggplot(data = all_trait, aes(x=Treatment, y=SPUR_Squelette_Corrige_mm, color=Treatment)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none") +
  ylim(13000,80000)+
  annotate("text",
           x=c(1:2),
           y=15000,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$Treatment=="C"),"SPUR_Squelette_Corrige_mm"])), length(na.omit(all_trait[which(all_trait$Treatment=="S"),"SPUR_Squelette_Corrige_mm"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=79000,
           vjust=1,
           label="***",
           size=6)+
  labs(y = "Root length (mm)",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.1)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_SURF_Surface_Projetee_mm2 <- ggplot(data = all_trait, aes(x=Treatment, y=SURF_Surface_Projetee_mm2, color=Treatment)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none") +
  ylim(13000,80000)+
  annotate("text",
           x=c(1:2),
           y=15000,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$Treatment=="C"),"SURF_Surface_Projetee_mm2"])), length(na.omit(all_trait[which(all_trait$Treatment=="S"),"SURF_Surface_Projetee_mm2"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=79000,
           vjust=1,
           label="**",
           size=6)+
  labs(y = expression("Root surface (mm"^2*")"),
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.1)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))
  

plot_grid(plotlist = list(plt_leaf, plt_tiller, plt_leafN, plt_Shoot_DW, plt_Root_DW, plt_Total_DW, plt_RSratio, plt_SPUR_Squelette_Corrige_mm, plt_SURF_Surface_Projetee_mm2),
          nrow=3,
          ncol = 3, 
          labels=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"),
          label_fontface = "bold",
          hjust = 0, label_x = 0.01)
ggsave("outputs/plots/pot_level_treatment_effect.pdf", dpi=300, height=8, width=8)

## Plotting the effect of the stand type (pure vs mixed) on above an belowground traits
plt_leaf <- ggplot(data = all_trait, aes(x=stand, y=Leaf_nb_main_stem, color=stand)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_stand, guide="none") +
  annotate("text",
           x=c(1:2),
           y=13,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$stand=="pure"),"Leaf_nb_main_stem"])), length(na.omit(all_trait[which(all_trait$stand=="mixed"),"Leaf_nb_main_stem"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=44,
           vjust=1,
           label="n.s.",
           size=4)+
  ylim(12,45)+
  labs(y = "# leaves on the main stems",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_tiller <-  ggplot(data = all_trait, aes(x=stand, y=Tiller_nb, color=stand)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_stand, guide="none") +
  ylim(-1,35)+
  annotate("text",
           x=c(1:2),
           y=0,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$stand=="pure"),"Tiller_nb"])), length(na.omit(all_trait[which(all_trait$stand=="mixed"),"Tiller_nb"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=34,
           vjust=1,
           label="n.s.",
           size=4)+
  labs(y = "# tillers",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_leafN <- ggplot(data = all_trait, aes(x=stand, y=leaf_N, color=stand)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_stand, guide="none") +
  ylim(2.2,4)+
  annotate("text",
           x=c(1:2),
           y=2.25,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$stand=="pure"),"leaf_N"])), length(na.omit(all_trait[which(all_trait$stand=="mixed"),"leaf_N"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=3.9,
           vjust=1,
           label="n.s.",
           size=4)+
  labs(y = "Leaf N (%)",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_Shoot_DW <- ggplot(data = all_trait, aes(x=stand, y=Shoot_DW, color=stand)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_stand, guide="none") +
  ylim(0,3300)+
  annotate("text",
           x=c(1:2),
           y=125,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$stand=="pure"),"Shoot_DW"])), length(na.omit(all_trait[which(all_trait$stand=="mixed"),"Shoot_DW"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=3250,
           vjust=1,
           label="n.s.",
           size=4)+
  labs(y = "Shoot biomass (mg)",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_Root_DW <- ggplot(data = all_trait, aes(x=stand, y=Root_DW, color=stand)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_stand, guide="none") +
  ylim(90,1600)+
  annotate("text",
           x=c(1:2),
           y=130,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$stand=="pure"),"Root_DW"])), length(na.omit(all_trait[which(all_trait$stand=="mixed"),"Root_DW"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=1550,
           vjust=1,
           label="n.s.",
           size=4)+
  labs(y = "Root biomass (mg)",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))


plt_Total_DW <- ggplot(data = all_trait, aes(x=stand, y=Total_DW, color=stand)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_stand, guide="none") +
  ylim(250,4800) +
  annotate("text",
           x=c(1:2),
           y=350,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$stand=="pure"),"Total_DW"])), length(na.omit(all_trait[which(all_trait$stand=="mixed"),"Total_DW"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=4700,
           vjust=1,
           label="n.s.",
           size=4)+
  labs(y = "Total biomass (mg)",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_RSratio <- ggplot(data = all_trait, aes(x=stand, y=Root_Shoot_Ratio, color=stand)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  ylim(0.10,1.25)+
  scale_color_manual(values=cols_stand, guide="none") +
  annotate("text",
           x=c(1:2),
           y=0.12,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$stand=="pure"),"Root_Shoot_Ratio"])), length(na.omit(all_trait[which(all_trait$stand=="mixed"),"Root_Shoot_Ratio"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=1.22,
           vjust=1,
           label="n.s.",
           size=4)+
  labs(y = "Root:Shoot ratio",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_SPUR_Squelette_Corrige_mm <- ggplot(data = all_trait, aes(x=stand, y=SPUR_Squelette_Corrige_mm, color=stand)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_stand, guide="none") +
  ylim(13000,80000)+
  annotate("text",
           x=c(1:2),
           y=15000,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$stand=="pure"),"SPUR_Squelette_Corrige_mm"])), length(na.omit(all_trait[which(all_trait$stand=="mixed"),"SPUR_Squelette_Corrige_mm"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=79000,
           vjust=1,
           label="n.s.",
           size=4)+
  labs(y = "Root length (mm)",
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_SURF_Surface_Projetee_mm2 <- ggplot(data = all_trait, aes(x=stand, y=SURF_Surface_Projetee_mm2, color=stand)) + 
  geom_violin(trim = F, size=0.9) + 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_stand, guide="none") +
  ylim(13000,80000)+
  annotate("text",
           x=c(1:2),
           y=15000,
           vjust=1,
           label=c(length(na.omit(all_trait[which(all_trait$stand=="pure"),"SURF_Surface_Projetee_mm2"])), length(na.omit(all_trait[which(all_trait$stand=="mixed"),"SURF_Surface_Projetee_mm2"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=79000,
           vjust=1,
           label="n.s.",
           size=4)+
  labs(y = expression("Root surface (mm"^2*")"),
       x = "") +
  theme_classic() +
  theme(axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))


plot_grid(plotlist = list(plt_leaf, plt_tiller, plt_leafN, plt_Shoot_DW, plt_Root_DW, plt_Total_DW, plt_RSratio, plt_SPUR_Squelette_Corrige_mm, plt_SURF_Surface_Projetee_mm2),
          nrow=3,
          ncol = 3)

ggsave("outputs/plots/pot_level_stand_effect.pdf", dpi=300, height=8, width=8)


## Plotting the effect of the block on above an belowground traits
all_trait$Block <- as.factor(all_trait$Block)

plt_leaf <- ggplot(data = all_trait, aes(x=Block, y=Leaf_nb_main_stem)) + 
  geom_boxplot() + 
  labs(y = "# leaves on the main stems",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_tiller <- ggplot(data = all_trait, aes(x=Block, y=Tiller_nb)) + 
  geom_boxplot() + 
  labs(y = "# tillers",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_leafN <- ggplot(data = all_trait, aes(x=Block, y=leaf_N)) + 
  geom_boxplot() + 
  labs(y = "Leaf N (%)",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_RSratio <- ggplot(data = all_trait, aes(x=Block, y=Root_Shoot_Ratio)) + 
  geom_boxplot() + 
  labs(y = "Root:Shoot ratio",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_Shoot_DW <- ggplot(data = all_trait, aes(x=Block, y=Shoot_DW)) + 
  geom_boxplot() + 
  labs(y = "Shoot biomass (mg)",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_Root_DW <- ggplot(data = all_trait, aes(x=Block, y=Root_DW)) + 
  geom_boxplot() + 
  labs(y = "Root biomass (mg)",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_Total_DW <- ggplot(data = all_trait, aes(x=Block, y=Total_DW)) + 
  geom_boxplot() + 
  labs(y = "Total biomass (mg)",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_SPUR_Squelette_Corrige_mm <- ggplot(data = all_trait, aes(x=Block, y=SPUR_Squelette_Corrige_mm)) + 
  geom_boxplot() + 
  labs(y = "Root lenght (mm)",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_SURF_Surface_Projetee_mm2 <- ggplot(data = all_trait, aes(x=Block, y=SURF_Surface_Projetee_mm2)) + 
  geom_boxplot() + 
  labs(y = expression("Root surface (mm"^2*")"),
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plot_grid(plotlist = list(plt_leaf, plt_tiller, plt_leafN, plt_Shoot_DW, plt_Root_DW, plt_Total_DW, plt_RSratio, plt_SPUR_Squelette_Corrige_mm, plt_SURF_Surface_Projetee_mm2),
          nrow=3,
          ncol = 3)
ggsave("outputs/plots/pot_level_block_effect.pdf", dpi=300, height=8, width=8)


## Plotting the effect of the sampling date on aboveground traits
all_trait$Sampling_ID <-as.factor(all_trait$Sampling_ID)
all_trait$Day_leafN_measurement <-as.factor(all_trait$Day_leafN_measurement)

plt_leaf <- ggplot(data = all_trait, aes(x=Sampling_ID, y=Leaf_nb_main_stem)) + 
  geom_boxplot() + 
  labs(y = "Leaf nb on the main stem",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_tiller <- ggplot(data = all_trait, aes(x=Sampling_ID, y=Tiller_nb)) + 
  geom_boxplot() + 
  labs(y = "Tiller nb",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_leafN <- ggplot(data = na.omit(all_trait), aes(x=Day_leafN_measurement, y=leaf_N)) + 
  geom_boxplot() + 
  labs(y = "Leaf N (%)",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_RSratio <- ggplot(data = all_trait, aes(x=Sampling_ID, y=Root_Shoot_Ratio)) + 
  geom_boxplot() + 
  labs(y = "Root:Shoot ratio",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_Shoot_DW <- ggplot(data = all_trait, aes(x=Sampling_ID, y=Shoot_DW)) + 
  geom_boxplot() + 
  labs(y = "Shoot biomass (mg)",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_Root_DW <- ggplot(data = all_trait, aes(x=Sampling_ID, y=Root_DW)) + 
  geom_boxplot() + 
  labs(y = "Root biomass (mg)",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plt_Total_DW <- ggplot(data = all_trait, aes(x=Sampling_ID, y=Total_DW)) + 
  geom_boxplot() + 
  labs(y = "Total biomass (mg)",
       x = "") +
  theme(panel.grid.major = element_blank(), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plot_grid(plotlist = list(plt_leaf, plt_tiller, plt_leafN, plt_Shoot_DW, plt_Root_DW, plt_Total_DW, plt_RSratio),
          nrow=3,
          ncol = 3)
ggsave("outputs/plots/pot_level_sampling_date_effect.pdf", dpi=300, height=8, width=8)


###########################################
#### III. COMPARING MIXTURES AND MONOCULTURES WITH RYT
###########################################

## Summing all variables per RT and per genotype
har_ag2 <- aggregate(.~RT_ID+Focal+Neighbour+stand+Treatment+Block+Sampling_ID, data=har[,-which(colnames(har)%in%c("Plant_ID","Root_Shoot_Ratio","Comments", "leaf_N","Day_leafN_measurement"))], FUN=sum)
har_ag2$Root_Shoot_Ratio <- har_ag2$Root_DW/har_ag2$Shoot_DW

## Separating monoculture and mixture data
har_ag_M <- droplevels(har_ag2[which(har_ag2$stand=="pure"),])
har_ag_P <- droplevels(har_ag2[which(har_ag2$stand=="mixed"),])

## computing monoculture BLUPS
har_ag_M$Focal <- as.factor(har_ag_M$Focal)
blup_monoc <- data.frame(Genotype=levels(har_ag_M$Focal), Treatment=rep(c("C","S"), each=36))

## Shoot_DW
mod <- lmer(Shoot_DW ~  Sampling_ID + Block + Treatment + (1|Focal), data=har_ag_M)
summary(mod)

blup_monoc[which(blup_monoc$Treatment=="C"),"Shoot_DW"] <- fixef(mod)[1]+ranef(mod)$Focal[1]
blup_monoc[which(blup_monoc$Treatment=="S"),"Shoot_DW"] <- fixef(mod)[1]+fixef(mod)[5]+ranef(mod)$Focal[1]

## Root_DW
mod <- lmer(Root_DW ~  Sampling_ID + Block + Treatment + (1|Focal), data=har_ag_M)
summary(mod)

blup_monoc[which(blup_monoc$Treatment=="C"),"Root_DW"] <- fixef(mod)[1]+ranef(mod)$Focal[1]
blup_monoc[which(blup_monoc$Treatment=="S"),"Root_DW"] <- fixef(mod)[1]+fixef(mod)[5]+ranef(mod)$Focal[1]

## Total_DW
mod <- lmer(Total_DW ~  Sampling_ID + Block + Treatment + (1|Focal), data=har_ag_M)
summary(mod)

blup_monoc[which(blup_monoc$Treatment=="C"),"Total_DW"] <- fixef(mod)[1]+ranef(mod)$Focal[1]
blup_monoc[which(blup_monoc$Treatment=="S"),"Total_DW"] <- fixef(mod)[1]+fixef(mod)[5]+ranef(mod)$Focal[1]


## computing mixture BLUPS
har_ag_P$pair <- as.factor(paste(har_ag_P$Focal, har_ag_P$Neighbour, sep=";"))
blup_mixt <- data.frame(Pair=levels(har_ag_P$pair), Treatment=rep(c("C","S"), each=108))

## Shoot_DW
mod <- lmer(Shoot_DW ~ Sampling_ID + Block + Treatment + (1+Treatment|pair), data=har_ag_P)
summary(mod)

blup_mixt[which(blup_mixt$Treatment=="C"),"Shoot_DW"] <- fixef(mod)[1]+ranef(mod)$pair[1]
blup_mixt[which(blup_mixt$Treatment=="S"),"Shoot_DW"] <- fixef(mod)[1]+fixef(mod)[5]+ranef(mod)$pair[1]+ranef(mod)$pair[2]

## Root_DW
mod <- lmer(Root_DW ~ Sampling_ID + Block + Treatment + (1+Treatment|pair), data=har_ag_P)
summary(mod)

blup_mixt[which(blup_mixt$Treatment=="C"),"Root_DW"] <- fixef(mod)[1]+ranef(mod)$pair[1]
blup_mixt[which(blup_mixt$Treatment=="S"),"Root_DW"] <- fixef(mod)[1]+fixef(mod)[5]+ranef(mod)$pair[1]+ranef(mod)$pair[2]

## Total_DW
mod <- lmer(Total_DW ~ Sampling_ID + Block + Treatment + (1+Treatment|pair), data=har_ag_P)
summary(mod)

blup_mixt[which(blup_mixt$Treatment=="C"),"Total_DW"] <- fixef(mod)[1]+ranef(mod)$pair[1]
blup_mixt[which(blup_mixt$Treatment=="S"),"Total_DW"] <- fixef(mod)[1]+fixef(mod)[5]+ranef(mod)$pair[1]+ranef(mod)$pair[2]


## Computing RYs
for (i in c(1:nrow(blup_mixt))) {
  foc <- as.character(strsplit(blup_mixt[i,"Pair"], ";")[[1]][1])
  t <- as.character(blup_mixt[i,"Treatment"]) 
  
  for (v in c("Shoot_DW","Root_DW","Total_DW")) {
    blup_mixt[i,paste(v,"RY",sep="_")] <- blup_mixt[i,v]/blup_monoc[which(blup_monoc$Genotype==foc & blup_monoc$Treatment==t),v]
  }
}


## Computing RYT
for (i in c(1:nrow(blup_mixt))) {
  
  foc_oriented <- blup_mixt[i,"Pair"]
  nei_oriented <- paste(strsplit(blup_mixt[i,"Pair"], ";")[[1]][2],strsplit(blup_mixt[i,"Pair"], ";")[[1]][1], sep=";")
  
  foc <- strsplit(blup_mixt[i,"Pair"], ";")[[1]][1]
  nei <- strsplit(blup_mixt[i,"Pair"], ";")[[1]][2]
  
  t <- as.character(blup_mixt[i,"Treatment"]) 
  
  blup_mixt[i, "Pair_unoriented"] <- paste(sort(c(foc,nei)),collapse=";")
  
  for (v in c("Shoot_DW","Root_DW","Total_DW")) {
    
    blup_mixt[i, paste(v, "RYT",sep="_")] <- blup_mixt[which(blup_mixt$Pair==foc_oriented & blup_mixt$Treatment==t),paste(v,"RY",sep="_")] + blup_mixt[which(blup_mixt$Pair==nei_oriented & blup_mixt$Treatment==t),paste(v,"RY",sep="_")]
    
  }
}

RYT <- unique(blup_mixt[,c("Pair_unoriented","Treatment","Shoot_DW_RYT","Root_DW_RYT","Total_DW_RYT")])

### Data distribution
missing <- format(round(sum(is.na(RYT$Shoot_DW_RYT))/nrow(RYT)*100,2), nsmall=2)
mu <- format(round(mean(RYT$Shoot_DW_RYT, na.rm=T),2), nsmall=2)
plt_RYT_Shoot <- ggplot(RYT, aes(x=Shoot_DW_RYT)) +
  geom_histogram(binwidth=0.05, colour="black", fill="grey")+
  labs(x="",
       y="Nb obs", 
       title="Shoot biomass (RYT)")+
  annotate("text",
           label=c(paste("N obs = ",nrow(RYT)), paste("%NA = ", missing, sep=""), paste("mean = ",mu,"***", sep="")),
           x=0.75,
           y=c(50,0.95*50,0.90*50),
           hjust=0)+
  geom_vline(xintercept = 1, linetype=2)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))


missing <- format(round(sum(is.na(RYT$Root_DW_RYT))/nrow(RYT)*100,2), nsmall=2)
mu <- format(round(mean(RYT$Root_DW_RYT, na.rm=T),2), nsmall=2)
plt_RYT_Root <- ggplot(RYT, aes(x=Root_DW_RYT)) +
  geom_histogram(binwidth=0.05, colour="black", fill="grey")+
  labs(x="",
       y="Nb obs", 
       title="Root biomass (RYT)")+
  annotate("text",
           label=c(paste("N obs = ",nrow(RYT)), paste("%NA = ", missing, sep=""), paste("mean = ",mu,"***", sep="")),
           x=0.75,
           y=c(50,0.95*50,0.90*50),
           hjust=0)+
  geom_vline(xintercept = 1, linetype=2)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

missing <- format(round(sum(is.na(RYT$Total_DW_RYT))/nrow(RYT)*100,2), nsmall=2)
mu <- format(round(mean(RYT$Total_DW_RYT, na.rm=T),2), nsmall=2)
plt_RYT_Total <- ggplot(RYT, aes(x=Total_DW_RYT)) +
  geom_histogram(binwidth=0.05, colour="black", fill="grey")+
  labs(x="",
       y="Nb obs", 
       title="Total biomass (RYT)")+
  annotate("text",
           label=c(paste("N obs = ",nrow(RYT)), paste("%NA = ", missing, sep=""), paste("mean = ",mu,"***", sep="")),
           x=0.75,
           y=c(50,0.95*50,0.90*50),
           hjust=0)+
  geom_vline(xintercept = 1, linetype=2)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

plot_grid(plotlist = list(plt_RYT_Shoot, plt_RYT_Root, plt_RYT_Total),
          nrow=1,
          ncol = 3)
ggsave("outputs/plots/RYT_distribution.pdf", dpi=300, height=4, width=9)

## Statistical analysis
RYT$Treatment <- as.factor(RYT$Treatment)

mod_Shoot_RYT <- lmer(Shoot_DW_RYT ~  Treatment + (1|Pair_unoriented), data=RYT)
summary(mod_Shoot_RYT) 
anova(mod_Shoot_RYT, ddf = "Kenward-Roger")
## --> Significant effect  of the treatment on aboveground biomass RYT

mod_Root_RYT <- lmer(Root_DW_RYT ~  Treatment + (1|Pair_unoriented), data=RYT)
summary(mod_Root_RYT) 
anova(mod_Root_RYT, ddf = "Kenward-Roger")
## --> Significant effect of the treatment on belowground biomass RYT

mod_Total_RYT <- lmer(Total_DW_RYT ~  Treatment + (1|Pair_unoriented), data=RYT)
summary(mod_Total_RYT) 
anova(mod_Total_RYT, ddf = "Kenward-Roger")
## --> Significant effect of the treatment on total biomass RYT


## Plotting the effect of the treatment on RYT
mu <- format(round(mean(RYT$Shoot_DW_RYT, na.rm=T),2), nsmall=2)
plt_RYT_Shoot <- ggplot(data = RYT, aes(x=Treatment, y=Shoot_DW_RYT, color=Treatment)) + 
  geom_violin(trim = F, size=0.9)+ 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none")+
  ylim(0.6,1.3)+
  labs(y = "RYT",
       x = "",
       title = paste0("Shoot biomass ", paste("(mean = ",mu,"**)", sep=""))) +
  annotate("text",
           x=c(1:2),
           y=0.6,
           vjust=1,
           label=c(length(na.omit(RYT[which(RYT$Treatment=="C"),"Shoot_DW_RYT"])), length(na.omit(RYT[which(RYT$Treatment=="S"),"Shoot_DW_RYT"]))),
             size=3.5)+
  annotate("text",
           x=1.5,
           y=1.3,
           vjust=1,
           label="***",
           size=6) + 
  geom_hline(yintercept=1, linetype=2)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))


mu <- format(round(mean(RYT$Root_DW_RYT, na.rm=T),2), nsmall=2)
plt_RYT_Root <- ggplot(data = RYT, aes(x=Treatment, y=Root_DW_RYT, color=Treatment)) + 
  geom_violin(trim = F, size=0.9)+ 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none")+
  ylim(0.6,1.3)+
  labs(y = "RYT",
       x = "",
       title = paste0("Root biomass ", paste("(mean = ",mu,"***)", sep=""))) +
  annotate("text",
           x=c(1:2),
           y=0.6,
           vjust=1,
           label=c(length(na.omit(RYT[which(RYT$Treatment=="C"),"Shoot_DW_RYT"])), length(na.omit(RYT[which(RYT$Treatment=="S"),"Shoot_DW_RYT"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=1.3,
           vjust=1,
           label="***",
           size=6) + 
  geom_hline(yintercept=1, linetype=2)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))


mu <- format(round(mean(RYT$Total_DW_RYT, na.rm=T),2), nsmall=2)
plt_RYT_Total <- ggplot(data = RYT, aes(x=Treatment, y=Total_DW_RYT, color=Treatment)) + 
  geom_violin(trim = F, size=0.9)+ 
  stat_summary(fun.data="mean_sdl",
               fun.args=list(mult=1), 
               geom="pointrange") +
  scale_color_manual(values=cols_Treatment, guide="none")+
  ylim(0.6,1.3)+
  labs(y = "RYT",
       x = "",
       title = paste0("Total biomass ", paste("(mean = ",mu,"***)", sep=""))) +
  annotate("text",
           x=c(1:2),
           y=0.6,
           vjust=1,
           label=c(length(na.omit(RYT[which(RYT$Treatment=="C"),"Shoot_DW_RYT"])), length(na.omit(RYT[which(RYT$Treatment=="S"),"Shoot_DW_RYT"]))),
           size=3.5)+
  annotate("text",
           x=1.5,
           y=1.3,
           vjust=1,
           label="***",
           size=6) + 
  geom_hline(yintercept=1, linetype=2)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x= element_text(size=rel(1.5), margin=margin(t=1)), axis.text.y= element_text(size=rel(1.5)), axis.title.x= element_text(size=rel(1.2)), axis.title.y= element_text(size=rel(1.2)))

plot_grid(plotlist = list(plt_RYT_Shoot, plt_RYT_Root, plt_RYT_Total),
          nrow=1,
          ncol = 3)
ggsave("outputs/plots/pot_level_RYT_treatment_effect.pdf", dpi=300, height=4, width=10)



#######################################
### IV. EXPLAINING MIXTURES'RYT WITH MONOCULTURE TRAITS 
#######################################

## Computing monoculure average trait values in each treatment
trait_monoc <- all_trait[which(all_trait$stand=="pure"),]
for (i in c(1:nrow(trait_monoc))){
  trait_monoc[i,"Focal"] <- strsplit(as.character(trait_monoc[i,"pair"]), ";")[[1]][1]
}
trait_monoc$Focal <- as.factor(trait_monoc$Focal)
blup_trait_monoc <- data.frame(Genotype=levels(trait_monoc$Focal), Treatment=rep(c("C","S"), each=36))

## Leaf Nb
mod <- lmer(Leaf_nb_main_stem ~ Sampling_ID + Block + Treatment + (1|Focal), data=trait_monoc)

blup_trait_monoc[which(blup_trait_monoc$Treatment=="C"),"Leaf_nb_main_stem"] <- fixef(mod)[1]+ranef(mod)$Focal[1]
blup_trait_monoc[which(blup_trait_monoc$Treatment=="S"),"Leaf_nb_main_stem"] <- fixef(mod)[1]+ranef(mod)$Focal[1]+fixef(mod)[7]

## Tiller Nb
mod <- lmer(Tiller_nb ~ Sampling_ID + Block + Treatment + (1+Treatment|Focal), data=trait_monoc)

blup_trait_monoc[which(blup_trait_monoc$Treatment=="C"),"Tiller_nb"] <- fixef(mod)[1]+ranef(mod)$Focal[1]
blup_trait_monoc[which(blup_trait_monoc$Treatment=="S"),"Tiller_nb"] <- fixef(mod)[1]+ranef(mod)$Focal[1]+fixef(mod)[7]+ranef(mod)$Focal[2]

## Leaf N
mod <- lmer(leaf_N ~ Day_leafN_measurement + Block + Treatment + (1+Treatment|Focal), data=trait_monoc)

blup_trait_monoc[which(blup_trait_monoc$Treatment=="C"),"leaf_N"] <- fixef(mod)[1]+ranef(mod)$Focal[1]
blup_trait_monoc[which(blup_trait_monoc$Treatment=="S"),"leaf_N"] <- fixef(mod)[1]+ranef(mod)$Focal[1]+fixef(mod)[7]+ranef(mod)$Focal[2]

## Root prof
mod <- lmer(BE_BoitEng_Hauteur_mm ~ Block + Treatment + (1 |Focal), data=trait_monoc)
summary(mod)

blup_trait_monoc[which(blup_trait_monoc$Treatment=="C"),"BE_BoitEng_Hauteur_mm"] <- fixef(mod)[1]+ranef(mod)$Focal[1]
blup_trait_monoc[which(blup_trait_monoc$Treatment=="S"),"BE_BoitEng_Hauteur_mm"] <- fixef(mod)[1]+ranef(mod)$Focal[1]+fixef(mod)[4]

## Root length
mod <- lmer(SPUR_Squelette_Corrige_mm ~ Block + Treatment + (1+Treatment|Focal), data=trait_monoc)
summary(mod)

blup_trait_monoc[which(blup_trait_monoc$Treatment=="C"),"SPUR_Squelette_Corrige_mm"] <- fixef(mod)[1]+ranef(mod)$Focal[1]
blup_trait_monoc[which(blup_trait_monoc$Treatment=="S"),"SPUR_Squelette_Corrige_mm"] <- fixef(mod)[1]+ranef(mod)$Focal[1]+fixef(mod)[4]+ranef(mod)$Focal[2]

## Root area
mod <- lmer(SURF_Surface_Projetee_mm2 ~ Block + Treatment + (1+Treatment|Focal), data=trait_monoc)
summary(mod)

blup_trait_monoc[which(blup_trait_monoc$Treatment=="C"),"SURF_Surface_Projetee_mm2"] <- fixef(mod)[1]+ranef(mod)$Focal[1]
blup_trait_monoc[which(blup_trait_monoc$Treatment=="S"),"SURF_Surface_Projetee_mm2"] <- fixef(mod)[1]+ranef(mod)$Focal[1]+fixef(mod)[4]+ranef(mod)$Focal[2]

## computing trait averages and trait distances for each mixture
for (i in c(1:nrow(RYT))) {
  geno1 <- strsplit(RYT[i,"Pair_unoriented"], ";")[[1]][1]
  geno2 <- strsplit(RYT[i,"Pair_unoriented"], ";")[[1]][2]
  trt <- RYT[i,"Treatment"]
  
  for (v in c("Leaf_nb_main_stem", "Tiller_nb", "leaf_N","BE_BoitEng_Hauteur_mm", "SPUR_Squelette_Corrige_mm", "SURF_Surface_Projetee_mm2")){
    RYT[i,paste(v, "avg", sep="_")] <- mean(blup_trait_monoc[which(blup_trait_monoc$Treatment==trt & blup_trait_monoc$Genotype%in%c(geno1,geno2)), v])
    RYT[i,paste(v, "dist", sep="_")] <- abs(blup_trait_monoc[which(blup_trait_monoc$Treatment==trt & blup_trait_monoc$Genotype==geno1), v]-blup_trait_monoc[which(blup_trait_monoc$Treatment==trt & blup_trait_monoc$Genotype==geno2), v])
  }
}

##################
#### Running model averaging ####
##################

######  GENERAL PARAMETERS ####
N <- 10L # number of models to be retained among the best models
w <- 0.2 # barplot width space
conversion <- 0.005 # conversion factor to define colors on a rgb scale between [0,1]
cols <- c(rgb(0,150*conversion,50*conversion),rgb(150*conversion,75*conversion,0))  # vector of colors to differenciate aboveground traits (green) and belowground traits (brown)


## Shoot RYT - C Treatment

## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(RYT[which(RYT$Treatment=="C"),c(3,6:11,14:17)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(Shoot_DW_RYT ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05)[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05)),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(row.names(effect_sizes), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]

## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]

## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)
info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

######
### Plots
########

## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(2,2,2,1,2,1,1)]
bgvec <- cols[c(4,2,2,1,4,4,1)]
densvec <- c(20,-9,-9,-9,20,20,-9)
anglevec <- rep(45,length(colvec))

pdf("outputs/plots/Model_selection_Shoot_DW_RYT_C.pdf", height=3, width=5, pointsize = 0.2)

par(mfrow=c(1,2))
par(mar=c(4,5,1,1), oma=c(5,2,1,2), lwd=1.5)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-1,1), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)
abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)

axis(2, las=1, labels=c("Root surf.","Root surf.","Root len.","# leaves","Root len.","# tiller", "Leaf N"), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.5) 
axis(1, at=seq(-1,1,by=0.5), cex.lab=1.3, cex.axis=1.3)

## 2nb plot: variable importance
par(mar=c(4,1,1,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1, density = rev(densvec), angle=rev(anglevec), ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)

## plot legend

xlegend = -0.45
ylegend =-(nrow(effect_sizes)+2)*(1-0.7)

legend(x=xlegend,y=ylegend,legend=c(" "," "), pch=c(19,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA)
legend(xlegend+0.02,ylegend*1.15, legend="Mean", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.4,y=ylegend,legend=c(" "," "), pch=c(21,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA, density=c(0,20))
legend(xlegend+0.02+0.4,ylegend*1.15, legend="Dist", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.7,y=ylegend,legend=c(" "," "), fill=c(cols[1],cols[2]), bty="n", col=c("black","black"), border=c(NA,NA), xjust=0, cex=1.5, xpd=NA)
legend(x=xlegend+0.7,y=ylegend*1.01,legend=c("Aboveground","Belowground"), fill=c(NA,NA), bty="n", col=c(NA,NA), border=c(NA,NA), xjust=0,cex=1.4, xpd=NA)

dev.off()

## Shoot RYT - S Treatment

## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(RYT[which(RYT$Treatment=="S"),c(3,6:11,14:17)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(Shoot_DW_RYT ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05)[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05)),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(row.names(effect_sizes), decreasing = F),]
effect_sizes <- effect_sizes[order(abs(effect_sizes$Estimate), decreasing = T),]

## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]

## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

######
### Plots
########

## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(2,2,1,1,2,2,1,1)]
bgvec <- cols[c(2,2,4,1,4,4,4,1)]
densvec <- c(-9,-9,20,-9,20,20,20,-9)
anglevec <- rep(45,length(colvec))

pdf("outputs/plots/Model_selection_Shoot_DW_RYT_S.pdf", height=3, width=5, pointsize = 0.2)

par(mfrow=c(1,2))
par(mar=c(4,5,1,1), oma=c(5,2,1,2), lwd=1.5)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-2,2), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)
abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)
axis(2, las=1, labels=c("Root surf.","Root len.","# leaves","# leaves","Root surf.","Root len.","# tillers", "Leaf N"), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.5) 
axis(1, at=seq(-2,2,by=0.5), cex.lab=1.3, cex.axis=1.3)

## 2nb plot: variable importance
par(mar=c(4,1,1,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1, density = rev(densvec), angle=rev(anglevec), ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)

## plot legend

xlegend = -0.7
ylegend =-(nrow(effect_sizes)+2)*(1-0.7)

legend(x=xlegend,y=ylegend,legend=c(" "," "), pch=c(19,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA)
legend(xlegend+0.02,ylegend*1.15, legend="Mean", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.4,y=ylegend,legend=c(" "," "), pch=c(21,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA, density=c(0,20))
legend(xlegend+0.02+0.4,ylegend*1.15, legend="Dist", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.7,y=ylegend,legend=c(" "," "), fill=c(cols[1],cols[2]), bty="n", col=c("black","black"), border=c(NA,NA), xjust=0, cex=1.5, xpd=NA)
legend(x=xlegend+0.7,y=ylegend*1.01,legend=c("Aboveground","Belowground"), fill=c(NA,NA), bty="n", col=c(NA,NA), border=c(NA,NA), xjust=0,cex=1.4, xpd=NA)

dev.off()


## Root DW - C Treatment

## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(RYT[which(RYT$Treatment=="C"),c(4,6:11,14:17)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(Root_DW_RYT ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05)[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05)),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(row.names(effect_sizes), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]

## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]

## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

######
### Plots
########

## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(1,2,1,1,2,2,2)]
bgvec <- cols[c(1,2,4,1,4,4,2)]
densvec <- c(-9,-9,20,-9,20,20,-9)
anglevec <- rep(45,length(colvec))


pdf("outputs/plots/Model_selection_Root_DW_RYT_C.pdf", height=3, width=5, pointsize = 0.2)

par(mfrow=c(1,2))
par(mar=c(4,5,1,1), oma=c(5,2,1,2), lwd=1.5)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-1,1), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)
abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)
axis(2, las=1, labels=c("# leaves", "Root surf.", "Leaf N", "Leaf N", "Root surf.","Root len.", "Root len."), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.5) 
axis(1, at=seq(-1,1,by=0.5), cex.lab=1.3, cex.axis=1.3)

## 2nb plot: variable importance
par(mar=c(4,1,1,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1, density = rev(densvec), angle=rev(anglevec), ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)

## plot legend

xlegend = -0.7
ylegend =-(nrow(effect_sizes)+2)*(1-0.7)

legend(x=xlegend,y=ylegend,legend=c(" "," "), pch=c(19,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA)
legend(xlegend+0.02,ylegend*1.15, legend="Mean", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.4,y=ylegend,legend=c(" "," "), pch=c(21,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA, density=c(0,20))
legend(xlegend+0.02+0.4,ylegend*1.15, legend="Dist", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.7,y=ylegend,legend=c(" "," "), fill=c(cols[1],cols[2]), bty="n", col=c("black","black"), border=c(NA,NA), xjust=0, cex=1.5, xpd=NA)
legend(x=xlegend+0.7,y=ylegend*1.01,legend=c("Aboveground","Belowground"), fill=c(NA,NA), bty="n", col=c(NA,NA), border=c(NA,NA), xjust=0,cex=1.4, xpd=NA)

dev.off()


## Root RYT - S Treatment

## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(RYT[which(RYT$Treatment=="S"),c(4,6:11,14:17)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(Root_DW_RYT ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05)[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05)),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(row.names(effect_sizes), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]

## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]

## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

######
### Plots
########

## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(2,1,1,2,1,2,2,1,1)]
bgvec <- cols[c(2,1,1,2,4,4,4,1,1)]
densvec <- c(-9,-9,-9,-9,20,20,20,-9,20)
anglevec <- rep(45,length(colvec))


pdf("outputs/plots/Model_selection_Root_DW_RYT_S.pdf", height=3, width=5, pointsize = 0.2)

par(mfrow=c(1,2))
par(mar=c(4,5,1,1), oma=c(5,2,1,2), lwd=1.5)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-1,1), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)
abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)
axis(2, las=1, labels=c("Root surf.", "# tillers", "Leaf N","Root len.", "# leaves","Root len.","Root surf.","# leaves", "# tillers"), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.5) 
axis(1, at=seq(-1,1,by=0.5), cex.lab=1.3, cex.axis=1.3)

## 2nb plot: variable importance
par(mar=c(4,1,1,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1, density = rev(densvec), angle=rev(anglevec), ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)

## plot legend

xlegend = -0.7
ylegend =-(nrow(effect_sizes)+2)*(1-0.7)

legend(x=xlegend,y=ylegend,legend=c(" "," "), pch=c(19,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA)
legend(xlegend+0.02,ylegend*1.15, legend="Mean", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.4,y=ylegend,legend=c(" "," "), pch=c(21,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA, density=c(0,20))
legend(xlegend+0.02+0.4,ylegend*1.15, legend="Dist", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.7,y=ylegend,legend=c(" "," "), fill=c(cols[1],cols[2]), bty="n", col=c("black","black"), border=c(NA,NA), xjust=0, cex=1.5, xpd=NA)
legend(x=xlegend+0.7,y=ylegend*1.01,legend=c("Aboveground","Belowground"), fill=c(NA,NA), bty="n", col=c(NA,NA), border=c(NA,NA), xjust=0,cex=1.4, xpd=NA)

dev.off()

## Total RYT - C Treatment

## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(RYT[which(RYT$Treatment=="C"),c(5:11,14:17)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(Total_DW_RYT ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05)[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05)),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(row.names(effect_sizes), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]

## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]

## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

######
### Plots
########

## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(2,2,1,2,1,1,2)]
bgvec <- cols[c(4,2,1,2,1,1,4)]
densvec <- c(20,-9,-9,-9,-9,-9,20)
anglevec <- rep(45,length(colvec))

pdf("outputs/plots/Model_selection_Total_DW_RYT_C.pdf", height=3, width=5, pointsize = 0.2)

par(mfrow=c(1,2))
par(mar=c(4,5,1,1), oma=c(5,2,1,2), lwd=1.5)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-2,2), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)
abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)
axis(2, las=1, labels=c("Root surf.","Root surf.","# leaves", "Root len.", "Leaf N", "# tillers", "Root len."), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.5) 
axis(1, at=seq(-2,2,by=0.5), cex.lab=1.3, cex.axis=1.3)

## 2nb plot: variable importance
par(mar=c(4,1,1,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1, density = rev(densvec), angle=rev(anglevec), ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)

## plot legend

xlegend = -0.7
ylegend =-(nrow(effect_sizes)+2)*(1-0.7)

legend(x=xlegend,y=ylegend,legend=c(" "," "), pch=c(19,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA)
legend(xlegend+0.02,ylegend*1.15, legend="Mean", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.4,y=ylegend,legend=c(" "," "), pch=c(21,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA, density=c(0,20))
legend(xlegend+0.02+0.4,ylegend*1.15, legend="Dist", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.7,y=ylegend,legend=c(" "," "), fill=c(cols[1],cols[2]), bty="n", col=c("black","black"), border=c(NA,NA), xjust=0, cex=1.5, xpd=NA)
legend(x=xlegend+0.7,y=ylegend*1.01,legend=c("Aboveground","Belowground"), fill=c(NA,NA), bty="n", col=c(NA,NA), border=c(NA,NA), xjust=0,cex=1.4, xpd=NA)


dev.off()

## Total RYT - S Treatment

## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(RYT[which(RYT$Treatment=="S"),c(5:11,14:17)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(Total_DW_RYT ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05)[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05)),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(row.names(effect_sizes), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]

## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]

## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

######
### Plots
########

## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(2,2,1,1,1,1,2)]
bgvec <- cols[c(2,2,1,1,4,4,4)]
densvec <- c(-9,-9,-9,-9,20,20,20)
anglevec <- rep(45,length(colvec))

pdf("outputs/plots/Model_selection_Total_DW_RYT_S.pdf", height=3, width=5, pointsize = 0.2)

par(mfrow=c(1,2))
par(mar=c(4,5,1,1), oma=c(5,2,1,2), lwd=1.5)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-2,2), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)
abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)
axis(2, las=1, labels=c("Root surf.", "Root len.", "# tillers", "# leaves","# leaves","Leaf N","Root surf."), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.5) 
axis(1, at=seq(-2,2,by=0.5), cex.lab=1.3, cex.axis=1.3)

## 2nb plot: variable importance
par(mar=c(4,1,1,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1, density = rev(densvec), angle=rev(anglevec), ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)

## plot legend

xlegend = -0.7
ylegend =-(nrow(effect_sizes)+2)*(1-0.7)

legend(x=xlegend,y=ylegend,legend=c(" "," "), pch=c(19,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA)
legend(xlegend+0.02,ylegend*1.15, legend="Mean", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.4,y=ylegend,legend=c(" "," "), pch=c(21,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA, density=c(0,20))
legend(xlegend+0.02+0.4,ylegend*1.15, legend="Dist", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.7,y=ylegend,legend=c(" "," "), fill=c(cols[1],cols[2]), bty="n", col=c("black","black"), border=c(NA,NA), xjust=0, cex=1.5, xpd=NA)
legend(x=xlegend+0.7,y=ylegend*1.01,legend=c("Aboveground","Belowground"), fill=c(NA,NA), bty="n", col=c(NA,NA), border=c(NA,NA), xjust=0,cex=1.4, xpd=NA)

dev.off()

########################
### FOLLOW UP ANALYSIS ON RYT-ROOT TRAITS RELATIONSHIPS 
########################

## Checking the relationship btw avg root surface and RYT on Total biomass RYT
ggplot(RYT, aes(x=SURF_Surface_Projetee_mm2_avg, y=Total_DW_RYT, color=Treatment, shape=Treatment))+
  geom_point()+
  scale_color_manual(values=c("blue","red"))+
  scale_shape_manual(values=c(16,17))+
  labs(x = expression("Average root projected area (mm"^2*")"),
       y = "Total biomass RYT")+
  geom_smooth(method = "lm", fill = NA, show.legend = F)+
  stat_cor(method = "pearson", show.legend = F)+
  theme_bw()
ggsave("outputs/plots/RYT_vs_root_surface_measured_in_monocultures.pdf", dpi=300, height=6, width=8)


### Checking the relationship btw avg root surface and RYT on Total biomass RYT, using avg root surface from mixture rhyzotrons
##1. computing avg root surface in mixture rhyzotrons
root_traits_mix <- all_root_trait[which(all_root_trait$stand=="mixed"),]
for (i in c(1:nrow(root_traits_mix))) {
  root_traits_mix[i,"Pair_unoriented"] <- paste(sort(c(root_traits_mix[i,"Focal"],root_traits_mix[i,"Neighbour"])), collapse = ";")
}

root_traits_mix$Pair_unoriented <- as.factor(root_traits_mix$Pair_unoriented)
blup_SURF_mix <- data.frame(Pair_unoriented=levels(root_traits_mix$Pair_unoriented), Treatment=rep(c("C","S"), each=54))

mod <- lmer(SURF_Surface_Projetee_mm2 ~ Block + Treatment +  (1|Pair_unoriented), data=root_traits_mix)

blup_SURF_mix[which(blup_SURF_mix$Treatment=="C"),"SURF_Surface_Projetee_mm2"] <- fixef(mod)[1]+ranef(mod)$Pair_unoriented[1]
blup_SURF_mix[which(blup_SURF_mix$Treatment=="S"),"SURF_Surface_Projetee_mm2"] <- fixef(mod)[1]+ranef(mod)$Pair_unoriented[1]+fixef(mod)[3]

##2. merging root surface in mixture rhyzotrons with the global data set
RYT <- merge(RYT, blup_SURF_mix, by=c("Pair_unoriented","Treatment"))
RYT$SURF_Surface_Projetee_mm2_diff <- (RYT$SURF_Surface_Projetee_mm2-RYT$SURF_Surface_Projetee_mm2_avg)/RYT$SURF_Surface_Projetee_mm2_avg

##3. plotting the relationship between RYT and root surface measured in mixtures
ggplot(RYT, aes(x=SURF_Surface_Projetee_mm2, y=Total_DW_RYT, color=Treatment, shape=Treatment))+
  geom_point()+
  scale_color_manual(values=c("blue","red"))+
  scale_shape_manual(values=c(16,17))+
  labs(x = expression("Average root projected area (mm"^2*")"),
       y = "Total biomass RYT")+
  geom_smooth(method = "lm", fill = NA, show.legend = F)+
  stat_cor(method = "pearson", show.legend = F)+
  theme_bw()
ggsave("outputs/plots/RYT_vs_root_surface_measured_in_mixtures.pdf", dpi=300, height=6, width=8)


#### checking the relationship between RYT and the average productivity of the two genotypes in monoculture
for (i in c(1:nrow(RYT))) {
  geno1 <- strsplit(RYT[i,"Pair_unoriented"], ";")[[1]][1]
  geno2 <- strsplit(RYT[i,"Pair_unoriented"], ";")[[1]][2]
  trt <- RYT[i,"Treatment"]
  
  for (v in c("Shoot_DW", "Root_DW", "Total_DW")){
    RYT[i,paste(v, "avg", sep="_")] <- mean(blup_monoc[which(blup_monoc$Treatment==trt & blup_monoc$Genotype%in%c(geno1,geno2)), v])
  }
}

modWS <- lm(Total_DW_RYT~Total_DW_avg*Treatment, data=RYT)
anova(modWS)
summary(modWS)

ggplot(RYT, aes(x=Total_DW_avg, y=Total_DW_RYT, color=Treatment, shape=Treatment))+
  geom_point()+
  scale_color_manual(values=c("blue","red"))+
  scale_shape_manual(values=c(16,17))+
  labs(x = "Average biomass in monoculture (mg)",
       y = "Total biomass RYT")+
  geom_smooth(method = "lm", fill = NA, show.legend = F)+
  stat_cor(method = "pearson", show.legend = F)+
  theme_bw()
ggsave("outputs/plots/RYT_vs_monoculture_prod.pdf", dpi=300, height=6, width=8)
## --->  High RYT are obtained when mixing genotypes with low productivity in monoculture

#### checking the relationship monoculture biomass production and root surface
monoc_trait_prod <- merge(blup_monoc, blup_trait_monoc, by=c("Genotype","Treatment"))
monoc_trait_prod$Treatment <- as.factor(monoc_trait_prod$Treatment)
monoc_trait_prod$Treatment <- factor(monoc_trait_prod$Treatment, levels=c("C","S"))

ggplot(monoc_trait_prod, aes(x=SURF_Surface_Projetee_mm2, y=Total_DW, color=Treatment, shape=Treatment))+
  geom_point()+
  scale_color_manual(values=c("blue","red"))+
  scale_shape_manual(values=c(16,17))+
  labs(x = expression("Root projected area (mm"^2*")"),
       y = "Total biomass (mg)")+
  geom_smooth(method = "lm", fill = NA, show.legend = F)+
  stat_cor(method = "pearson", show.legend = F)+
  theme_bw()
ggsave("outputs/plots/monoculture_prod_vs_root_surface.pdf", dpi=300, height=6, width=8)
## ---> Higher root surface in monoculture is associated with increased biomas (both in the WW and WS treatment, and both above and belowground)

### checking the relationship between genotypes'RY and (i) biomass productivity in monoculture, (ii) root projected area in monoculture
for (i in c(1:nrow(blup_mixt))) {
  blup_mixt[i,"focal"] <- as.character(strsplit(blup_mixt[i,"Pair"], ";")[[1]][1])
  blup_mixt[i,"neighbour"] <- as.character(strsplit(blup_mixt[i,"Pair"], ";")[[1]][2])
  blup_mixt[i,"monoc_biomass"] <- blup_monoc[which(blup_monoc$Genotype==blup_mixt[i,"focal"] & blup_monoc$Treatment==blup_mixt[i,"Treatment"]),"Total_DW"]
  blup_mixt[i,"monoc_root_surf"] <- blup_trait_monoc[which(blup_trait_monoc$Genotype==blup_mixt[i,"focal"] & blup_trait_monoc$Treatment==blup_mixt[i,"Treatment"]),"SURF_Surface_Projetee_mm2"]
}

blup_mixt$Treatment <- as.factor(blup_mixt$Treatment)
blup_mixt$Treatment <- factor(blup_mixt$Treatment, levels=c("C","S"))
pl1 <- ggplot(blup_mixt, aes(x=monoc_root_surf, y=Total_DW_RY, color=Treatment, shape=Treatment))+
  geom_point()+
  scale_color_manual(values=c("blue","red"))+
  scale_shape_manual(values=c(16,17))+
  labs(x =  expression("monoculture root projected area (mm"^2*")"),
       y = "biomass RY")+
  geom_smooth(method = "lm", fill = NA, show.legend = F)+
  stat_cor(method = "pearson", show.legend = F)+
  theme_bw()

pl2 <- ggplot(blup_mixt, aes(x=monoc_biomass, y=Total_DW_RY, color=Treatment, shape=Treatment))+
  geom_point()+
  scale_color_manual(values=c("blue","red"))+
  scale_shape_manual(values=c(16,17))+
  labs(x = "monoculture biomass (mg)",
       y = "biomass RY")+
  geom_smooth(method = "lm", fill = NA, show.legend = F)+
  stat_cor(method = "pearson", show.legend = F)+
  theme_bw()

ggarrange(plotlist = list(pl1,pl2),
          nrow=1,
          ncol=2,
          common.legend = T)
ggsave("outputs/plots/RY_vs_monoculture_biomass_&_root_projected_area.pdf", dpi=300, height=5, width=9)

