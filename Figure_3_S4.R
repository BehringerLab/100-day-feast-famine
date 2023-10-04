library(stringr)
library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(grid)
library(ggpp)

#Read in massive Biofilm Data File
FF_100_Biofilm<-read.table("~/GitHub/datasets/BiofilmPopulations_100FF.txt",sep='\t',header=FALSE)
FF_100_Biofilm

#Break Biofilm Data down by batch
##Test Data
#Plate1<- FF_100_Biofilm[which(FF_100_Biofilm$V1=='Plate1'),] 
##Real Data
Plate2<- FF_100_Biofilm[which(FF_100_Biofilm$V1=='Plate2'),] 
Plate3<- FF_100_Biofilm[which(FF_100_Biofilm$V1=='Plate3'),] 
Plate4<- FF_100_Biofilm[which(FF_100_Biofilm$V1=='Plate4'),] 
Plate5<- FF_100_Biofilm[which(FF_100_Biofilm$V1=='Plate5'),] 
Plate6<- FF_100_Biofilm[which(FF_100_Biofilm$V1=='Plate6'),] 
Plate7<- FF_100_Biofilm[which(FF_100_Biofilm$V1=='Plate7'),] 
Plate8<- FF_100_Biofilm[which(FF_100_Biofilm$V1=='Plate8'),] 

##Fix headers on individual plate data

names(Plate2)<- as.matrix(Plate2[1,])
Plate2<- Plate2[-1,]
Plate2[]<-lapply(Plate2,function(x) type.convert(as.character(x)))
Plate2

names(Plate3)<- as.matrix(Plate3[1,])
Plate3<- Plate3[-1,]
Plate3[]<-lapply(Plate3,function(x) type.convert(as.character(x)))
Plate3

names(Plate4)<- as.matrix(Plate4[1,])
Plate4<- Plate4[-1,]
Plate4[]<-lapply(Plate4,function(x) type.convert(as.character(x)))
Plate4

names(Plate5)<- as.matrix(Plate5[1,])
Plate5<- Plate5[-1,]
Plate5[]<-lapply(Plate5,function(x) type.convert(as.character(x)))
Plate5

names(Plate6)<- as.matrix(Plate6[1,])
Plate6<- Plate6[-1,]
Plate6[]<-lapply(Plate6,function(x) type.convert(as.character(x)))
Plate6

names(Plate7)<- as.matrix(Plate7[1,])
Plate7<- Plate7[-1,]
Plate7[]<-lapply(Plate7,function(x) type.convert(as.character(x)))
Plate7

names(Plate8)<- as.matrix(Plate8[1,])
Plate8<- Plate8[-1,]
Plate8[]<-lapply(Plate8,function(x) type.convert(as.character(x)))
Plate8

##Melt individual plate data
Plate2Melt<-reshape2::melt(Plate2,id=c("Plate2", "Replicate"))
names(Plate2Melt)<-c("Plate","Replicate","Sample","OD550")
Plate2Melt

Plate3Melt<-reshape2::melt(Plate3,id=c("Plate3", "Replicate"))
names(Plate3Melt)<-c("Plate","Replicate","Sample","OD550")
Plate3Melt

Plate4Melt<-reshape2::melt(Plate4,id=c("Plate4", "Replicate"))
names(Plate4Melt)<-c("Plate","Replicate","Sample","OD550")
Plate4Melt

Plate5Melt<-reshape2::melt(Plate5,id=c("Plate5", "Replicate"))
names(Plate5Melt)<-c("Plate","Replicate","Sample","OD550")
Plate5Melt

Plate6Melt<-reshape2::melt(Plate6,id=c("Plate6", "Replicate"))
names(Plate6Melt)<-c("Plate","Replicate","Sample","OD550")
Plate6Melt

Plate7Melt<-reshape2::melt(Plate7,id=c("Plate7", "Replicate"))
names(Plate7Melt)<-c("Plate","Replicate","Sample","OD550")
Plate7Melt

Plate8Melt<-reshape2::melt(Plate8,id=c("Plate8", "Replicate"))
names(Plate8Melt)<-c("Plate","Replicate","Sample","OD550")
Plate8Melt

###Pull Blanks take average and subtract from wells to remove background by batch, fix zeros for statistics

Plate2BlankOnly<- Plate2Melt[which(Plate2Melt$Sample %like% "Blank",),] 
Plate3BlankOnly<- Plate3Melt[which(Plate3Melt$Sample %like% "Blank",),] 
Plate4BlankOnly<- Plate4Melt[which(Plate4Melt$Sample %like% "Blank",),] 
Plate5BlankOnly<- Plate5Melt[which(Plate5Melt$Sample %like% "Blank",),] 
Plate6BlankOnly<- Plate6Melt[which(Plate6Melt$Sample %like% "Blank",),]
Plate7BlankOnly<- Plate7Melt[which(Plate7Melt$Sample %like% "Blank",),]
Plate8BlankOnly<- Plate8Melt[which(Plate8Melt$Sample %like% "Blank",),]

Plate2Melt$BlankAvg<-median(Plate2BlankOnly$OD550)
Plate3Melt$BlankAvg<-median(Plate3BlankOnly$OD550)
Plate4Melt$BlankAvg<-median(Plate4BlankOnly$OD550)
Plate5Melt$BlankAvg<-median(Plate5BlankOnly$OD550)
Plate6Melt$BlankAvg<-median(Plate6BlankOnly$OD550)
Plate7Melt$BlankAvg<-median(Plate7BlankOnly$OD550)
Plate8Melt$BlankAvg<-median(Plate8BlankOnly$OD550)

Plate2Melt$OD550_Corrected<-Plate2Melt$OD550-Plate2Melt$BlankAvg
Plate3Melt$OD550_Corrected<-Plate3Melt$OD550-Plate3Melt$BlankAvg
Plate4Melt$OD550_Corrected<-Plate4Melt$OD550-Plate4Melt$BlankAvg
Plate5Melt$OD550_Corrected<-Plate5Melt$OD550-Plate5Melt$BlankAvg
Plate6Melt$OD550_Corrected<-Plate6Melt$OD550-Plate6Melt$BlankAvg
Plate7Melt$OD550_Corrected<-Plate7Melt$OD550-Plate7Melt$BlankAvg
Plate8Melt$OD550_Corrected<-Plate8Melt$OD550-Plate8Melt$BlankAvg

#Plate1Melt$OD550_Corrected[Plate1Melt$OD550_Corrected<=0] <- 0.00001
Plate2Melt$OD550_Corrected[Plate2Melt$OD550_Corrected<=0] <- 0.00001
Plate3Melt$OD550_Corrected[Plate3Melt$OD550_Corrected<=0] <- 0.00001
Plate4Melt$OD550_Corrected[Plate4Melt$OD550_Corrected<=0] <- 0.00001
Plate5Melt$OD550_Corrected[Plate5Melt$OD550_Corrected<=0] <- 0.00001
Plate6Melt$OD550_Corrected[Plate6Melt$OD550_Corrected<=0] <- 0.00001
Plate7Melt$OD550_Corrected[Plate7Melt$OD550_Corrected<=0] <- 0.00001
Plate8Melt$OD550_Corrected[Plate8Melt$OD550_Corrected<=0] <- 0.00001

#Grab WT controls from each batch for normalization

Plate2WTOnly<- Plate2Melt[which(Plate2Melt$Sample %like% "WT",),] 
Plate3WTOnly<- Plate3Melt[which(Plate3Melt$Sample %like% "WT",),] 
Plate4WTOnly<- Plate4Melt[which(Plate4Melt$Sample %like% "WT",),] 
Plate5WTOnly<- Plate5Melt[which(Plate5Melt$Sample %like% "WT",),] 
Plate6WTOnly<- Plate6Melt[which(Plate6Melt$Sample %like% "WT",),]
Plate7WTOnly<- Plate7Melt[which(Plate7Melt$Sample %like% "WT",),]
Plate8WTOnly<- Plate8Melt[which(Plate8Melt$Sample %like% "WT",),]


Plate2Melt$WTAvg<-median(Plate2WTOnly$OD550_Corrected)
Plate3Melt$WTAvg<-median(Plate3WTOnly$OD550_Corrected)
Plate4Melt$WTAvg<-median(Plate4WTOnly$OD550_Corrected)
Plate5Melt$WTAvg<-median(Plate5WTOnly$OD550_Corrected)
Plate6Melt$WTAvg<-median(Plate6WTOnly$OD550_Corrected)
Plate7Melt$WTAvg<-median(Plate7WTOnly$OD550_Corrected)
Plate8Melt$WTAvg<-median(Plate8WTOnly$OD550_Corrected)


Plate2Melt$OD550_Normal<-Plate2Melt$OD550_Corrected/Plate2Melt$WTAvg
Plate3Melt$OD550_Normal<-Plate3Melt$OD550_Corrected/Plate3Melt$WTAvg
Plate4Melt$OD550_Normal<-Plate4Melt$OD550_Corrected/Plate4Melt$WTAvg
Plate5Melt$OD550_Normal<-Plate5Melt$OD550_Corrected/Plate5Melt$WTAvg
Plate6Melt$OD550_Normal<-Plate6Melt$OD550_Corrected/Plate6Melt$WTAvg
Plate7Melt$OD550_Normal<-Plate7Melt$OD550_Corrected/Plate7Melt$WTAvg
Plate8Melt$OD550_Normal<-Plate8Melt$OD550_Corrected/Plate8Melt$WTAvg

##Combine plates
BiofilmPlates_withBlanks<-rbind(Plate2Melt,Plate3Melt,Plate4Melt,Plate5Melt,Plate6Melt,Plate7Melt,Plate8Melt)
BiofilmPlates_withBlanks


#Remove Blanks for plotting
BiofilmPlates<- BiofilmPlates_withBlanks[which(! BiofilmPlates_withBlanks$Sample %like% "Blank",),] 
BiofilmPlates

#Add and annotate Treatment value for Ancestor, Engineered Deletions, and Timepoint
BiofilmPlates$Treatment<-"EngClone"

D300<-BiofilmPlates$Sample %like% "300"
BiofilmPlates$Treatment[D300] <- "D300"
D900<-BiofilmPlates$Sample %like% "900"
BiofilmPlates$Treatment[D900] <- "D900"
WT<-BiofilmPlates$Sample %like% "WT"
BiofilmPlates$Treatment[WT] <- "WT"


#Create Data Table with just measurements of evolved populations (Time=900 days), Engineered Deletions, and controls
BiofilmPlatesD900<- BiofilmPlates[which(BiofilmPlates$Treatment %like% "EngClone"|BiofilmPlates$Treatment %like% "WT"|BiofilmPlates$Treatment %like% "D900"),] 

#Organize Levels for plotting
BiofilmPlatesD900$Sample<-factor(BiofilmPlatesD900$Sample, levels=c("WT","fimE","paaX","ydcI","hdfR","lrhA","hdfR/fimE","lrhA/fimE","P501_900","P503_900","P504_900","P506_900","P508_900","P509_900","P511_900","P512_900","P513_900","P515_900","P516_900","P518_900","P520_900","P521_900","P523_900","P524_900"))

##Split Engineered Deletions and Evolved Populations
BiofilmPlatesD900noEng<-BiofilmPlatesD900[which(BiofilmPlatesD900$Treatment=="WT"|BiofilmPlatesD900$Treatment=="D900"),]
BiofilmPlatesD900noEvo<-BiofilmPlatesD900[which(BiofilmPlatesD900$Sample=="WT"|BiofilmPlatesD900$Sample=="fimE"|BiofilmPlatesD900$Sample=="paaX"|BiofilmPlatesD900$Sample=="ydcI"|BiofilmPlatesD900$Sample=="hdfR"|BiofilmPlatesD900$Sample=="hdfR/fimE"),]

#Label genetic background for Evolved Populations and Control
Ans<-BiofilmPlatesD900noEng$Sample %in% c("WT")
WT<-BiofilmPlatesD900noEng$Sample %in% c("P501_900","P503_900","P508_900","P512_900","P513_900","P515_900","P520_900","P524_900")
MMR <-BiofilmPlatesD900noEng$Sample %in% c("P504_900","P506_900","P509_900","P511_900","P516_900","P518_900","P521_900","P523_900")
BiofilmPlatesD900noEng$Background[Ans] <- "Ancestor"
BiofilmPlatesD900noEng$Background[WT] <- "WT"
BiofilmPlatesD900noEng$Background[MMR] <- "MMR"


###Perform Wilcox Test
compare_means(data=BiofilmPlatesD900noEng, OD550_Normal~Background,ref.group="Ancestor", method="wilcox.test")

##Annotate Significance from Wilcox Test
BF_WT_Control<-grobTree(text_grob("***",x=0.49,y=0.9,hjust=0),gp=gpar(col="black",fontsize=10))
BF_WT_MMR<-grobTree(text_grob("***",x=0.805,y=0.9,hjust=0),gp=gpar(col="black",fontsize=10))


#Plot Figure 3A
BF_Plot_noEng<-ggplot(BiofilmPlatesD900noEng, aes(x=Background,y=OD550_Normal, col=Background,fill=Background))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.2)+
  stat_summary(fun.y="median", geom="point")+
  annotation_custom(BF_WT_Control)+
  annotation_custom(BF_WT_MMR)+
  ylab("Relative Biofilm")+
  xlab("Sample")+
  scale_y_continuous(limits=c(0,9),breaks=c(0,2,4,6,8))+
  geom_hline(yintercept=1,linetype="dashed",color="#AAAAAA")+
  scale_color_manual(values=c("#AAAAAA","#8a3da1","#bd08a0"))+
  scale_fill_manual(values=c("#AAAAAA22","#8a3da122","#bd08a022"))+
  theme_bw()+
  theme(legend.position="none", axis.title.y = element_text(face="bold",size=10),axis.title.x=element_blank(),legend.background = element_rect(color="black"), legend.title = element_text(face="bold",size=8), legend.title.align=0.5, axis.text.x = element_blank())
BF_Plot_noEng


#Perform Wilcox test for Engineered Deletions
compare_means(data=BiofilmPlatesD900noEvo, OD550_Normal~Sample,ref.group="WT", method="wilcox.test")

##Annotate Significance for Engineered Deletions
BF_WT_fimE<-grobTree(text_grob("***",x=0.25,y=0.9,hjust=0),gp=gpar(col="black",fontsize=10))
#BF_WT_paaX<-grobTree(text_grob("*",x=0.58,y=0.9,hjust=0),gp=gpar(col="black",fontsize=10))
BF_hdfRfimE<-grobTree(text_grob("***",x=0.89,y=0.9,hjust=0),gp=gpar(col="black",fontsize=10))


#Plot Figure 3C
BF_Plot_noEvo<-ggplot(BiofilmPlatesD900noEvo, aes(x=Sample,y=OD550_Normal, col=Sample, fill=Sample))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.2)+
  stat_summary(fun.y="median", geom="point")+
  annotation_custom(BF_WT_fimE)+
  annotation_custom(BF_hdfRfimE)+
  ylab("Relative Biofilm")+
  xlab("Sample")+
  scale_y_continuous(limits=c(0,9),breaks=c(0,2,4,6,8))+
  geom_hline(yintercept=1,linetype="dashed",color="#AAAAAA")+
  scale_color_manual(values=c("#AAAAAA","#00AFBB", "#E7B800", "#FC4E07","#4E84C4","#008080"),labels=c("Ancestor","fimE","paaX","ydcI","hdfR","hdfR/fimE"))+
  scale_fill_manual(values=c("#AAAAAA22","#00AFBB22", "#E7B80022", "#FC4E0722","#4E84C422","#00808022"),labels=c("Ancestor","fimE","paaX","ydcI","hdfR","hdfR/fimE"))+
  theme_bw()+
  theme(legend.position="none", axis.title.y = element_text(face="bold",size=10),axis.title.x=element_blank(),legend.background = element_rect(color="black"), legend.title = element_text(face="bold",size=8), legend.title.align=0.5, axis.text.x = element_blank())
BF_Plot_noEvo

####Pull in Motility Data
Motility<-read.table("~/GitHub/datasets/MotilityFinal.txt",sep="\t",header=TRUE)

#Organize Motility Data for Plotting
Motility$Population<-factor(Motility$Population, levels=c("WT","fimE","paaX","ydcI","hdfR","lrhA","hdfR/fimE","501","503","504","506","508","509","511","512","513","515","516","518","520","521","523","524"))

#Calculate Motility Statistics for Ancestor Control and normalize all samples by WT ancestor
Motility_WT<-Motility[which(Motility$Population=="WT"),]
Motility_WT
as.data.frame(Motility_WT %>% group_by(Population) %>% summarise(swim.avg=mean(Diff_Radius),swim.med=median(Diff_Radius), swim.nor=median(Mot_Normal)))
Motility$Mean_WT_AllTimepoints<-3.75394
Motility$NewMot_Normal<-Motility$Diff_Radius/Motility$Mean_WT_AllTimepoints

###Split data into Evolved populations and Engineered deletions
Motility_noEng<-Motility[which(Motility$Timepoint=="WT"|Motility$Timepoint=="Day_900"),]
Motility_noEvo<-Motility[which(Motility$Population=="WT"|Motility$Population=="fimE"|Motility$Population=="paaX"|Motility$Population=="ydcI"|Motility$Population=="hdfR"|Motility$Population=="hdfR/fimE"),]

#Calculate and annotate significance for Evolved Populations
compare_means(NewMot_Normal~Background,Motility_noEng, ref.group="Control", method="t.test")

Mot_noEngStat<-grobTree(textGrob("*",x=0.495,y=0.85,hjust=0))

#Plot Evolved Populations
Mot_Plot_noEng<-ggplot(data=Motility_noEng, aes(x=Background,y=NewMot_Normal,col=Background,fill=Background))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.2)+
  stat_summary(fun.y=median, geom="point", pch=18, size=4) +
  annotation_custom(Mot_noEngStat)+
  ylab("Relative Motility")+
  scale_y_continuous(limits=c(0,3.5))+
  geom_hline(yintercept=1,linetype="dashed",color="#AAAAAA")+
  scale_color_manual(values=c("#AAAAAA","#8a3da1","#bd08a0"))+
  scale_fill_manual(values=c("#AAAAAA22","#8a3da122","#bd08a022"))+
  theme_bw()+
  theme(legend.position="none", axis.title = element_text(face="bold",size=10),legend.background = element_rect(color="black"), legend.title = element_text(face="bold",size=8), legend.title.align=0.5, axis.text.x = element_text(size=10))
Mot_Plot_noEng

###Calculate and annotate significance for Eng. Deletions
compare_means(Mot_Normal~Treatment,Motility_noEvo, ref.group="WT",method="wilcox.test")

Mot_paaX<-grobTree(textGrob("**",x=0.41,y=0.85,hjust=0))
Mot_ydcI<-grobTree(textGrob("*",x=0.575,y=0.85,hjust=0))
Mot_hdfR<-grobTree(textGrob("*",x=0.735,y=0.85,hjust=0))

##Plot Engineered Deletions
Mot_Plot_noEvo<-ggplot(data=Motility_noEvo, aes(x=Population,y=Mot_Normal,col=Population,fill=Population))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.2)+
  stat_summary(fun.y=median, geom="point") +
  annotation_custom(Mot_paaX)+
  annotation_custom(Mot_ydcI)+
  annotation_custom(Mot_hdfR)+
  ylab("Relative Motility")+
  xlab("Clone")+
  scale_y_continuous(limits=c(0,3.5))+
  scale_x_discrete(labels=c("Ancestor","fimE","paaX","ydcI","hdfR","hdfR/fimE"))+
  geom_hline(yintercept=1,linetype="dashed",color="#AAAAAA")+
  scale_color_manual(values=c("#AAAAAA","#00AFBB", "#E7B800", "#FC4E07","#4E84C4","#008080"),labels=c("Ancestor","fimE","paaX","ydcI","hdfR","hdfR/fimE"))+
  scale_fill_manual(values=c("#AAAAAA22","#00AFBB22", "#E7B80022", "#FC4E0722","#4E84C422","#00808022"),labels=c("Ancestor","fimE","paaX","ydcI","hdfR","hdfR/fimE"))+
  theme_bw()+
  theme(legend.position="none", axis.title = element_text(face="bold",size=10),legend.background = element_rect(color="black"), legend.title = element_text(face="bold",size=8), legend.title.align=0.5, axis.text.x = element_text(size=10))
Mot_Plot_noEvo


#Combine plots and make Final Figure 3
Left_plot_Top<-plot_grid(BF_Plot_noEng,Mot_Plot_noEng, nrow=2, labels=c("A","B"),rel_heights = c(0.85,1))
Left_plot_Bottom<-plot_grid(BF_Plot_noEvo,Mot_Plot_noEvo, nrow=2, labels=c("C","D"),rel_heights = c(0.85,1))

plot_grid(Left_plot_Top,Left_plot_Bottom, rel_widths = c(0.8, 1))


###Expand Plots for Supplemental Figure S4
Biofilm_AllPopulations<-ggplot(BiofilmPlatesD900noEng, aes(x=Sample,y=OD550_Normal, col=Background,fill=Background))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.2)+
  stat_summary(fun.y="median", geom="point")+
  #annotation_custom(BF_WT_Control)+
  #annotation_custom(BF_WT_MMR)+
  ylab("Relative Biofilm")+
  xlab("Sample")+
  scale_y_continuous(limits=c(0,9),breaks=c(0,2,4,6,8))+
  geom_hline(yintercept=1,linetype="dashed",color="#AAAAAA")+
  scale_color_manual(values=c("#AAAAAA","#8a3da1","#bd08a0"))+
  scale_fill_manual(values=c("#AAAAAA22","#8a3da122","#bd08a022"))+
  theme_bw()+
  theme(legend.position="none", axis.title.y = element_text(face="bold",size=10),axis.title.x=element_blank(),legend.background = element_rect(color="black"), legend.title = element_text(face="bold",size=8), legend.title.align=0.5, axis.text.x = element_blank())


Motility_AllPopulations<-ggplot(data=Motility_noEng, aes(x=Population,y=NewMot_Normal,col=Background,fill=Background))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.2)+
  stat_summary(fun.y=median, geom="point", pch=18, size=4) +
  annotation_custom(Mot_noEngStat)+
  ylab("Relative Motility")+
  scale_y_continuous(limits=c(0,3.5))+
  geom_hline(yintercept=1,linetype="dashed",color="#AAAAAA")+
  scale_color_manual(values=c("#AAAAAA","#8a3da1","#bd08a0"))+
  scale_fill_manual(values=c("#AAAAAA22","#8a3da122","#bd08a022"))+
  theme_bw()+
  theme(legend.position="none", axis.title = element_text(face="bold",size=10),legend.background = element_rect(color="black"), legend.title = element_text(face="bold",size=8), legend.title.align=0.5, axis.text.x = element_text(size=10))

FigureS4<-plot_grid(Biofilm_AllPopulations,Motility_AllPopulations,nrow=2)
