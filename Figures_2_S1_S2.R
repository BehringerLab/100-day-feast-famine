###Script to Analyze and Produce Figure 2 and Supplementary Figures 1 and 2
library(ggplot2)
library(ggrepel)
library(cowplot)
library(grid)
#library(ggpubr)
library(scatterpie)
library(tidyr)
library(dplyr)
library(ggpp)
detach(package:ggpubr,unload=TRUE)

##Call in timeline of genes overrepresented for mutations
timelineAll<-read.table("~/GitHub/datasets/TimelineAll_A.txt",sep="\t",header=TRUE,check.names = FALSE)

#Summarize Mutation counts
timelineAll %>% group_by(GeneType) %>% count(MutGene)
aggregate(data = timelineAll,MutGene ~ GeneType,function(MutGene) length(unique(MutGene)))


## Calculate ANOVA to test for mutational order
fit0 = lm(Time~MutGene, timelineAll)
TimelineAll_AOV<-anova(fit0)
TimelineAll_AOV

###Plot Mutation Timeline
TimelinePlot_All<-ggplot(timelineAll,aes(x=Time/100,y=MutGene,group=MutGene,color=GeneType))+
  stat_summary(geom="point",fun="mean",size=3,shape=20)+
  stat_summary(geom="errorbar",fun.data="mean_cl_boot")+
  theme_classic()+
  scale_x_continuous(limits=c(0,9),breaks=c(1,2,3,4,5,6,7,8,9))+
  xlab("Time (d x 100)")+
  ylab("")+
  scale_color_manual(values=c("#00798c","#d1495b","#edae49","#66a182","#2e4057"))+
  scale_y_discrete(limits=c("fimE","mreB","araJ","pgsA","hns","paaX","cpdA","ydcI","hdfR","proQ","rho","cspC","rpoB","nudF","cdsA","crp","sspA","clpA","sucC","dacA","rpoC","spoT","dadX","sstT","dctA","dsdC","putA","ybaL","speC","argR","infB","putP","pitA","ppx","rplN","gabP","hisP"))+
  theme(axis.text.y = element_text(face="italic"), axis.title.x=element_text(size=10,face="bold"),legend.position = c(0.8,0.22),legend.title = element_blank())
TimelinePlot_All


##Format dataset to plot pie charts
OverRepMut.d<-timelineAll %>% group_by(GeneType,MutClass) %>% summarise(cnt = n()) %>% mutate(freq = round(cnt / sum(cnt), 3),pos=cumsum(freq)-0.5*freq)
OverRepMut.d<-as.data.frame(OverRepMut.d)
OverRepMut.d$pos2<- 1-OverRepMut.d$pos


OverRepPieChart<-ggplot(OverRepMut.d, aes(x="", y=freq, fill=MutClass)) +
  geom_bar(stat="identity", width=1, color="black") +
  geom_text(aes(x = "", y = pos2,label=cnt),size=2,position = position_nudge_center(x=0.7))+
  coord_polar(theta='y')+
  scale_fill_manual(values=c("pink","yellow","lightgreen","lightblue","#666666","lightgrey","orange","purple"),name="Mutation:",labels=c("indel","IS-element","large deletion","nonsense","PV","SNP"))+
  theme_void()+theme(strip.text.x = element_text(face="bold"),legend.position ="bottom", legend.background = element_rect(color="black"), legend.title = element_blank(), legend.title.align=0.5, legend.text=element_text(size=8),legend.key.size = unit(0.2, 'cm'))+guides(fill=guide_legend(nrow=3,byrow=TRUE))+facet_wrap(~GeneType, nrow=5)

OverRepPieChart

Figure2<-plot_grid(TimelinePlot_All,OverRepPieChart)
Figure2


##Supplemental Figure 1
#Split mutation list by WT and MMR- Background
TL_WT<-timelineAll[which(timelineAll$Background=="WT"),]
TL_MMR<-timelineAll[which(timelineAll$Background=="MMR-"),]

#Calculate ANOVA to test for mutational order
fit0WT = lm(Time~MutGene, TL_WT)
fit0MMR = lm(Time~MutGene, TL_MMR)

anova(fit0WT)
anova(fit0MMR)

###Use means to make discrete y axis
TL_WT %>% group_by(MutGene) %>% summarise(mean.WT=mean(Time)) %>% arrange(mean.WT) %>% print(n=40)

TimelinePlot_WT<-ggplot(TL_WT,aes(x=Time,y=MutGene,group=MutGene,color=GeneType,shape=GeneType,fill=GeneType))+
  stat_summary(geom="errorbar",fun.data="mean_cl_boot")+
  stat_summary(geom="point",fun="mean",size=3, color="black")+
  theme_classic()+
  scale_x_continuous(limits=c(0,900),breaks=c(90,200,300,400,500,600,700,800,900))+
  xlab("Time (d)")+
  ylab("")+
  scale_color_manual(values=c("#00798c","#d1495b","#edae49","#66a182","#2e4057"), name="Gene class")+
  scale_fill_manual(values=c("#00798c","#d1495b","#edae49","#66a182","#2e4057"), name="Gene class")+
  scale_shape_manual(values=c(21,22,23,24,25), name="Gene class")+
  scale_y_discrete(limits=c("proQ","fimE","paaX","hdfR","mreB","hns","ydcI","cdsA","crp","cpdA","sspA","rho","rpoC","sucC","dacA","clpA","putA","spoT","sstT","ybaL","infB","argR","cspC","dsdC","gabP","ppx","rplN","rpoB","putP","pitA","dctA","dadX"))+
  theme(axis.text.y = element_text(face="italic"),axis.text.x = element_text(angle=45,hjust=1), axis.title.x=element_text(size=10,face="bold"),legend.position = c(0.8,0.15),legend.title = element_blank())

TimelinePlot_WT

###Use means to make discrete y axis
TL_MMR %>% group_by(MutGene) %>% summarise(mean.WT=mean(Time)) %>% arrange(mean.WT) %>% print(n=40)

TimelinePlot_MMR<-ggplot(TL_MMR,aes(x=Time,y=MutGene,group=MutGene,color=GeneType,shape=GeneType,fill=GeneType))+
  stat_summary(geom="errorbar",fun.data="mean_cl_boot")+
  stat_summary(geom="point",fun="mean",size=3, color="black")+
  theme_classic()+
  scale_x_continuous(limits=c(0,900),breaks=c(90,200,300,400,500,600,700,800,900))+
  xlab("Time (d)")+
  ylab("")+
  scale_color_manual(values=c("#00798c","#d1495b","#edae49","#66a182","#2e4057"), name="Gene class")+
  scale_fill_manual(values=c("#00798c","#d1495b","#edae49","#66a182","#2e4057"), name="Gene class")+
  scale_shape_manual(values=c(21,22,23,24,25), name="Gene class")+
  scale_y_discrete(limits=c("dsdC","fimE","mreB","cpdA","araJ","pgsA","ydcI","hns","cspC","spoT","rho","gltA","paaX","hdfR","sucC","rpoB","clpA","dacA","argR","nudF","dadX","proQ","crp","rpoC","dctA","sspA","sstT","cdsA","pitA","putA","ybaL","speC","infB","putP","ppx","rplN","gabP","hisP"))+
  theme(axis.text.y = element_text(face="italic"),axis.text.x = element_text(angle=45,hjust=1), axis.title.x=element_text(size=10,face="bold"),legend.position = c(0.8,0.15),legend.title = element_blank())


TimelinePlot_MMR

FigureS1<-plot_grid(TimelinePlot_WT,TimelinePlot_MMR, labels=c("A","B"))
FigureS1

###Supplemental Figure 2
####Subset Regulators
timelineAll.Reg<-timelineAll[which(timelineAll$GeneType=="Regulator"),]
timelineAll.Reg %>% group_by(Mut_Type) %>% summarise(mut.mean=mean(Time))
timelineAll.Reg[which(timelineAll.Reg$Mut_Type=="NS"),]
#Regulator ANOVA
fit_reg = lm(Time~MutClass, timelineAll.Reg)
anova(fit_reg)


TimelinePlot_Regulators<-ggplot(timelineAll.Reg,aes(x=Time,y=MutGene,group=MutGene))+
  stat_summary(geom="point",fun="mean",size=3,shape=20,color="#66a182")+
  stat_summary(geom="errorbar",fun.data="mean_cl_boot",color="#66a182")+
  theme_classic()+
  scale_x_continuous(breaks=c(90,200,300,400,500,600,700,800,900))+
  xlab("Time (d)")+
  ylab("")+
  scale_y_discrete(limits=c("fimE","hns","paaX","ydcI","hdfR","rho","cspC","rpoB","crp","sspA","rpoC","spoT","dsdC","putA","argR", "infB","rplN"))+
  theme(axis.text.y = element_text(face="italic"), axis.title.x=element_text(size=10,face="bold"))

TimelinePlot_Regulators

TimelinePlot_Regulator_Mut<-ggplot(timelineAll.Reg,aes(x=Time,y=MutClass,group=MutClass))+
  stat_summary(geom="point",fun="mean",size=3,shape=20,color="#66a182")+
  stat_summary(geom="errorbar",fun.data="mean_cl_boot",color="#66a182")+
  theme_classic()+
  scale_x_continuous(limits=c(0,900),breaks=c(90,200,300,400,500,600,700,800,900))+
  xlab("Time (d)")+
  ylab("")+
  scale_y_discrete(limits=c("IS","NS","indel","snp"),labels=c("IS-element","Nonsense","Indel","Nonsynonymous/PV"))+
  theme(axis.text.y = element_text(face="italic"), axis.title.x=element_text(size=10,face="bold"))

TimelinePlot_Regulator_Mut

####Subset Regulators
timelineAll.Trans<-timelineAll[which(timelineAll$GeneType=="Transporter"),]

fit_trans = lm(Time~MutGene, timelineAll.Trans)

anova(fit_trans)

TimelinePlot_Transporters<-ggplot(timelineAll.Trans,aes(x=Time,y=MutGene,group=MutGene))+
  stat_summary(geom="point",fun="mean",size=3,shape=20,color="#2e4057")+
  stat_summary(geom="errorbar",fun.data="mean_cl_boot",color="#2e4057")+
  theme_classic()+
  scale_x_continuous(limits=c(0,900),breaks=c(90,200,300,400,500,600,700,800,900))+
  xlab("Time (d)")+
  ylab("")+
  scale_y_discrete(limits=c("araJ","sstT","dctA","ybaL","putP","pitA","gabP","hisP"))+
  theme(axis.text.y = element_text(face="italic"), axis.title.x=element_text(size=10,face="bold"))

TimelinePlot_Transporters

TimelinePlot_Transporters_Mut<-ggplot(timelineAll.Trans,aes(x=Time,y=Mut_Type,group=Mut_Type))+
  stat_summary(geom="point",fun="mean",size=3,shape=20,color="#2e4057")+
  stat_summary(geom="errorbar",fun.data="mean_cl_boot",color="#2e4057")+
  theme_classic()+
  scale_x_continuous(limits=c(0,900),breaks=c(90,200,300,400,500,600,700,800,900))+
  xlab("Time (d)")+
  ylab("")+
  scale_y_discrete(limits=c("indel","NS","snp"), labels=c("Indel","Nonsense","Nonsynonymous/PV"))+
  theme(axis.text.y = element_text(face="italic"), axis.title.x=element_text(size=10,face="bold"))

TimelinePlot_Transporters_Mut

##Envelope
timelineAll.Env<-timelineAll[which(timelineAll$GeneType=="Membrane"),]
timelineAll.Env
fit_Envelope= lm(Time~MutGene, timelineAll.Env)

anova(fit_Envelope)

TimelinePlot_Envelope<-ggplot(timelineAll.Env,aes(x=Time,y=MutGene,group=MutGene))+
  stat_summary(geom="point",fun="mean",size=3,shape=20,color="#d1495b")+
  stat_summary(geom="errorbar",fun.data="mean_cl_boot",color="#d1495b")+
  theme_classic()+
  scale_x_continuous(limits=c(0,900),breaks=c(90,200,300,400,500,600,700,800,900))+
  xlab("Time (d)")+
  ylab("")+
  scale_y_discrete(limits=c("mreB","pgsA","cdsA","dacA"))+
  theme(axis.text.y = element_text(face="italic"), axis.title.x=element_text(size=10,face="bold"))

TimelinePlot_Envelope

TimelinePlot_Envelope_Mut<-ggplot(timelineAll.Env,aes(x=Time,y=MutClass,group=MutClass))+
  stat_summary(geom="point",fun="mean",size=3,shape=20,color="#d1495b")+
  stat_summary(geom="errorbar",fun.data="mean_cl_boot",color="#d1495b")+
  theme_classic()+
  scale_x_continuous(limits=c(0,900),breaks=c(90,200,300,400,500,600,700,800,900))+
  xlab("Time (d)")+
  ylab("")+
  scale_y_discrete(limits=c("indel","snp","NS","IS","Large_Deletion"),labels=c("Indel","Nonsynonymous/PV","Nonsense","IS-element","Large Deletion"))+
  theme(axis.text.y = element_text(face="italic"), axis.title.x=element_text(size=10,face="bold"))

FigureS2<-plot_grid(TimelinePlot_Regulators,TimelinePlot_Regulator_Mut,TimelinePlot_Transporters,TimelinePlot_Transporters_Mut,TimelinePlot_Envelope,TimelinePlot_Envelope_Mut,nrow=3)
