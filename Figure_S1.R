library(ggvenn)
Daily_venn<-c("yftL","sspA","fimH","fimE","pldA","kgtP","gsiA","rph","cytR","mngA","mutL","acrB","speG","cobU","dosP","gatB","gatZ","mutT","xdhA","ydcI","yqiA","ybaL","yfaQ","yegH","gatA","mppA","uspF","glcF","infB","nlpD","rseB","tdcB")
Hund_Venn<-c("fimE","mreB","araJ","pgsA","hns","paaX","cpdA","ydcI","hdfR","proQ","rho","cspC","rpoB","nudF","cdsA","crp","sspA","clpA","sucC","dacA","rpoC","spoT","dadX","sstT","dctA","dsdC","putA","ybaL","speC","argR","infB","putP","pitA","ppx","rplN","gabP","hisP")

df2<-list(Daily_venn,Hund_Venn)

names(df2)<-c(c("1-day" , "100-day"))
venn.intersect<-Reduce(intersect,df2)

venn.plot<-ggvenn(df2,
       columns = c("1-day" ,"100-day"),
       show_percentage=FALSE,
       set_name_size = 5,
       fill_color = c("#EFC000", "#7E409D"))

FigureS1<-venn.plot+
  annotate("segment", x = 0, y = -0.4,xend = 0, yend = -1, color="black")+
  annotate("text", x = 0.25, y = -1.1, label = venn.intersect[1], fontface="italic")+
  annotate("text", x = 0, y = -1.1, label = venn.intersect[2], fontface="italic")+
  annotate("text", x = -0.25, y = -1.1, label = venn.intersect[3], fontface="italic")+
  annotate("text", x = 0.5, y = -1.1, label = venn.intersect[4], fontface="italic")+
  annotate("text", x = -0.5, y = -1.1, label = venn.intersect[5], fontface="italic")+
  scale_y_continuous(limits=c(-1.3,1.3))
