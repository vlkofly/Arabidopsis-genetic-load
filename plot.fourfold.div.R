setwd("/home/vlkofly/Arabidopsis_science/arenosa/scantools_diversity")
load(".RData")
library(tidyverse)
library(ggrepel)
div<-read.table("diversity.tsv",sep="\t",header=T)
summary(div)
div$ploidy<-as.factor(as.character(div$ploidy))



svg("divaa.alldataset_ploidy.svg")
ggplot(data=div[div$species == "AAtoAA" & div$locus == "fourfold",],aes(x=ploidy,y=Pi,fill=ploidy))+
  geom_boxplot()+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,size=2,position=position_jitterdodge(),alpha=0.8)+
  theme(axis.text.x=element_text(size=20,angle = 45),
        axis.title = element_text(size=23),
        axis.text.y=element_text(size=23),plot.title=element_text(size=25),
        legend.text = element_text(size=18),
        strip.text.x = element_text(size = 18))
dev.off()

load("~/Arabidopsis_science/arenosa/indels/analysis_oct2023_allchrom5x/.RData") # I will get the target pops from here.
th<-theme(axis.text.x=element_blank(),
          axis.title = element_text(size=23),
          axis.title.x = element_blank(),
          axis.text.y=element_text(size=23,color="black"),plot.title=element_text(size=25),
          #legend.text = element_text(size=18),
          strip.text.x = element_text(size = 18),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black",linetype="solid")
)

div.plot<-ggplot(data=div[div$species == "AAtoAA" & div$locus == "fourfold" & div$pop %in% pop20,],aes(x=ploidy,y=Pi*100,fill=ploidy))+
  geom_boxplot()+
  ylab("\u03C0 (x 0.01)")+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,size=3,position=position_jitterdodge(),alpha=0.8)+
  geom_text(aes(label=pop),size=4,check_overlap = T,nudge_x = 0.2)+
  th+ theme (legend.position = "none")
  
svg("diversity.ffinalw3.lab.svg",width=3)
plot(div.plot)
dev.off()



#eventually we decided to show plots for all 65 populations
aa.div<-div[div$species == "AAtoAA" & div$locus == "fourfold",]
summary(aa.div)

svg("diversity.65.nolab.svg",width=3)
ggplot(aa.div[aa.div$locus=="fourfold",],aes(x=ploidy,y=Pi*100,fill=ploidy))+
geom_boxplot()+
  ylab("\u03C0 (x 0.01)")+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,size=3,position=position_jitterdodge(),alpha=0.8)+
  #geom_text(aes(label=pop),size=4,check_overlap = T,nudge_x = 0.2)+ #nolab version
  th+ theme (legend.position = "none")
dev.off()



pt.all<-read.delim("popstat.17.05.23.tsv")
pt.all<-pt.all[pt.all$species == "AAtoAA",]
pt.all$ploidy<-as.factor(as.character(pt.all$ploidy))
pt<-pt.all[pt.all$species %in% c("AAtoAA") & pt.all$pop %in% pop20,]
pt$ploidy<-as.factor(as.character(pt$ploidy))

summary(pt)
# just check if the pi plot is the same, yes, so now plot TajD
tajd.plot<-ggplot(data=pt,aes(x=ploidy,y=Taj.D.4d,fill=ploidy))+
  geom_boxplot()+
  ylab("Tajima's D")+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_text_repel(aes(label=pop),size=3)+
  geom_point(pch=21,size=3,position=position_jitterdodge(),alpha=0.8)+
  th + theme (legend.position = "none")


svg("TajD.nolab.w3.svg",width=3)
print(tajd.plot)
dev.off()

summary(pt.all)

svg("TajD.65.svg",width=3)
ggplot(data=pt.all,aes(x=ploidy,y=Taj.D.4d,fill=ploidy))+
  geom_boxplot()+
  ylab("Tajima's D")+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  #geom_text_repel(aes(label=pop),size=3)+
  geom_point(pch=21,size=3,position=position_jitterdodge(),alpha=0.8)+
  th + theme (legend.position = "none")
dev.off()

# Is it statistically significant?
wilcox.test(pt[pt$ploidy==2,"Taj.D.4d"],pt[pt$ploidy==4,"Taj.D.4d"])
wilcox.test(pt[pt$ploidy==2,'Pi.4d'],pt[pt$ploidy==4,"Pi.4d"])
t.test(pt[pt$ploidy==2,'Pi.4d'],pt[pt$ploidy==4,"Pi.4d"])

wilcox.test(pt.all[pt.all$ploidy==2,"Taj.D.4d"],pt.all[pt.all$ploidy==4,"Taj.D.4d"])
wilcox.test(pt.all[pt.all$ploidy==2,'Pi.4d'],pt.all[pt.all$ploidy==4,"Pi.4d"])
pt.all%>%
group_by(ploidy)%>% summarise(mean(Pi.4d))

0.025/0.021
dp.glm<-glm(Pi.4d~ploidy+X4d.depth,data=pt)
summary(dp.glm)
anova(dp.glm,test="Chisq")

# Is it significant also across the whole dataset?
dp.all.glm<-glm(Pi.4d~ploidy+X4d.depth,data=pt.all)
summary(dp.all.glm)
summary(pt.all)


# Is it significant also with the pops used in Monnahan?
pop300<-c("BEL","BGS","BIH","BRD","CHO","DRA","FOJ","GOR","GUL","HAR","HNE","HNI","HOC","KAS","KOS","KOW","KZL","LAC","MIE","PRE","RFT","RZA","SCH","SNO","SPI","STE","SWA","SZI","TBG","TKO","TRD","TRE","TRT","TZI","VEL","VID","WEK","ZAP")
length(pop300)
pt.300<- pt.all %>%
  filter(pop %in% pop300) %>%
  droplevels()
summary(pt.300)
dim(pt.300)
dp.300.glm<-glm(Pi.4d~ploidy+X4d.depth,data=pt.300)
summary(dp.300.glm)
# Is it significant without interploidy admixed populations removed from focal dataset?
pt$pop
popadmix<-c("ZIT","PHT","TIS","BAL","BUT","CAR","HRN","INE","PAT") #HRA,KAM? in the landscape they are not but in my table they are yesWC
pt.noadmix<-pt %>%
  filter(! pop %in% popadmix) %>%
  droplevels()
summary(pt.noadmix) # only 6 tetraploids
dim(pt.noadmix)
pt.noadmix$pop
dp.noadmix.glm<-glm(Pi.4d~ploidy+X4d.depth,data=pt.noadmix)
summary(dp.noadmix.glm)

# nonadmixed tetraploids

wilcox.test(pt.all[pt.all$pop %in% popadmix,"Pi.4d"], pt.noadmix[pt.noadmix$ploidy == 4,"Pi.4d"]) # admixed vs. non-admixed tetraploids

#' Plot fourfold and zerofold ratio 

summary(pt)
pt$piNS<-pt$Pi.0d/pt$Pi.4d

piNS.plot<-ggplot(data=pt,aes(x=ploidy,y=piNS,fill=ploidy))+
  geom_boxplot()+
  ylab("0-fold \u03C0 / 4-fold \u03C0 ")+
  xlab("Ploidy")+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  #geom_text_repel(aes(label=pop),size=3)+
  geom_point(pch=21,position=position_jitterdodge(),alpha=0.5)+
  th + theme(axis.text.x = element_text(size=23,color="black"),axis.title.x = element_text(size=23))
#+ theme (legend.position = "none")

svg("piNS.w3.svg",width=3)
print(piNS.plot)
dev.off()

svg("labels.svg")
ggplot(data=pt,aes(x=ploidy,y=piNS,fill=ploidy))+
geom_text(y=0.2,label="A B C D E F DFE")+
  th
dev.off()

wilcox.test(pt[pt$ploidy==2,'piNS'],pt[pt$ploidy==4,"piNS"])
