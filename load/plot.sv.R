setwd("/home/vlkofly/Arabidopsis_science/arenosa/svs")
library(ggrepel)
library(tidyverse)
library(scales)

gen.indpop20<-read.table("arenosa_sv_counts_v2.txt",header=T,sep="\t")
summary(gen.indpop20)
gen.indpop20$ploidy <- as.factor(as.character(gen.indpop20$ploidy))
gen.indpop20$id
gen.indpop20$pop<-gsub("([[:alpha:]]*)_.*",replacement = "\\1",gen.indpop20$id)
gen.indpop20$pop<-gsub("AA126-BD-E","DRE",gen.indpop20$pop)
gen.indpop20$pop<-gsub("AA268","MON",gen.indpop20$pop)
gen.indpop20$pop<-gsub("AA358","BOR",gen.indpop20$pop)
gen.indpop20$pop<-gsub("AA397.*","OPP",gen.indpop20$pop)
gen.indpop20$pop<-gsub("AA452.*","KRY",gen.indpop20$pop)
gen.indpop20$pop<-gsub("AA506.*","MIL",gen.indpop20$pop)
gen.indpop20$pop -> sv.pops


saveRDS(sv.pops,"sv.pops.rds")

# change the names of populations
th<-theme(axis.title = element_text(size=23),
          axis.text.y=element_text(size=23,color="black"),
          axis.text.x=element_text(size=23,color="black"),
          plot.title=element_text(size=25),
          legend.text = element_text(size=18),
          strip.text.x = element_text(size = 18),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black",linetype="solid")
)




gen.indpop20$scaled_interg_n <-scale(gen.indpop20$interg_n)
gen.indpop20$scaled_exon_n <-scale(gen.indpop20$exon_n)

lm.sv.2<-lm(scaled_exon_n~scaled_interg_n,data=gen.indpop20[gen.indpop20$ploidy == 2,])
lm.sv.4<-lm(scaled_exon_n~scaled_interg_n,data=gen.indpop20[gen.indpop20$ploidy == 4,])
plot(scaled_interg_n~scaled_exon_n,data=gen.indpop20[gen.indpop20$ploidy == 4,])
summary(lm(scaled_interg_n~scaled_exon_n*ploidy,data=gen.indpop20))

plot(lm.sv.2)
summary(lm.sv.2)
summary(lm.sv.4)

slodip<-coef(lm.sv.2)[2]
slotet<-coef(lm.sv.4)[2]

slodip<-round(slodip,digits=2)
slotet<-round(slotet,digits=2) # why is the slope smaller in tetraploids? wtf
slodip
slotet
dim(gen.indpop20)

fig2a<-ggplot(gen.indpop20,aes(x=interg_n,y=exon_n,fill=ploidy,color=ploidy))+
  ylab("Exonic SVs")+
  xlab("Intergenic SVs")+ # call it synonymous sites and for additive load I will call it synonymous allelels
  geom_smooth(method="glm",linewidth=1,method.args = list(family = "gaussian"))+
  geom_point(pch=21,size=4,alpha=0.7,color="black")+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  scale_x_continuous(breaks = breaks_pretty(n=4),label = label_number(scale = 1e-3,suffix = "K"))+
  scale_y_continuous(breaks = breaks_pretty(n=4),label = label_number(scale = 1e-3,suffix = "K"))+
  geom_text_repel(aes(label=pop),min.segment.length = 0, size=3,point.padding = 0.01,box.padding = 0.6,max.overlaps = 15,direction = "both",color="black")+
  #geom_text(aes(2e3,0.5e3,label=paste("\u03B2","=",as.character(slodip),sep="")),color="#1e90ff",size=7)+
  #geom_text(aes(2.3e3,0.5e3,label=paste("/",as.character(slotet),sep="")),color="#ffa500",size=7)+
  th

svg("sv2.load.svg")
print(fig2a)
dev.off()

saveRDS(fig2a,"sv.load.fig.rds")
#' Test the significance of interaction

sv0.glm<-glm(exon_n~interg_n+ploidy+degen.indpop20h,family="poisson",data=gen.indpop20)
sv1.glm<-glm(exon_n~interg_n*ploidy+degen.indpop20h,family="poisson",data=gen.indpop20)
anova(sv0.glm,sv1.glm,test="LRT")
summary(sv1.glm)


#' get the total number of svs
gen.indpop20$total_n<-gen.indpop20$interg_n+gen.indpop20$exon_n + gen.indpop20$intron_n
ggplot(gen.indpop20,aes(x=total_n,y=degen.indpop20h,colour=ploidy))+
  geom_point(size=4,pch=16)+
  geom_text(aes(label=pop),size=2.5,check_overlap = F,nudge_x = +6)+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  #geom_hline(yintercegen.indpop20=20,colour="red")+
  #geom_hline(yintercegen.indpop20=55,colour="red")+
  #geom_hline(yintercegen.indpop20=22,colour="blue")+
  xlab("Number of SVs")+
  ylab("Degen.indpop20h of coverage")+
  th
 

wilcox.test(gen.indpop20[gen.indpop20$ploidy==2,"total_n"],gen.indpop20[gen.indpop20$ploidy==4,"total_n"])
gen.indpop20 %>%
  group_by(ploidy)%>% summarise(mean(total_n))
