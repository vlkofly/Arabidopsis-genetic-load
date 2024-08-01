# Analyse polyDFE results
getwd()
setwd("/home/vlkofly/Arabidopsis_science/arenosa/polyDFE_arenosa/downsample8")
source("/home/vlkofly/Programs/polyDFE/postprocessing.R")

library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggrepel)
library(scales)
files<-list.files(pattern="\\.all.models.txt$")
filesA<-list.files(pattern="\\A.all.models.txt$")
filesC<-list.files(pattern="\\C.all.models.txt$")
length(filesA)
length(filesC)
filesAC<-c(filesA,filesC) # work only with A and C models

pops.f<-unique(gsub("_.*","",files)) # upgrade
pops.f<-as.character(pops$V1)
pops<-pops[-7] # HRA removed
#trial variable
p="Alpine4X"

#' These commands are to get familiar with the polyDFE data structure and functions
parseOutput(filesC[2])
A.params<-data.frame(matrix(ncol=3,nrow=0))
for (f in filesA) {
babA<-parseOutput(f)
p<-sapply(babA, function(e) e$values)[[3]][c(6:8)] # now I need average from all populations
A.params<-rbind(A.params,p)
}

names(A.params)<-c("S_bar","b","S_max")
head(A.params)
summary(A.params)

# get the shape parameter b from C+r-eps
C.params<-data.frame(matrix(ncol=3,nrow=0))
for (f in filesC) {
  babC<-parseOutput(f)
  p<-sapply(babC, function(e) e$values)[[3]][c(6:8)] # now I need average from all populations, the 3 model is C+r-eps
  C.params<-rbind(C.params,p)
}
names(C.params)<-c("S_bar","b","S_max")
head(C.params)
summary(C.params)
svg("C+r-eps.shape.param.svg")
ggplot(C.params,aes(x=b))+
  ggtitle("Histogram of shape parameter of DFE C+r-eps")+
  geom_histogram(fill="white",colour=1)+
 theme(axis.text.x=element_text(size=15),
            axis.title = element_text(size=23),
            axis.text.y=element_text(size=15),plot.title=element_text(size=15),
            legend.text = element_text(size=18))
dev.off()



#analyse DFE across populations modified for focal dataset
#' vlkofly@vlkofly-ThinkPad-X13-Gen-1:~/Arabidopsis_science/arenosa/polyDFE_arenosa/downsample8$ grep -f ../../snp_load/snp.norm/pop20.txt ploidy.pop.tsv

#' I will get 21 populations. Should I use only these? The missing populations are:
#'   DRE2, FEL4, MIL2, OLO2, RIB2, SEN4

#' select threshold and select 27 pops to have the same number as in sfs

#' So based on SNP depth I will include 6 following populations:
pop20<-read.table("../../snp_load/snp.norm/pop20.txt") # change this to my focal populations
extra.pops<-c("VEL","TRE","GUL","KAS","ING","LAC","SUB") # KRM has strange DFE, so there will be less diploids
remove.pops<-c("DRE","FEL","MIL","OLO","RIB","SEN","ZID") # ZID is missing in the output
pop20<-pop20[!pop20$V1 %in% remove.pops,]
pops<-c(pop20,extra.pops) # OK 27 populations here, but there is X diploids and X tetraploids
p="BAB"
#' But there are several population output files missing
#' I will rather select 8 diploids and 8 tetraploids
#' 
all.models<-data.frame(matrix(ncol=15,nrow=0))

###' The main function starts here ###
#' I change the script on 7.2.2024 to focus only on a subset of populations that were used to estimate genetic load
#' I need to do likelihood ratio test for each population to get the best model for each if it is possible
#' Ideally it will be the same 
#' Then I need to do bootstrapping of DFE
#' Let's take only A and C models and do this compareModels(estBAB.C,estBAB.A,nested=T)$LRT

popdfe<-function(p){
  pf<-grep(p,filesAC,value = T)
  # sometimes there are no three models
  m<-c(parseOutput(pf[1]),parseOutput(pf[2])) # parse AC models
  print(p)
  #print(m)
  grads<-sapply(m,function(e) e$criteria) # get the gradient, it should be far from 1
  names(grads)<-sapply(m,getModelName)
  model<-sapply(m,getModelName)
  #alpha<-unlist(sapply(m,function(e) e$alpha))
  S_bar<-sapply(m,function(e) e$values)[[7]][6] #  scale parameter of C +r -eps
  b<-sapply(m,function(e) e$values)[[7]][7]  # shape parameter
  #S_max<-sapply(m,function(e) e$values)[[3]][8]  # works only for A
  cm<-as.data.frame(compareModels(m)$AIC) # get log likelihoods and AIC
  # what if I want to compare only A+r-eps and C+r-eps
  # compareModels(m[7],m[8],nested = T)$LRT
  cm<-cbind(model,cm,grads,S_bar,b,"pop" = rep(p,8))
  #cnC<-colnames(getDiscretizedDFE(m[[8]])) # get both discrete DFEs
  #cnA<-colnames(getDiscretizedDFE(m[[4]])) # 
  dfe<-as.data.frame(t(sapply(m,getDiscretizedDFE)))
  colnames(dfe)<-cn
  cm<-cbind(cm,dfe)
  all.models<<-rbind(all.models,cm)
  #cm[order(- cm$`log lk`),]
  best<-cm[which.max(cm$`log lk`),] # get the best model according likelihood
  # get index of the model so later you can analyze the DFE
  best.i<- which(rownames(cm) == (rownames(best))) # index of the best model
  dfe<-(getDiscretizedDFE(m[[best.i]]))
  
  # for the purpos of first look lets print the discrete dfe, alpha and model
  n<- paste(p,rownames(best),"grad:", round(best$grads,4),"logL:", round(best$`log lk`,1))
  barplot(dfe,main = n)
  #ggplot(dfe,aes(rownames(dfe)))+
   # geom_bar(stat="identity")
}


popdfe("BAB")
pdf("focal.pops/discreteDFE.pdf",onefile = T)
lapply(pops,popdfe)

dev.off()

###gather all models
lapply(pops,popdfe)

#and analyse it
all.models$model<-as.factor(rownames(all.models))
all.models$model<-as.factor(gsub("eps[0-9]*",replacement = "eps",all.models$model))
summary(all.models)

pop.anal<-unique(all.models$pop)
intersect(pop.anal,pops)
setdiff(pops,pop.anal) # there are 8 populations missing
?setdiff

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
          axis.line = element_line(color="black",linetype="solid"))

svg("likelihoods.all.models.svg")
ggplot(all.models,aes(x=model,y=`log lk`))+
  geom_boxplot()+
  th
dev.off()

all.r.models<-all.models[(grep("\\+ r",all.models$model)),]
svg("focal.pops/likelihoods.r.models.svg")
ggplot(all.r.models,aes(x=model,y=`log lk`))+
  geom_boxplot()+
  geom_point()+
  th
dev.off()
svg("gradient.r.zoom.models.svg")
ggplot(all.r.models[all.r.models$grads,],aes(x=model,y=grads))+
  geom_boxplot()+
  th
dev.off()

summary(all.r.models$grads)



# proportion of beneficial
summary(all.r.models[all.r.models$model == 'A + r - eps',])

# difference between diploids and tetraploids
#ploidy.pop<-read.table("../first.look/ploidy.pop.tsv",header=T,sep="\t")
ploidy.pop<-read.table("ploidy.pop.tsv",header=T,sep="\t")

summary(ploidy.pop)
all.r.models<-(merge(all.r.models,ploidy.pop))

summary(all.r.models)
all.r.models$ploidy<-as.factor(as.character(all.r.models$ploidy))
svg("focal.pops/likelihoods.ploidy.r.models.svg")
ggplot(all.r.models,aes(x=model,y=`log lk`,fill=ploidy))+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,size=1,position=position_jitterdodge(),alpha=0.5)+
  geom_boxplot()+
  th
dev.off()
# plot one boxplot with ploidy point colored
svg("focal.pops/likelihoods.ploidy.r.models.labels.svg")
ggplot(all.r.models,aes(x=model,y=`log lk`,fill=ploidy))+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_boxplot()+
  geom_point(pch=21,size=1,position=position_jitterdodge(),alpha=0.5)+
  geom_text_repel(aes(label=pop),size=3)+
  th+
  theme(axis.text.x=element_text(size=15))
dev.off()


# outlier diploid population
all.r.models[all.r.models$`log lk` < -300 & all.r.models$ploidy == "2",]

svg("focal.pops/gradient.ploidy.r.models.svg")
ggplot(all.r.models,aes(x=model,y=(grads),fill=ploidy))+
  geom_boxplot()+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  th+
  theme(axis.text.x=element_text(size=15))
dev.off()
all.r.models[all.r.models$grads > 0.5,]

#gradients per model
tapply(all.r.models$grads,all.r.models$model,summary)








### discrete DFE 
# first I will need to melt the data set to get the discrete DFE into one column 
bestdfe<-all.r.models[all.r.models$model == "C + r - eps",]
bestdfe<-all.r.models[all.r.models$model == "A + r - eps",]
summary(bestdfe)
meltdfe<-melt(bestdfe[,c(9:16)],id.vars = "ploidy")
summary(meltdfe)

th3<-theme(axis.text.x=element_text(size=13),
          axis.title = element_text(size=23),
          axis.text.y=element_text(size=20),plot.title=element_text(size=25),
          legend.text = element_text(size=18),
          )


svg("focal.pops/DFE.C.svg") #C vary based on best model
ggplot(meltdfe, aes(x=variable,y=value,colour=ploidy))+ # 
  geom_boxplot()+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(position=position_jitterdodge(),alpha=0.5)+
  xlab("Ne*s")+
  ylab("proportion")+
  th+
  theme(axis.text.x=element_text(size=15))
  
dev.off()

svg("focal.pops/DFE.C.del.var.svg") #This will be the figure that goes to the panel
ggplot(meltdfe[!meltdfe$variable %in% c("(0, 1)","(1, 10)","10 <"),], aes(x=variable,y=value,fill=ploidy))+ # 
  geom_boxplot()+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,size=1,position=position_jitterdodge(),alpha=0.5)+
  xlab("Ne*s")+
  ylab("Proportion")+
  th+
  theme(axis.text.x=element_text(size=20))

dev.off()

svg("focal.pops/DFE.C.ben.lab.svg") #This will be the figure that goes to the panel
ggplot(meltdfe[meltdfe$variable %in% c("(0, 1)","(1, 10)","10 <"),], aes(x=variable,y=value,fill=ploidy))+
  geom_boxplot()+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,size=1,position=position_jitterdodge(),alpha=0.5)+
  
  xlab("Ne*s")+
  
  ylab("Proportion")+
  th+
  theme(axis.text.x=element_text(size=20))

dev.off()

svg("focal.pops/DFE.del.del.lab.svg") 
ggplot(bestdfe, aes(x=ploidy,y=`< -100`,fill=ploidy))+ # 
  geom_boxplot()+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,position=position_jitterdodge(),alpha=0.5)+
  xlab("Ploidy")+
  ylab("Proportion of highly deleterious")+
  geom_text_repel(aes(label=pop),size=3)+
  th+
  theme(axis.text.x=element_text(size=15))

dev.off()

svg("focal.pops/DFE.del.100-10.lab.svg") 
ggplot(bestdfe, aes(x=ploidy,y=`(-100, -10)`,fill=ploidy))+ # 
  geom_boxplot()+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,position=position_jitterdodge(),alpha=0.5)+
  xlab("Ploidy")+
  ylab("Proportion of (-100, -10) deleterious")+
  geom_text_repel(aes(label=pop),size=3)+
  th+
  theme(axis.text.x=element_text(size=15))

dev.off()

svg("focal.pops/DFE.del.-1-0.lab.svg") 
ggplot(bestdfe, aes(x=ploidy,y=`(-1, 0)`,fill=ploidy))+ # 
  geom_boxplot()+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,position=position_jitterdodge(),alpha=0.5)+
  xlab("Ploidy")+
  ylab("Proportion of mild (-1, 0) deleterious")+
  geom_text_repel(aes(label=pop),size=3)+
  th+
  theme(axis.text.x=element_text(size=15))

dev.off()

svg("focal.pops/DFE.ben01.lab.svg") 
ggplot(bestdfe, aes(x=ploidy,y=`(1, 10)`,fill=ploidy))+ # 
  geom_boxplot()+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,position=position_jitterdodge(),alpha=0.5)+
  xlab("Ploidy")+
  ylab("Proportion of mildly benficial")+
  geom_text_repel(aes(label=pop),size=3)+
  th+
  theme(axis.text.x=element_text(size=15))

dev.off()

# is there any difference between individual categories between diploids and tetraploids?
wilcox.test(meltdfe[meltdfe$variable == "< -100" & meltdfe$ploidy == 2,"value"],meltdfe[meltdfe$variable == "< -100" & meltdfe$ploidy == 4,"value"])
### W = 131, p-value = 0.008893 yes there is difference tetraploids have less highely deleterious 
wilcox.test(meltdfe[meltdfe$variable == "(-100, -10)" & meltdfe$ploidy == 2,"value"],meltdfe[meltdfe$variable == "(-100, -10)" & meltdfe$ploidy == 4,"value"])
### W = 44, p-value = 0.0595 no there is no difference tetraploids have less highely deleterious 

wilcox.test(meltdfe[meltdfe$variable == "(-10, -1)" & meltdfe$ploidy == 2,"value"],meltdfe[meltdfe$variable == "(-10, -1)" & meltdfe$ploidy == 4,"value"])
### W = 34, p-value = 0.01461 yes there is difference tetraploids have less highely deleterious 
# the difference in fitness effect is non-significant, while in the highly positive mutation is significant

svg("DFE.C.beneficial.points.svg")
ggplot(meltdfe[meltdfe$variable %in% c("(0, 1)","(1, 10)","10 <"),], aes(x=variable,y=value,colour=ploidy))+ # 
  geom_boxplot()+
  geom_point(pch=2,size=0.5,position=position_jitterdodge(),alpha=0.5)+
  #geom_text(aes(label=pop),size=2.5,check_overlap = F,nudge_y = -2)+ there is no info about the population
  xlab("Ne*s")+
  ylab("proportion")+
  th3
dev.off()


all.r.models[all.r.models$model == "A + r - eps",c("pop","alpha")]


# check shape parameter

str(bestdfe)

svg("focal.pops/shape.param.ploidy.w3.svg",width=3) #
ggplot(bestdfe,aes(y=b,x=ploidy,fill=ploidy))+
  ylab("Shape parameter (b)")+
  geom_boxplot()+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  #scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,position=position_jitterdodge(),alpha=0.5)+
  #geom_text_repel(aes(label=pop),size=3)+
  th
dev.off()

# is it statistically significant?
wilcox.test(bestdfe[bestdfe$ploidy == 2,"b"],bestdfe[bestdfe$ploidy == 4,"b"])
#W = 59, p-value = 0.2746

# mean selection coefficients
svg("focal.pops/S.ploidy.w3.svg",width=3)
ggplot(bestdfe,aes(y=S_bar,x=ploidy,fill=ploidy))+
  ylab("Selection coefficient (S)")+
  geom_boxplot()+
  #geom_text(aes(label=pop),size=2.5,check_overlap = F)+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_y_continuous(breaks = breaks_pretty(n=4),label = label_number(scale = 1e-3,suffix = "K"))+
  geom_point(pch=21,position=position_jitterdodge(),alpha=0.5)+
  th
dev.off()

# is it statistically significant?
wilcox.test(bestdfe[bestdfe$ploidy == 2,"S_bar"],bestdfe[bestdfe$ploidy == 4,"S_bar"])
#W = 45, p-value = 0.0672



#bootstrap the data
unique (all.r.models$pop)


##### solution of the error 
getModelName<-function(estimates)
{
  estimates$values<-unlist(estimates$values) # added this into the function
  modelName = estimates$model
  if (estimates$model == "A")
  {
    if (!"S_max" %in% estimates$estimated && estimates$values["S_max"] == 0) # the problem is in the accession of named vector
    {
      modelName = paste(modelName, "del")
    }
  }
  if (estimates$model == "B" || estimates$model == "C")
  {
    if (!"p_b" %in% estimates$estimated && estimates$values["p_b"] == 0)
    {
      modelName = paste(modelName, "del")
    }
  }
  pos = grep("r ", names(estimates$values))
  if (!"r" %in% estimates$estimated 
      && isTRUE(all.equal(estimates$values[pos], 
                          rep(1, length(pos)), 
                          check.attributes = FALSE)))
  {
    modelName = paste(modelName, "- r")
  } else
  {
    modelName = paste(modelName, "+ r")
  }
  if (!"eps_an" %in% estimates$estimated && estimates$values["eps_an"] == 0)
  {
    modelName = paste(modelName, "- eps")
  } else
  {
    modelName = paste(modelName, "+ eps")
  }
  
  return(modelName)
}
