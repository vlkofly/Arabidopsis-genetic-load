##'  Analysis of sfs based on fourfold and zerofold sites 31.1.2024 
#'  
#'
#' The sfs generated from subset vcf where tetraploids have only four individuals and diploids 8 individuals
#' /storage/brno12-cerit/home/filip_kolar/brno3/ERC/arenosa_mapped_to_arenosa/filtering_results/annotate.snps/sep.pop.annotated/gen.counts.norm/subset
#' used script get.sfs.sh
#' parallel "cat {}*zerofold.sfs > {}.zerofold.sfs" :::: pops.txt this generates genome wide sfs from scaffolds
#' parallel "cat {}*fourfold.sfs > {}.fourfold.sfs" :::: pops.txt
#' 
#' scp filip_kolar@storage-brno12-cerit.metacentrum.cz:brno3/ERC/arenosa_mapped_to_arenosa/filtering_results/annotate.snps/sep.pop.annotated/gen.counts.norm/subset/sfs_zerofold/pop/* .

#' I will read the data but I will merge the scaffolds so I have one file per population and then I can sum 
#' each column to get the final sfs.
#' This is how the sfs files look like
#' 379465,10678,5236,3400,2422,1980,1658,1427,1262,1074,1029,1002,1091
#' 

#' 
#' ##' I need to do two things within this script:
#' b) fold the sfs and plot it per ploidy
# 
#' if I wanted to do sfs with the same length for diploids and tetraploids I would need to subset zerofold and fourfold sites from these subset vcf files:
#' /storage/brno12-cerit/home/filip_kolar/brno3/ERC/arenosa_mapped_to_arenosa/filtering_results/annotate.snps/sep.pop.annotated/gen.counts.norm/subset
#'

setwd("/home/vlkofly/Arabidopsis_science/arenosa/tuomas_sfs/sfs16chroms/")
#Functions for folding the SFS and estimating Tajima's D
fold_sfs <- function(sfs){
  out <- NA
  j <- 1
  for(i in 1:(length(sfs-1))){
    if(i <= length(sfs)/2){
      if(i == length(sfs)-i) out[j] <- sfs[i]
      else out[j] <- sfs[i] + sfs[length(sfs)-i]
      j <- j + 1
    }
  }
  out
}

fold_sfs(sfsHIGH[-1])


est_var <- function(n, S){
  a1 <- sum(1/1:(n-1))
  a2 <- sum(1/(1:(n-1))^2)
  b1 <- (n+1)/(3.0*(n-1))
  b2 <- 2*((n*n)+n+3.0)/(9*n*(n-1))
  c1 <- b1-1/a1
  c2 <- b2-(n+2)/(a1*n)+(a2/(a1*a1))
  e1 <- c1/a1
  e2 <- c2/((a1*a1)+a2)
  sqrt(e1*S+e2*S*(S-1))
}

est_tajD <- function(sfs){
  n <- length(sfs)-1
  S <- sum(sfs[-c(1,length(sfs))])
  tW <- S/sum(1/1:(n-1))
  w <- seq(0,n)/n
  tP <- n/(n-1)*2*sum(sfs*w*(1-w))
  (tP-tW)/est_var(n,S)
}
sfs<-sfsHIGH # you should include invariants here
est_pi <- function(sfs){
  n <- length(sfs)-1
  w <- seq(0,n)/n
  tP <- n/(n-1)*2*sum(sfs*w*(1-w))
  tP/sum(sfs)
}

est_wat <- function(sfs){
  n <- length(sfs)-1
  S <- sum(sfs[-c(1,length(sfs))])
  tW <- S/sum(1/1:(n-1))
  tW/sum(sfs)
}


pop20<-as.character(gen.indpop20$pop)


pop="BAB"

fold_sfs(sfsHIGH)
library(tidyverse)
plotsfs<- function (pop) {
  # HIGH effect sites sfs = zerofold
  sfsHIGH<- read.table(paste("zerofold/",pop,".zerofold.sfs",sep=""),sep=",")
  sfsHIGH<- apply(sfsHIGH,2,sum) # sum the frequencies across genome
  sfsHIGH<-sfsHIGH[-1] # remove invariants
  sfsHIGH.fold<-fold_sfs(sfsHIGH) # fold
  sfsHIGH.fold.prop<-sfsHIGH.fold/sum(sfsHIGH.fold)
  sfsHIGH<-cbind ("cat" = rep("zerofold",length(sfsHIGH.fold)),"tot"=sfsHIGH.fold,"prop"=sfsHIGH.fold.prop,"freq"=seq(1,length(sfsHIGH.fold)))
  # WEAK effect sites sfs = fourfold
  sfsWEAK<- read.table(paste("fourfold/",pop,".fourfold.sfs",sep=""),sep=",")
  sfsWEAK<- apply(sfsWEAK,2,sum) # sum the frequencies across genome
  #invWEAK<-sfsWEAK[1]
  sfsWEAK<-sfsWEAK[-1] # remove invariants
  sfsWEAK.fold<-fold_sfs(sfsWEAK) # fold
  sfsWEAK.fold.prop<-sfsWEAK.fold/sum(sfsWEAK.fold)
  sfsWEAK<-cbind ("cat" = rep("fourfold",length(sfsWEAK.fold)),"tot"=sfsWEAK.fold,"prop"=sfsWEAK.fold.prop,"freq"=seq(1,length(sfsWEAK.fold)))

  
  sfs.eff<-as.data.frame(rbind(sfsHIGH,sfsWEAK))
  sfs.eff$freq<-as.factor(as.integer(sfs.eff$freq))
  sfs.eff$prop<-as.numeric(sfs.eff$prop)
  sfs.eff$cat<-as.factor(sfs.eff$cat)
  sfs.eff$pop<-rep(pop,dim(sfs.eff)[1])
  th<-theme(axis.text.x=element_text(size=15),
            axis.title = element_text(size=23),
            axis.text.y=element_text(size=23),plot.title=element_text(size=25),
            legend.text = element_text(size=18))
  ggplot(sfs.eff,aes(x=freq,y=prop,fill=cat))+
    geom_bar(position = "dodge",stat="identity")+
    ggtitle(pop)+
    th
  return(sfs.eff) # switch this on if you want to generate the diploid tetraploid sfs
  
}
#plotsfs("BAB")
pop<-read.table("pops.txt")
pop<-pop$V1
# plot folded sfs of all populations
pdf("sfs.folded.pdf",onefile = T)
lapply(pop,plotsfs)
dev.off()





#'**folded sfs from the pop20**
#'
pop20<-(read.table("pop20.txt",header = F))[,1]
pdf("sfs.folded.pop20.pdf",onefile = T)
lapply(pop20,plotsfs)
dev.off()


# Now I need to generate average tetraploid sfs and average diploid sfs
# to do that I need matrix of all tet and dip sfs switch the function so it returns the sfs.eff table
plotsfs("BAB")
load("../snp_load/snp.sep/.RData") # load snp load data to get ploidy
pop20.dip<-as.character(gen.indpop20[gen.indpop20$ploidy == 2,"pop"])
pop20.tet<-as.character(gen.indpop20[gen.indpop20$ploidy == 4,"pop"])

# before this you must switch on the return function in the plotsfs function
dip.sfs<-bind_rows(lapply(pop20.dip,plotsfs))
tet.sfs<-bind_rows(lapply(pop20.tet,plotsfs))

dip.sfs$freq<-as.factor((dip.sfs$freq))
tet.sfs$freq<-as.factor((tet.sfs$freq))

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

svg("diploid.sfs.blue.points.svg",width=10)
ggplot(dip.sfs,aes(x=freq,y=prop,fill=cat))+
  geom_boxplot()+
  #geom_point(pch=21,size=2,alpha=0.5,position = position_jitterdodge())+
  ggtitle("diploid")+
  ylab("Proportion of variants")+
  xlab("Number of alleles")+
  scale_color_manual(values=c("#1c96edff","#0511efff"))+
  scale_fill_manual(values=c("#1c96edff","#0511efff"))+
  scale_y_continuous(limits = c(0,0.45))+
  geom_point(pch=21,size=1,alpha=0.3,position = position_jitterdodge())+
  th
dev.off()

svg("tetraploid.sfs.orange.points.svg",width=10)
ggplot(tet.sfs,aes(x=freq,y=prop,fill=cat))+
  ggtitle("tetraploid")+
  ylab("Proportion of variants")+
  xlab("Number of alleles")+
  scale_color_manual(values=c("#ffa500","#ff3302ff"))+
  scale_fill_manual(values=c("#ffa500","#ff3302ff"))+
  scale_y_continuous(limits = c(0,0.45))+
  geom_boxplot()+
  geom_point(pch=21,size=1,alpha=0.3,position = position_jitterdodge())+
  #geom_bar(position = position_jitterdodge(),stat = "identity",alpha=0.3,width=0.1)+
  th
dev.off()
pop="BAB"

# plot together diploid and tet fourfold and tetraploid sfs
head(dip.sfs)
dim(tet.sfs)
dim(dip.sfs)
# add ploidy category
sfs.all<-rbind(tet.sfs,dip.sfs)
ploidy<-gen.indpop20[,c("pop","ploidy")]
sfs.all<-merge(sfs.all,ploidy)

svg("fourfold.sfs.points.svg",width=10)
ggplot(sfs.all[sfs.all$cat == "fourfold",],aes(x=freq,y=prop,fill=ploidy))+
  ggtitle("Fourfold")+
  ylab("Proportion of variants")+
  xlab("Allele frequency")+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  scale_y_continuous(limits = c(0,0.45))+
  geom_boxplot()+
  geom_point(pch=21,size=1,alpha=0.3,position = position_jitterdodge())+
  #geom_bar(position = position_jitterdodge(),stat = "identity",alpha=0.3,width=0.1)+
  th
dev.off()

svg("zerofold.sfs.points.svg",width=10)
ggplot(sfs.all[sfs.all$cat == "zerofold",],aes(x=freq,y=prop,fill=ploidy))+
  ggtitle("Zerofold")+
  ylab("Proportion of variants")+
  xlab("Allele frequency")+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  scale_y_continuous(limits = c(0,0.45))+
  geom_boxplot()+
  geom_point(pch=21,size=1,alpha=0.3,position = position_jitterdodge())+
  #geom_bar(position = position_jitterdodge(),stat = "identity",alpha=0.3,width=0.1)+
  th
dev.off()

# calculate pi from sfs
calcpi<- function (pop) {
  sfsHIGH<- read.table(paste("zerofold/",pop,".zerofold.sfs",sep=""),sep=",")
  sfsHIGH<- apply(sfsHIGH,2,sum) # sum the frequencies across genome
  sfsWEAK<- read.table(paste("fourfold/",pop,".fourfold.sfs",sep=""),sep=",")
  sfsWEAK<- apply(sfsWEAK,2,sum) # sum the frequencies across genome
  piHIGH<-est_pi(sfsHIGH)
  piWEAK<-est_pi(sfsWEAK)
  pi<-as.data.frame(cbind(pop,piHIGH,piWEAK))
  return(pi)
}

pi.pop20<-bind_rows(lapply(pop20,calcpi))
pi.pop20<-merge(pi.pop20,distinct(indpop[,c("pop","ploidy")]))
pi.pop20$piWEAK<-as.numeric(pi.pop20$piWEAK)
pi.pop20$piHIGH<-as.numeric(pi.pop20$piHIGH)
names(pi.pop20)


svg("fourfold.div.16.chrom.svg")
ggplot(pi.pop20,aes(x=ploidy,y=piWEAK,fill=ploidy))+
geom_boxplot()+
  ylab("\u03C0 (x 0.01)")+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,size=3,position=position_jitterdodge(),alpha=0.8)+
  geom_text(aes(label=pop),size=4,check_overlap = T,nudge_x = 0.2)+
  th+ theme (legend.position = "none")
dev.off()

g<-glm(piWEAK~ploidy,pi.pop20,family="gaussian")
summary(g)
wilcox.test(pi.pop20[pi.pop20$ploidy == "2","piWEAK"],pi.pop20[pi.pop20$ploidy == "4","piWEAK"])

ggplot(pi.pop20,aes(x=ploidy,y=piHIGH,fill=ploidy))+
  geom_boxplot()+
  ylab("\u03C0 (x 0.01)")+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(pch=21,size=3,position=position_jitterdodge(),alpha=0.8)+
  geom_text(aes(label=pop),size=4,check_overlap = T,nudge_x = 0.2)+
  th+ theme (legend.position = "none")

