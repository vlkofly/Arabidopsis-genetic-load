##############################
# environment setup 
##############################
# change your path to your working dir with the vcf file, load in libraries

library(ggplot2)
library(adegenet)
library(adegraphics)
library(vcfR)
library(pegas)
library(StAMPP)
library(gplots)
library(ape)

##############################
### IMPORT SNP data from VCF 
##############################
# prepare the data
# download from /storage/brno12-cerit/home/filip_kolar/brno3/ERC/arenosa_mapped_to_arenosa/filtering_results/masked
#' zgrep CHROM scaffold_1arenosa_snp.fourfold.biallelic.dp8.nc.m0.5.vcf.gz | tr \\t \\n | grep _ | grep -f pop20.txt | wc -l > pop20.samples.txt
#' 238 samples
#' bcftools view -S pop20.samples.txt scaffold_1arenosa_snp.fourfold.biallelic.dp8.nc.m0.5.vcf.gz | bcftools view --threads 2 -Oz -c 1 > scf1.pop20.vcf.gz
#' bcftools concat --threads 2 -Oz scf1.pop20.vcf.gz scf2.pop20.vcf.gz > scf1.scf2.pop20.vcf.gz get a bit more sites
#' size reduced three times
setwd("/home/vlkofly/Arabidopsis_science/arenosa/pca")
getwd()

vcf <- read.vcfR("scf1.scf2.pop20.vcf.gz",nrows = 100000)
head(vcf)               #check
vcf@fix[1:10,1:5]      #check

### convert to genlight
aa.genlight.all <- vcfR2genlight.tetra(vcf)
locNames(aa.genlight.all) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")   # add real SNP.names
pop(aa.genlight.all)<-substr(indNames(aa.genlight.all),1,3)               # add pop names: here pop names are first 5 chars of ind name
ls()
#rm("vcf") # remove the vcf to save space
toremove<-c("VLH_01dl","STR_01dl","PEK_04tl","PIZ_02dl",toremove) # depth


aa.genlight <- aa.genlight.all
rm("aa.genlight.all")
(indNames(aa.genlight.all)  %in% toremove)
aa.genlight <- new("genlight", (as.matrix(aa.genlight)[, colSums(as.matrix(aa.genlight) != 0, na.rm = T)>0]))   #remove all nonvar sites
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3)

#get only DIPLOIDS or TETRAPLOIDS
aa.genlight.badremove <- aa.genlight # save original for subsetting
include.list<-grep("4",ploidy(aa.genlight.badremove)) # select here
aa.genlight <- new("genlight", (as.matrix(aa.genlight.badremove)[include.list, ])) # repeat this
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3)

ploidy.pop<-as.data.frame(cbind(substr(indNames(aa.genlight),1,3),ploidy(aa.genlight)))
summary(ploidy.pop)
ploidy.pop$V2<-as.numeric(as.character(ploidy.pop$V2))
ploidy.pop<-ploidy.pop[ploidy.pop$V2 %in% c("2","4"),]


# get only: 
toMatch <- c("") #crown DIPLOIDS

exclude.list<-grep(paste(toMatch,collapse="|"),pop(aa.genlight))                

aa.genlight <- new("genlight", (as.matrix(aa.genlight)[-(exclude.list), ])) # if excluding
#aa.genlight <- new("genlight", (as.matrix(aa.genlight)[(exclude.list), ]))# if inculding selected populations
aa.genlight <- new("genlight", (as.matrix(aa.genlight)[, colSums(as.matrix(aa.genlight) != 0, na.rm = T)>0]))   #remove all nonvar sites
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3)

# get only good populations

setwd("goodpop/") # change folder for a given analysis

#check
aa.genlight
indNames(aa.genlight)
grep("MIL",indNames(aa.genlight))
ploidy(aa.genlight)
#as.matrix(aa.genlight)[1:16,1:10]
# check number of samples per population
pop(aa.genlight)
popsize<-as.data.frame(table(pop(aa.genlight)))
popsize[order(popsize$Freq),]
write.table(popsize,"popsize.tsv",sep="\t",quote=F)

##############################
#  data checks & statistics
##############################


###set plotting environment
th<-theme(axis.title = element_text(size=23),
          axis.text=element_text(size=23,color="black"),
          plot.title=element_text(size=25),
          legend.text = element_text(size=18),
          strip.text.x = element_text(size = 18),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black",linetype="solid"))

## 1) plot the entire matrix
png("glplot.scf1.scf2.png",width = 1000, height = 500, units = "px")
glPlot (aa.genlight)  # takes some time
dev.off()



# I would remove some of the shitty individuals
### 2) N missing SNPs per sample # run it again for all godind dataset
x <- summary(t(as.matrix(aa.genlight)))
write.table(substr(x[7,],9,100), file = "missing.persample.txt", sep = "\t", quote=F, col.names=F)  # NAs, if present, are in seventh row of summary

### 3) check persample read depth and plot it
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
sample.dp<-as.data.frame(colMeans(dp,na.rm=T))
sample.dp$sample.name<-row.names(sample.dp)
dim(sample.dp)
write.table(sample.dp,"sample.depth.snps.tsv",sep="\t",row.names=T,col.names=F)
mean(colMeans(dp,na.rm=T)) #average coverage 
head(dp)
png(file="DP.png",width=1000, height=500)
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)", las=2, cex=0.7)
dev.off()



### get the list of samples with low coverage
names(sample.dp)[1]<-"dp"
head(sample.dp)
remove.sample<-sample.dp[sample.dp$dp < 14.5,"sample.name"]
summary(sample.dp)
cat(remove.sample)
rem.ind<-which((indNames(aa.genlight)  %in% remove.sample),arr.ind=TRUE)
aa.genlight <- new("genlight", (as.matrix(aa.genlight)[-(rem.ind), ]))
which((indNames(aa.genlight)  == "FEL_09ta"),arr.ind=TRUE)
aa.genlight <- new("genlight", (as.matrix(aa.genlight)[-(69), ]))
aa.genlight <- new("genlight", (as.matrix(aa.genlight)[, colSums(as.matrix(aa.genlight) != 0, na.rm = T)>0]))   #remove all nonvar sites
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3)

png("hist.sample.depth.zoom.png")
ggplot(sample.dp[sample.dp$dp < 20,],aes(x=dp))+
  geom_histogram(bins=40)+
  th
dev.off()

sample.dp[sample.dp$dp < 10 ,"sample"] # MIL_05 is the only low coverage individual
##rather plot histogram
dp.sub<-as.data.frame(dp[dp < 100])
names(dp.sub)<-"depth"
png("hist.depth.png")
ggplot(dp.sub,aes(x=depth))+
  geom_histogram()+
  th
dev.off()


png(file="DP_detail.png",width=1000, height=500)
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)", las=2, cex=0.7, ylim=c(0,50))
abline(h=4, col="red")
dev.off()

### 4) Calculate missingness percentage per sample
m<-as.data.frame(cbind("m" = as.numeric(substr(x[7,],9,100)),"sample.name" = colnames(x)))
m$sample.name<-gsub(" ","",as.character(m$sample.name))
m$m<-as.numeric(as.character(m$m))
head(m)
head(sample.dp)
head(sample.dp$sample.name)
length(m$m)
sample.stats<-merge(sample.dp,m,by = "sample.name")
names(sample.stats)[2]<-"dp"

sample.stats$misperc<-as.numeric(as.character(sample.stats$m))/100000
#sample.stats$depth<-as.numeric(as.character(sample.stats$`colMeans(dp, na.rm = T)`))
summary(sample.stats)
head(sample.stats)
# remove the removed individuals from sample.stats and get values of coverage that will go to the paper
final.sample.stats<-sample.stats[! sample.stats$sample.name %in% c(remove.sample,"FEL_09ta"),]
summary(final.sample.stats)
head(final.sample.stats)
# ok I want depth per ploidy, add ploidy
ploidy(aa.genlight)-> sample.ploidy
final.sample.stats<-cbind(final.sample.stats,sample.ploidy)
final.sample.stats$sample.ploidy<-as.factor(as.character(final.sample.stats$sample.ploidy))
library(tidyverse)
final.sample.stats %>% select(sample.ploidy,dp) %>% group_by (sample.ploidy) %>% summarise (meandp = mean(dp))

# is the difference in depth significant 

wilcox.test(final.sample.stats[final.sample.stats$sample.ploidy ==2,"dp"], final.sample.stats[final.sample.stats$sample.ploidy ==4,"dp"] )
# yes it is

svg("dp.sample.ploidy.svg",width=5)
ggplot(final.sample.stats,aes(x=sample.ploidy,y=dp,color=sample.ploidy))+
  geom_violin(draw_quantiles = 0.5)+
  scale_fill_manual(values=c("#1e90ff","#ffa500"))+
  scale_color_manual(values=c("#1e90ff","#ffa500"))+
  geom_point(aes(fill=sample.ploidy),pch=21,size=4,alpha=0.5,position=position_jitterdodge())+
  ylab("Depth of coverage")+
  th+
  theme(axis.text.x=element_text(size=23,color="black"),legend.text = element_blank(),
        legend.position = "none")
dev.off()

head(sample.stats[order(sample.stats$dp),])

png("hist.sample.depth.zoom.png")
ggplot(sample.dp[sample.dp$dp < 15,],aes(x=dp))+
  geom_histogram(bins=40)+
  th
dev.off()
# missingness and depth per population
sample.stats$pop<-substr(sample.stats$sample,1,3)
#merge populations
pop.stats<-aggregate(sample.stats[,c("dp","misperc")],by=list(sample.stats$pop),FUN=mean)
write.table(pop.stats,"pop.stats.tsv",sep="\t",row.names=F)
svg("depth_pop.svg",width=15)
ggplot(sample.stats,aes(x=pop,y=dp))+
  geom_boxplot()+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size=30),
        axis.line=element_line(color='black'),
        plot.title=element_text(size=17.5),
        axis.text.x=element_text(color='black', size=20, hjust=0.5, vjust=0.5,angle=90),
        axis.text.y=element_text(color='black', size=30, hjust=1, vjust=1),
        #text=element_text(family="serif"),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_blank())
dev.off()


### 5) plot total allele frequency spectrum (AFS) of entire dataset
mySum <- glSum(aa.genlight, alleleAsUnit = TRUE) 
head(mySum)
png("AFS.png",width = 1000, height = 500, units = "px")
barplot(table(mySum), col="blue", space=0,xlim=c(0,80), xlab="Allele counts",  main="Distribution of ALT allele counts in total dataset")
dev.off()

png("AFS.total.png",width = 1000, height = 500, units = "px")
barplot(table(mySum), col="blue", space=0, xlab="Allele counts",  main="Distribution of ALT allele counts in total dataset")
dev.off()

############
#   PCA 
############
# do PCA, retain first 300 axes (for later use in find.clusters)
pca.1 <- glPcaFast(aa.genlight, nf=300) # see the modified function glPcaFast


# proportion of explained variance by first three axes
round((pca.1$eig[1]/sum(pca.1$eig)*100),1) # proportion of variation explained by 1st axis
round((pca.1$eig[2]/sum(pca.1$eig)*100),1) # proportion of variation explained by 2nd axis
round((pca.1$eig[3]/sum(pca.1$eig)*100),1) # proportion of variation explained by 3rd axis


# just to see pops coloured in a palette
# how many populations? 
pop(aa.genlight)
col <- funky(length(levels(pop(aa.genlight))))
s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, col=transp(col,.6), 
        ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F)

# save nice figs
#pdf ("PCA_SNPs_ax12.pdf", width=14, height=7)
pdf ("PCA_SNPs_ax12.pdf", width=14, height=7)
g1 <- s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, xlab=paste(round((pca.1$eig[1]/sum(pca.1$eig)*100),1)," %"), ylab=paste(round((pca.1$eig[2]/sum(pca.1$eig)*100),1)," %"), col=transp(col,.6), 
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plot = FALSE)
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels = list(box = list(draw = FALSE), 
               optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

#pdf ("PCA_all_SNPs_ax13.pdf", width=14, height=7)
pdf ("PCA_SNPs_ax13.pdf", width=14, height=7)
g1 <- s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=3,xlab=paste(round((pca.1$eig[1]/sum(pca.1$eig)*100),1)," %"), ylab=paste(round((pca.1$eig[3]/sum(pca.1$eig)*100),1)," %"),  col=transp(col,.6), 
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plot = FALSE)
g2 <- s.label (pca.1$scores, xax=1, yax=3, ppoints.col = "red", plabels = list(box = list(draw = FALSE), 
                                                                               optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

#ploidy - differentiated plots
pdf ("PCA_all_ploidycol_SNPs_ax12.pdf", width=14, height=7)
g1 <- s.class(pca.1$scores, as.factor(as.vector(ploidy(aa.genlight))), xax=1, yax=2, xlab=paste(round((pca.1$eig[1]/sum(pca.1$eig)*100),1)," %"), ylab=paste(round((pca.1$eig[2]/sum(pca.1$eig)*100),1)," %"), col=transp(c("#FF0000","#0000FF")), 
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plab.cex = 0 , plot = FALSE)
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels = list(box = list(draw = FALSE), 
              optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

pdf ("PCA_all_ploidycol_SNPs_ax13.pdf", width=14, height=7)
g1 <- s.class(pca.1$scores, as.factor(as.vector(ploidy(aa.genlight))), xax=1, yax=3, xlab=paste(round((pca.1$eig[1]/sum(pca.1$eig)*100),1)," %"), ylab=paste(round((pca.1$eig[3]/sum(pca.1$eig)*100),1)," %"), col=transp(c("#FF0000","#0000FF")), 
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plab.cex = 0 , plot = FALSE)
g2 <- s.label (pca.1$scores, xax=1, yax=3, ppoints.col = "red", plabels = list(box = list(draw = FALSE), 
                                                                               optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()


# now create the final figure for my paper, I would rather plot it with ggplot, but lets try to use what I have and fix in the inkscape
# for the labels I need just one individual per population
fip<-grep("01",row.names(pca.1$scores))
length(fip)
pca.1$scores[fip,]

svg("pca.1.2.svg")
s.class(pca.1$scores, as.factor(as.vector(ploidy(aa.genlight))), xax=1, yax=2, xlab=paste(round((pca.1$eig[1]/sum(pca.1$eig)*100),1)," %"), ylab=paste(round((pca.1$eig[2]/sum(pca.1$eig)*100),1)," %"), col=transp(c("#1e90ff","#ffa500")), 
        ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=F, pgrid.draw =F, plab.cex = 0 , plot = TRUE)
dev.off()

svg("pca.1.2.labels.svg")
s.label (pca.1$scores[fip,], xax=1, yax=2,labels=levels(pop(aa.genlight)), plabels = list(box=list(draw=FALSE),optim=TRUE),  paxes.draw=F, pgrid.draw =F, plot = T)
dev.off()

svg("pca.1.3.svg")
s.class(pca.1$scores, as.factor(as.vector(ploidy(aa.genlight))), xax=1, yax=3, xlab=paste(round((pca.1$eig[1]/sum(pca.1$eig)*100),1)," %"), ylab=paste(round((pca.1$eig[2]/sum(pca.1$eig)*100),1)," %"), col=transp(c("#1e90ff","#ffa500")), 
        ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=F, pgrid.draw =F, plab.cex = 0 , plot = TRUE)
dev.off()

svg("pca.1.3.labels.svg")
s.label (pca.1$scores[fip,], xax=1, yax=3,labels=levels(pop(aa.genlight)), plabels = list(box=list(draw=FALSE),optim=TRUE),  paxes.draw=F, pgrid.draw =F, plot = T)
dev.off()

?s.class
  ##################################
# distance matrix for SplitsTree
##################################
# calculate matrix Nei's 1972 distance between indivs and pops
# export matrix - open it then in SplitsTree

aa.D.ind <- stamppNeisD(aa.genlight, pop = FALSE)                 # Nei's 1972 distance between indivs
stamppPhylip(aa.D.ind, file="al.indiv_Neis_distance.phy.dst")    

aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)                  # Nei's 1972 distance between pops
stamppPhylip(aa.D.pop, file="al.pops_Neis_distance.phy.dst")      


############################################
# additional analyses using these distances
############################################

### create a dist object
colnames(aa.D.ind) <- rownames(aa.D.ind) 
aa.D.ind.dist<-as.dist(aa.D.ind, diag=T)
attr(aa.D.ind.dist, "Labels")<-rownames(aa.D.ind)          # name the rows of a matrix  

### plot NJ tree
plot(nj(aa.D.ind))
write.tree(njs(aa.D.ind),file="NJ.stampp.tree.tre")
summary(aa.D.ind)
### heatmap of the pop distance matrix
pdf ("heatmap_Neis_Distances.pdf", width=14, height=7)
 heatmap.2(aa.D.pop, trace="none", cexRow=0.7, cexCol=0.7)
dev.off()

### calculate AMOVA (differentiation among populations)
pops<-as.factor(substr(rownames(aa.D.ind),1,3))
(res <- pegas::amova(aa.D.ind.dist ~ pops))



##################
#     DAPC
##################
## find clusters = K-means clustering
grp <- find.clusters(aa.genlight, max.n.clust=20, glPca = pca.1, perc.pca = 99, n.iter=1e6, n.start=1000)

# good to add more random starts (e.g. n.start=1000) to reach convergence in K-means)
# SELECT 100% of Pcs to retain, i.e., retaining all the variation (i.e. clasical K-means without PCA transformation)
# if set "perc.pca", just enter after the first question and wait
# select optimal number of K based on BIC (optimally where is the bend ... i.e., the line again starts to rise)

# optionally save the grouping
write.table(grp$grp, file="grouping_diploid_K7.txt", sep="\t", quote=F, col.names=F)

## run DAPC
dapc1<-dapc(aa.genlight, grp$grp, glPca = pca.1)
# two interactive questions 
#   1) how many PCs - aim is retaining a few PCs without sacrifying too much information, e.g. 50
#      - the more axes sellected the better the groups will be resolved (the less "admixture" will be in the compoplot)
#   2) how many discriminant functions retain (for small K, all could be retained)

## plot results
scatter(dapc1)           # scatterplot

# barplot - cluster memberships of individuals
#pdf("DAPC_all.pdf",width=20,height=5)
col <- funky(7)
compoplot(dapc1, cex.names = 0.4, legend=F, col=col)
#compoplot(dapc1, cex.names = 0.4, legend=F, col=c("#008B8B", "#8B008B", "#00008B", "#FF0000", "#00EEEE","#0000CD","#8B2323"))  # BSW K=7
#dev.off()


###################################################
####   distance-based analyses  #######
###################################################



# OPTIONAL read back the saved phylip matrices from .dst files
#aa.D.ind.imported <- as.matrix(read.table("aa.indiv_Neis_distance.phy.dst", skip=1, row.names=1))
#colnames(aa.D.ind.imported) <- rownames(aa.D.ind.imported) 
#aa.D.ind.noHARnoDFS <- aa.D.ind.imported[!rownames(aa.D.ind.imported) %in% c("DFS_001_1","HAR_005_1","HAR_006_1"), !colnames(aa.D.ind.imported) %in% c("DFS_001_1","HAR_005_1","HAR_006_1")]         

#aa.D.pop.imported <- as.matrix(read.table("aa.pops_Neis_distance.phy.dst", skip=1, row.names=1))
#colnames(aa.D.pop.imported) <- rownames(aa.D.pop.imported) 
#aa.D.pop.noHARnoDFS <- aa.D.pop.imported[!rownames(aa.D.pop.imported) %in% c("DFS","HAR"), !colnames(aa.D.pop.imported) %in% c("DFS","HAR")]         


### create the dist objects used in analyses below
colnames(aa.D.ind) <- rownames(aa.D.ind) 
#aa.D.ind.noHARnoDFS <- aa.D.ind[!rownames(aa.D.ind) %in% c("DFS_001_1","HAR_005_1","HAR_006_1"), !colnames(aa.D.ind) %in% c("DFS_001_1","HAR_005_1","HAR_006_1")]    # if needed to exclude     
aa.D.ind.dist<-as.dist(aa.D.ind, diag=T)
attr(aa.D.ind.dist, "Labels")<-rownames(aa.D.ind)          # name the rows of a matrix  
pops<-as.factor(substr(rownames(aa.D.ind),1,3))

colnames(aa.D.pop) <- rownames(aa.D.pop) 
#aa.D.pop.noHARnoDFS <- aa.D.pop[!rownames(aa.D.pop) %in% c("DFS","HAR"), !colnames(aa.D.pop) %in% c("DFS","HAR")]         
aa.D.pop.dist<-as.dist(aa.D.pop, diag=T)
attr(aa.D.pop.dist, "Labels")<-rownames(aa.D.pop)          # name the rows of a matrix  



# trees, heatmaps
# ----------------
# plot NJ tree
library(ape)
plot(nj(aa.D.ind))
write.tree(nj(aa.D.ind),file="NJ.stampp.tree.tre")

### hemap of the pop distance matrix
library(gplots)
heatmap.2(aa.D.pop, trace="none", cexRow=0.7, cexCol=0.7)


###plot the distribution of distances
# -----------------------------------
aa.dist.df <- data.frame(t(combn(labels(aa.D.ind.dist),2)), as.numeric(aa.D.ind.dist))
hist(aa.dist.df$as.numeric.aa.D.ind.dist.,breaks=30, col="gray", xlab="Nei's distance")

aa.dist.df.pop <- data.frame(t(combn(labels(aa.D.pop.dist),2)), as.numeric(aa.D.pop.dist))
hist(aa.dist.df.pop$as.numeric.aa.D.pop.dist., breaks=30, col="gray", xlab="Nei's distance")

# what are the most differentiated pairs?)
tail(aa.dist.df.pop[with(aa.dist.df.pop, order(aa.dist.df.pop$as.numeric.aa.D.pop.dist.)), ], n=20)

### PCoA using the Stampp distance
# --------------------------------
pcoa.1 <- dudi.pco (aa.D.ind.dist, scannf=F, nf=4)      # do PCA, retain first five axes
pcoa.1$eig[1]/sum(pcoa.1$eig) # proportion of variation explained by 1st axis
pcoa.1$eig[2]/sum(pcoa.1$eig) # proportion of variation explained by 2nd axis
pcoa.1$eig[3]/sum(pcoa.1$eig) # proportion of variation explained by 3rd axis

###  plotting
pdf ("all_ax12_PCoA.pdf", width=14, height=7)
col <- funky(10)
g1 <- s.class(pcoa.1$li, pop(aa.genlight),  xax=1, yax=2, col=transp(col,.6), ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plot = FALSE)
g2 <- s.label (pcoa.1$li, xax=1, yax=2, ppoints.col = "red", plabels = list(box = list(draw = FALSE), optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

#check the PCoA on the same distances in vegan
library(vegan)
ordiplot(cmdscale(aa.D.ind.dist), type="text")
ordiplot(cmdscale(aa.D.ind.dist), type="points")


# Isolation by distance 
#-----------------------
#coords file = three-column file with column names of coords belonging to each pop (in the same order as pops !!)
#pop  lat	lon
#BEL	46.16167	16.115
#BGS	47.62806	13.00167
#BIH	44.88181	15.89882
# ...
coords <- read.csv ("pop_coords_all.txt", sep ="\t")     # tab-separated file for all pops
xy.coords.only<- subset(coords, select=c("lat","lon"))
Dgeo <- dist(xy.coords.only)

#check plotting the points on a map
library(ggmap)
map <- get_map(location =  c(lon = 14, lat = 50), zoom = 5)
#mapPoints <- ggmap(map) + geom_point(aes(x = coords$lon, y = coords$lat, size = myValues2$N.SNPs, shape=21), colour = "black", fill = "red", data = myValues2) + scale_shape_identity()
mapPoints <- ggmap(map) + geom_point(data = coords, aes(x = lon, y = lat, colour="blue")) + geom_text(data = coords, aes(x = lon, y = lat, label = pop, colour = "red"), 
          size = 4, vjust = 0, hjust = -0.5)
mapPoints 

#test IBD
library(ade4)
IBD <- mantel.randtest(Dgeo,aa.D.pop.dist)
IBD
plot(Dgeo,aa.D.pop.dist, pch=20,cex=.5)
abline(lm(aa.D.pop.dist~Dgeo))

#plot and check for denser areas in the plot indicating sub-groups
library(MASS)
dens <- kde2d(Dgeo,aa.D.pop.dist, n=300, lims=c(-1, 22, 0, 0.04))
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, aa.D.pop.dist, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(aa.D.pop.dist~Dgeo))
title("Correlation of Genetic and Geograpgical distances")


### AMOVA
#-------
pops<-as.factor(substr(rownames(aa.D.ind),1,5))
# grouping file = one-column list of group names belonging to each individual (in the same orderas indivs  !!)
groups<-read.csv("all_grouping.txt", header=F)$V1   

# one-level AMOVA
(res <- pegas::amova(aa.D.ind.dist ~ pops))  # default nperm=1000, pegas package must be loaded

# hierarchical AMOVA 
(res <- pegas::amova(aa.D.ind.dist ~ groups/pops)) # hierarchical AMOVA








###################################
# !!! HIGHLY EXPERIMENTAL BELOW
##################################


####################
#   plotting AFS
####################

library(lattice)
###plot total AFS
mySum <- glSum(aa.genlight, alleleAsUnit = TRUE) 
barplot(table(mySum), col="blue", space=0, xlab="Allele counts",  main="Distribution of ALT allele counts in total dataset")

###plot AFS per one pop - !!! after seppop you must remove the nonvariant positions
#aa.genlight.sep <- seppop(aa.genlight, drop=TRUE)  # separate genlight per population
#aa.genlight.AA007 <- new("genlight", (as.matrix(aa.genlight.sep$AA007))[,colSums(as.matrix(aa.genlight.sep$AA007)) > 0]) # removing the nonvariant positions
#summary(colSums(as.matrix(aa.genlight.AA007)))  #  check if there are no zeros
# plot unfolded AFS - for one pop.
#mySum <- glSum(aa.genlight.AA007, alleleAsUnit = TRUE) # alleleAsUnit = TRUE: eg tetraploid genotype equals four samples, diploid two, ...
#barplot(table(mySum), col="blue", space=0, xlab="Allele counts",  main="Distribution of ALT allele counts in AA007")    # plot the original counts of each category

### plot AFS per pop in a batch
aa.genlight.sep <- seppop(aa.genlight, drop=TRUE)  # separate genlight per population
aa.genlight.sep.2 <- lapply (aa.genlight.sep, function (pop) {new("genlight", (as.matrix(pop))[,colSums(as.matrix(pop)) > 0])})  # removing the nonvariant positions
listnames<-names(aa.genlight.sep.2)
for (i in seq(listnames)) {pop(aa.genlight.sep.2[[i]])<-substr(indNames(aa.genlight.sep.2[[i]]),1,5)} ##add pop identity to list elements
aa.genlight.sep.2[which(names(aa.genlight.sep.2) %in% c("AA013"))] <- NULL          # remove AA013 pop which has only one indiv



# loop over each population in a list of populations and draw into one fig
# all pops in one multi-panel plot - barplots
pdf("./AFS_all_barplot.pdf", width=20, height=20)
par(mfrow=c(6,7),mar=c(2,2,2,0))
mySum <- lapply (aa.genlight.sep.2, function (pop) {
  barplot(table(glSum(pop, alleleAsUnit=T)), col="blue", space=0,
          xlab="Allele counts",  main=paste(levels(pop(pop)),sum(table(glSum(pop, alleleAsUnit=T))),"SNPs", sep=" "))
}) 
dev.off()

par(mfrow=c(1,1))


#########################
#   per-pop statistics   !!! highly experimental, must be double-checked !!!!
########################

# for one object
aa.genlight.sep$AA007
AA007<-as.matrix(aa.genlight.sep$AA007)

nchr <- nrow(AA007)*2  # for a DIPLOID pop.
pi.persite <- 2*(colSums(AA007)/nchr)*(1-(colSums(AA007)/nchr)*(nchr/(nchr-1)))   # pi = 2pq (?? why is there nchr/(nchr-1))
pi<-(mean(pi.persite))                                      
S<-ncol(AA007)  # N SNPs
tmp <- 1:(nchr - 1)
a1 <- sum(1/tmp)
a2 <- sum(1/tmp^2)
b1 <- (nchr + 1)/(3 * (nchr - 1))
b2 <- 2 * (nchr^2 + nchr + 3)/(9 * nchr * (nchr - 1))
c1 <- b1 - 1/a1
c2 <- b2 - (nchr + 2)/(a1 * nchr) + a2/a1^2
e1 <- c1/a1
e2 <- c2/(a1^2 + a2)
tajima.D <- (pi - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
Het.exp <- 2*(colSums(AA007)/nchr)*(1-(colSums(AA007)/nchr))

# loop over all pops
# I take the per-pop splitted genlights, PRIOR removing the within-pop nonvariant sites
aa.genlight.sep.calc <- aa.genlight.sep
aa.genlight.sep.calc[which(names(aa.genlight.sep.calc) %in% c("AA013"))] <- NULL          # remove AA013 pop which has only one indiv

myValues <- lapply (aa.genlight.sep.calc, function (pop) {
  print (levels(pop(pop)))  # print the pop ID being processed
  
  if (mean(ploidy(pop))==2) {
    nchr <- (nrow(as.matrix(pop))*2)  # for a DIPLOID pop.
  } else {nchr <- (nrow(as.matrix(pop))*4)}  # for a TETRAPLOID pop.
  pi.persite <- 2*(colSums(as.matrix(pop))/nchr)*(1-(colSums(as.matrix(pop))/nchr))*(nchr/(nchr-1))   # pi = 2pq
  pi<-(mean(pi.persite))                                      
  
  S<-ncol(as.matrix(pop))  # N SNPs
  tmp <- 1:(nchr - 1)
  a1 <- sum(1/tmp)
  a2 <- sum(1/tmp^2)
  b1 <- (nchr + 1)/(3 * (nchr - 1))
  b2 <- 2 * (nchr^2 + nchr + 3)/(9 * nchr * (nchr - 1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (nchr + 2)/(a1 * nchr) + a2/a1^2
  e1 <- c1/a1
  e2 <- c2/(a1^2 + a2)
  tajima.D <- (pi - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
  
  if (mean(ploidy(pop))==2) {
    nchr <- (nrow(as.matrix(pop))*2)  # for a DIPLOID pop.
    Het.exp.persite <- 2*(colSums(as.matrix(pop))/nchr)*(1-(colSums(as.matrix(pop))/nchr))
    Het.exp<-mean(Het.exp.persite)
  } else {nchr <- (nrow(as.matrix(pop))*4)
  Het.exp.persite <- 1-(colSums(as.matrix(pop))/nchr)^4-(1-(colSums(as.matrix(pop))/nchr))^4
  Het.exp<-mean(Het.exp.persite)
  }  # for a TETRAPLOID pop.
  
  N.SNPs <- sum(colSums(as.matrix(pop)) > 0)        # NSNPs = number of nonzero columns in a matrix
  N.SNPs.nosingl <- sum(colSums(as.matrix(pop)) > 1) # NSNPs excluding those with singletons = number of nonzero columns in a matrix
  
  return(c(pi,tajima.D,Het.exp,N.SNPs,N.SNPs.nosingl))
}) 


myValues2 <- t(as.data.frame(myValues))
colnames(myValues2) <- c("pi","Tajimas_D","Het.exp","N.SNPs","N.SNPs.nosingl")
myValues2 <- as.data.frame(myValues2)
write.table(myValues2,file="aa.popstats.txt", sep="\t")

# positive Tajima's D = excess of heterozygosity / freqs of polymorphic nucleotides are too nearly equal
#             - balancing selection, diversifying sel, recent admixture 
# negative Tajima's D = too many rare alleles / freqs of polymorphic nucleotides are too unequal
#             - selection against slightly deleterious alleles, pop. expansion

#plot this on map
coords <- read.csv ("pop_coords_noDFS.txt", sep ="\t")     # for all pops, without DFS

library(ggmap)
map <- get_map(location =  c(lon = 14, lat = 50), zoom = 5)
#mapPoints <- ggmap(map) + geom_point(aes(x = coords$lon, y = coords$lat, size = myValues2$N.SNPs, shape=21), colour = "black", fill = "red", data = myValues2) + scale_shape_identity()
mapPoints <- ggmap(map) + geom_point(aes(x = coords$lon, y = coords$lat, size = myValues2$pi, colour = factor(coords$ploidy)), data = myValues2) + scale_colour_manual(values = c("red","blue"))
mapPoints

pdf(file="./figs/map_pi.pdf")
mapPoints
dev.off()


#####################











##################### 
#SOME MODIFIED FUNCTIONS IF NEEDED

# a function for conversion from vcfR object to genlight in tetraploids
vcfR2genlight.tetra <- function (x, n.cores = 1) 
  {
    bi <- is.biallelic(x)
    if (sum(!bi) > 0) {
      msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
      msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
      msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
      warning(msg)
      x <- x[bi, ]
    }
    x <- addID(x)
    CHROM <- x@fix[, "CHROM"]
    POS <- x@fix[, "POS"]
    ID <- x@fix[, "ID"]
    x <- extract.gt(x)
    x[x == "0|0"] <- 0
    x[x == "0|1"] <- 1
    x[x == "1|0"] <- 1
    x[x == "1|1"] <- 2
    x[x == "0/0"] <- 0
    x[x == "0/1"] <- 1
    x[x == "1/0"] <- 1
    x[x == "1/1"] <- 2
    x[x == "1/1/1/1"] <- 4
    x[x == "0/1/1/1"] <- 3
    x[x == "0/0/1/1"] <- 2
    x[x == "0/0/0/1"] <- 1
    x[x == "0/0/0/0"] <- 0
    if (requireNamespace("adegenet")) {
      x <- new("genlight", t(x), n.cores = n.cores)
    }
    else {
      warning("adegenet not installed")
    }
    adegenet::chromosome(x) <- CHROM
    adegenet::position(x) <- POS
    adegenet::locNames(x) <- ID
    return(x)
  }

#######################


## -------------------------------------------------------
### a patch for MUCH MUCH faster PCA calculation on genlight objects
# see https://github.com/thibautjombart/adegenet/pull/150
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks..
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar
  # to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS
  # to support thousands of samples, this could be
  # replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}

# ---------------------------------------------------------



