args <- commandArgs(TRUE)

inputdata <- args[1] # first argument is the depth file
pop <- args[2] # second argument is pop
scf <- args[3]
inputdata
a<-read.table(inputdata,header=F,sep=" ",col.names=c("freq","dp","rr"))

a$dp<-as.numeric(as.character(a$dp))

summary(a)


png(paste(pop,scf,".dp.png",sep=""))
plot(freq~dp,a, main=pop)
dev.off()


png(paste(pop,scf,".dp.zoom.png",sep=""))
plot(freq~dp,a[a$dp<50,], main=pop)
abline(v=8,col="red")
dev.off()

