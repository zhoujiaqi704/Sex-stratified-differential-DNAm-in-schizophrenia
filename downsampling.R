##bootstrap 4000
set.seed(2402) # Setting the seed for replication purposes

n <- 24 # Sample size (male control)

P <- 4000 # Number of bootstrap samples 

variable <- pd_m_control$Sample_Name


# initialize a matrix to store the permutation data
PermSamples <- matrix(0, nrow=n, ncol=P)
# each column is a permutation sample of data

# now, get those permutation samples, using a loop
for(i in 1:P){
  PermSamples[,i] <- sample(variable, size= n, replace=FALSE)
}

# we can take a quick look at the first 5 columns of PermSamples
PermSamples[, 1:5]

library(dplyr)
case=as.data.frame(pd_m_case$Sample_Name)
colnames(case)[1]="Name"

case=bind_cols(replicate(4000, case, simplify = FALSE))
colnames(PermSamples)=colnames(case)
PermSamples.all=rbind(PermSamples,case)

#group info
#do parallel
library(doParallel);
options(stringsAsFactors = F)
cl <- makeCluster(4)
registerDoParallel(cl)

group<-factor(c(rep("Contrl",24),rep("Case",41)))

Perm.dataMethyl = vector(mode="list",length = 4000)

for (i in 1:ncol(PermSamples.all)){
  Perm.dataMethyl$beta[[i]] <- RegressSCZ_Male[,match(PermSamples.all[,i],colnames(RegressSCZ_Male))]
  Perm.dataMethyl$pd[[i]] <- data.frame(PermSamples.all[,i],group)
  }


#DMP
library(limma)
dmpTableSz=as.data.frame(matrix(NA,nrow=250028, ncol=4000),row.names=rownames(RegressSCZ_Male))
for (i in 1:length(Perm.dataMethyl$beta)) {
design= model.matrix(~group, data=Perm.dataMethyl$pd[[i]])#desigh matrix for ref-based
fit=lmFit(Perm.dataMethyl$beta[[i]],design)
fit = eBayes(fit)
pval = fit$p.value[,2]
dmpTableSz[,i] = pval
}

dmpTableSz_male.resampling4000=dmpTableSz
colnames(dmpTableSz_male.resampling4000)=paste("PVal",1:4000,sep="")
save(dmpTableSz_male.resampling4000,file="dmpTableSz_male.resampling4000.rda")


#summarize
PermP.male.resampling.num4000=as.data.frame(matrix(NA,nrow=4000,ncol=3))
for(i in 1 : ncol(dmpTableSz_male.resampling4000)) {
	PermP.male.resampling.num4000[i,1] <- sum(as.numeric(dmpTableSz_male.resampling4000[,i]<0.05))
	PermP.male.resampling.num4000[i,2] <- sum(as.numeric(dmpTableSz_male.resampling4000[,i]<1e-04))
	PermP.male.resampling.num4000[i,3] <- sum(as.numeric(dmpTableSz_male.resampling4000[,i]<1e-05))
}

colnames(PermP.male.resampling.num4000)=c("num.adjP","num.4e","num.5e")
rownames(PermP.male.resampling.num4000)=paste("subsampled", 1:4000, sep = "")

#
dmpTableSz_male.resampling.fdr4000=as.data.frame(matrix(NA,nrow=250028,ncol=4000))
for(i in 1: ncol(dmpTableSz_male.resampling4000)) {
	dmpTableSz_male.resampling.fdr4000[,i] <- p.adjust(dmpTableSz_male.resampling4000[,i],"BH")
}

rownames(dmpTableSz_male.resampling.fdr4000)=rownames(dmpTableSz_male.resampling4000)
colnames(dmpTableSz_male.resampling.fdr4000)=paste("fdr",1:4000,sep="")


PermP.male.resampling.num.fdr4000=as.data.frame(matrix(NA,nrow=4000,ncol=1))
for(i in 1 : ncol(dmpTableSz_male.resampling.fdr4000)) {
	PermP.male.resampling.num.fdr4000[i,1] <- sum(as.numeric(dmpTableSz_male.resampling.fdr4000[,i]<0.05))
}

colnames(PermP.male.resampling.num.fdr4000)=c("adjP")
PermP.male.resampling.num4000=cbind(PermP.male.resampling.num4000,PermP.male.resampling.num.fdr4000)
save(PermP.male.resampling.num4000,file="PermP.male.resampling.num4000.rda")


###plot
library("ggplot2")
a <- ggplot(PermP.male.resampling.num4000, aes(x = adjP))+ geom_histogram(color = "black", fill = "gray",bins=20) +
  geom_vline(aes(xintercept = 3),
             linetype = "dashed", size = 1,color="red")+
             labs(x="number of mDMPs (adj.p<0.05)",y="frequency",title="4000 times subsampling (males)")+
             theme(plot.title = element_text(hjust = 0.5),text=element_text(size=15))
ggsave("subsampling4000.comparison.adjP.pdf")


