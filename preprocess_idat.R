
library(ChAMP)
library(limma)
library(sva)
library(ggfortify)

myLoad <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") 

####normalization
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",method='BMIQ',cores=4)
champ.SVD(beta=myNorm,pd=myLoad$pd)

#high-quality
b<-read.csv("./high_qualityProbes.csv")
myNorm2<-myNorm[rownames(myNorm)%in%b[,1]==T,]

#pca
df<-t(myNorm)
jpeg(file = "pcaSCZ_Norm.jpg",width =2000,height = 2000,units = "px",res =300)
d<-autoplot(prcomp(df), data = myLoad$pd, colour ='Sample_Group',shape = "Sample_Plate")
plot(d)
dev.off()

###########batch and position correction
Combat1_position <- champ.runCombat(beta=myNorm2,pd=myLoad$pd,batchname=c("Array"))#position
Combat2_Batch <- champ.runCombat(beta=Combat1_position,pd=myLoad$pd,batchname=c("Sample_Plate"))#batch

champ.SVD(beta=Combat2_batch,pd=myLoad$pd)

#pca for combat
jpeg(file = "pcaSCZ_aftPosBat.jpg",width =2000,height = 2000,units = "px",res =300)
d<-autoplot(prcomp(t(Combat2_batch)), data = myLoad$pd, colour ='Sample_Group',shape = "Sample_Plate")
plot(d)
dev.off()

###cell type composion estimate (modified from champ.refbase)
source("./refbase.R")
cell<-DLPFC.refbase(beta=Combat2_Batch,
	                arraytype="450K")
cellTypes<-cell$CellFraction

##########sva
library(sva)
mod1 = model.matrix(~Sample_Group+sex+age_years+race+cellTypes$,data=myLoad$pd) #All cell proportion influence except the one with least cell proportion get corrected
mod0<-mod1[,-c(2,3)]
n.sv<-num.sv(Combat2_Batch, mod1, method = "be", vfilter = NULL, B = 20,seed = NULL)
sva_14 = sva(Combat2_Batch,mod1,mod0, n.sv =14 )$sv

colnames(sva_14)<-paste("sv", 1:14, sep = "")
pd<-cbind(myLoad$pd,sva_14)

###pvca
source("./pvca.R")
pvca(Combat2_Batch,pd)

#regression
design = model.matrix(~Sample_Group+sex+age_years+race+cellTypes$+sv$, data=pd)#desigh matrix for ref-based
Y<-Combat2_Batch
X<-as.matrix(design)
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
RegressSCZ_sexcombined= Y - t(X[,c(3:ncol(X))] %*% beta[c(3:nrow(beta)),]) # regress all except group
RegressSCZ_sexstratified= Y - t(X[,c(4:ncol(X))] %*% beta[c(4:nrow(beta)),]) # regress all except group sex



