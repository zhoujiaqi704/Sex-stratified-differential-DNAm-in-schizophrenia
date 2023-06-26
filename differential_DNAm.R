library(ChAMP)
library(limma)
library(RRHO)

load("./pd.rda")
load("./RegressSCZ_sexstratified.rda")
load("./RegressSCZ_sexcombined.rda")

#sex-stratified
pd_f<-pd[grep("F",pd$sex),]
pd_m<-pd[grep("M",pd$sex),]

RegressSCZ_Female<-RegressSCZ_sexstratified[,match(pd_f$Sample_Name,colnames(RegressSCZ_sexstratified))]
RegressSCZ_Male<-RegressSCZ_sexstratified[,match(pd_m$Sample_Name,colnames(RegressSCZ_sexstratified))]

DMPscz_female <- champ.DMP(beta = RegressSCZ_Female,pheno=pd_f$Sample_Group,adjPVal=1)#
DMPscz_male <- champ.DMP(beta = RegressSCZ_Male,pheno=pd_m$Sample_Group,adjPVal=1)#
DMPscz_sexcombined <- champ.DMP(beta = RegressSCZ_sexcombined,pheno=pd_sva14_178HQ$Sample_Group,adjPVal=1)#

#interaction
mod <- model.matrix(~Sample_Group+sex+Sample_Group*sex+age_years+race+cellTypes$+sv , data = pd)
fit <- lmFit(Combat2_Batch, mod)
fitEb <- eBayes(fit)
options(digits=4)
subset_ALL_interactionEffect <- topTable(fitEb, num=Inf, coef=24)

#RRHO
rrho.female=DMPscz_female[,c(1,4)]
rownames(rrho.female)=rownames(DMPscz_female)
rrho.female$DMPvalue=-log10(rrho.female$P.Value)
rrho.female.down=rrho.female[which(rrho.female$logFC<0),]
rrho.female.up=rrho.female[which(rrho.female$logFC>0),]
rrho.female.down$DMPvalue=-(rrho.female.down$DMPvalue)
rrho.female=rbind(rrho.female.down,rrho.female.up)

rrho.male=DMPscz_male[,c(1,4)]
rrho.male=cbind(rownames(DMPscz_male),rrho.male)
rrho.male$DMPvalue=-log10(rrho.male$P.Value)
rrho.male.down=rrho.male[which(rrho.male$logFC<0),]
rrho.male.up=rrho.male[which(rrho.male$logFC>0),]
rrho.male.down$DMPvalue=-(rrho.male.down$DMPvalue)
rrho.male=rbind(rrho.male.down,rrho.male.up)

rrho.f=rrho.female[,c(4,3)]
rrho.m=rrho.male[,c(1,4)]

RRHO.SCZ<-RRHO(rrho.f,rrho.m,labels=c('SCZ female','SCZ male'),alternative="enrichment",plots = TRUE,outputdir = getwd(),BY = TRUE,log10.ind=FALSE)

#function for principle component slope
pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  return(list(slope,intercept))
}
