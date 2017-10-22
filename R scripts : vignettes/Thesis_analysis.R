####NAT vs GTEx######

NATvsGTEx<-cbind(lusc.gtex.rnaseqdb, lusc.tcga.rnaseqdb.healthy)
rownames(NATvsGTEx)<-NATvsGTEx$Hugo_Symbol
NATvsGTEx$Hugo_Symbol<-NULL
NATvsGTEx$Entrez_Gene_Id<-NULL

NAT<-c()
  for(i in colnames(NATvsGTEx)){
    if(grepl("GTEX", i)){
      NAT<-c(NAT, "G")
    }
    else
      NAT<-c(NAT, "N")
  }

NATvsGTEx.N<-TCGAanalyze_Normalization(tabDF = NATvsGTEx,
                                            geneInfo = geneInfo,
                                            method = "gcContent")

NATvsGTEx.NF<- TCGAanalyze_Filtering(tabDF = NATvsGTEx.N,
                                          method = "quantile", 
                                          qnt.cut =  0.25)

designDEA.limma.NAT<-model.matrix(~0+NAT)
design.matrix.NAT<- model.matrix(~NAT)
v.NAT<-voom(NATvsGTEx.NF, design.matrix.NAT, plot=TRUE)

contr.NAT<-makeContrasts(contrasts="NATN-NATG", levels=designDEA.limma.NAT)

tableDEA.limma.NAT<-DE_limma(contr.NAT, v.NAT, designDEA.limma.NAT, 1, 0.01)
limmadown.NAT<-rownames(tableDEA.limma.NAT[which(tableDEA.limma.NAT$direction == "down"),])
limmaup.NAT<-rownames(tableDEA.limma.NAT[which(tableDEA.limma.NAT$direction == "up"),])
######No NAT########

##U: No unique batches
##NP: No correction for plates
dataPrep.noNAT<-merge(dataPrep.gtex.filt, dataPrep.tcga.tumor, by="row.names") ##All data: tumor and healthy
rownames(dataPrep.noNAT)<-dataPrep.noNAT$Row.names
dataPrep.noNAT$Row.names<-NULL


purityinfo.R<-TCGAtumor_purity(colnames(lusc.tcga.rnaseqdb.tumor), 0, 0, 0, 0, 0.6)
dataSmTP.lusc.pure.R<-purityinfo.R$pure
dataSmTP.lusc.pure.R<-setdiff(dataSmTP.lusc.pure.R, discordant.lusc)
lusc.tcga.rnaseqdb.tumor<-lusc.tcga.rnaseqdb.tumor[,dataSmTP.lusc.pure.R]

lusc.noNAT<-cbind(lusc.gtex.rnaseqdb, lusc.tcga.rnaseqdb.tumor)

rownames(lusc.noNAT)<-lusc.noNAT$Hugo_Symbol
lusc.noNAT$Hugo_Symbol<-NULL
lusc.noNAT$Entrez_Gene_Id<-NULL
ncol(lusc.noNAT)

####Normalization- GC content######
dataPrep.noNAT.N<-TCGAanalyze_Normalization(tabDF = dataPrep.noNAT,
                                           geneInfo = geneInfo,
                                           method = "gcContent")

dataPrep.noNAT.NF<- TCGAanalyze_Filtering(tabDF = dataPrep.noNAT.N,
                                         method = "quantile", 
                                         qnt.cut =  0.25)

lusc.noNAT.N<-TCGAanalyze_Normalization(tabDF = lusc.noNAT,
                                               geneInfo = geneInfo,
                                               method = "gcContent")

lusc.noNAT.NF<- TCGAanalyze_Filtering(tabDF = lusc.noNAT.N,
                                    method = "quantile", 
                                    qnt.cut =  0.25)

######DEA no NAT#######
cancer.noNAT<-rep(c("N", "C"), c(ncol(dataPrep.gtex.filt), ncol(dataPrep.tcga.tumor)))
batch.noNAT<-rep(c("GTEX", "TCGA"), c(ncol(dataPrep.gtex.filt), ncol(dataPrep.tcga.tumor)))
plates.noNAT<-c(gtex.batch, as.character(get_IDs(dataPrep.tcga.tumor)$plate))
plates.noNAT<-make.names(plates.noNAT)

cancer.noNAT.R<-rep(c("N", "C"), c(ncol(lusc.gtex.rnaseqdb)-2, ncol(lusc.tcga.rnaseqdb.tumor)))
batch.noNAT.R<-rep(c("GTEX", "TCGA"), c(ncol(lusc.gtex.rnaseqdb)-2, ncol(lusc.tcga.rnaseqdb.tumor)))
plates.R.noNAT<-c(gtex.batch.rnaseqdb, get_IDs(lusc.tcga.rnaseqdb.tumor)$plate)
plates.R.noNAT<-make.names(plates.R.noNAT)



designDEA.limma.noNAT <- model.matrix(~0+cancer.noNAT+plates.noNAT)
designDEA.limma.noNAT.R <- model.matrix(~0+cancer.noNAT.R+plates.R.noNAT)
designDEA.limma.noNAT.NP <- model.matrix(~0+cancer.noNAT)
designDEA.limma.noNAT.NP.R <- model.matrix(~0+cancer.noNAT.R)

design.matrix.noNAT<- model.matrix(~cancer.noNAT)
design.matrix.noNAT.R<-model.matrix(~cancer.noNAT.R)

ComBat(lusc.noNAT.NF, batch=plates.R.noNAT, mod=model.matrix(~1, data=cancer.noNAT.R), par.prior=TRUE, prior.plots = TRUE)


#dataPrep.noNAT.NF.comb<-ComBat(lusc.noNAT.NF,  batch=plates.R.noNAT, par.prior=TRUE, prior.plots = TRUE)

v.noNAT<-voom(dataPrep.noNAT.NF, design.matrix.noNAT, plot=TRUE)
v.noNAT.R<-voom(lusc.noNAT.NF, design.matrix.noNAT.R, plot=TRUE)
v.noNAT.R2<-voom(lusc.noNAT.NF, designDEA.limma.noNAT.NP.R, plot=TRUE)

contr.matrix.limma.noNAT<-makeContrasts(contrasts="cancer.noNATC-cancer.noNATN", levels=designDEA.limma.noNAT)
contr.matrix.limma.noNAT.R<-makeContrasts(contrasts="cancer.noNAT.RC-cancer.noNAT.RN", levels=designDEA.limma.noNAT.R)
contr.matrix.limma.noNAT.NP<-makeContrasts(contrasts="cancer.noNATC-cancer.noNATN", levels=designDEA.limma.noNAT.NP)
contr.matrix.limma.noNAT.NP.R<-makeContrasts(contrasts="cancer.noNAT.RC-cancer.noNAT.RN", levels=designDEA.limma.noNAT.NP.R)

tableDEA.limma.noNAT<-DE_limma(contr.matrix.limma.noNAT, v.noNAT, designDEA.limma.noNAT, 1, 0.01)
tableDEA.limma.noNAT.R<-DE_limma(contr.matrix.limma.noNAT.R, v.noNAT.R, designDEA.limma.noNAT.R, 1, 0.01)
tableDEA.limma.noNAT.NP<-DE_limma(contr.matrix.limma.noNAT.NP, v.noNAT, designDEA.limma.noNAT.NP, 1, 0.01)
tableDEA.limma.noNAT.NP.R<-DE_limma(contr.matrix.limma.noNAT.NP.R, v.noNAT.R, designDEA.limma.noNAT.NP.R, 1, 0.01)

tableDEA.limma.noNAT.R2<-DE_limma(contr.matrix.limma.noNAT.NP.R, v.noNAT.R2, designDEA.limma.noNAT.NP.R, 1, 0.01)

overlap(limmadown.noNAT, limmaup)
overlap(limmaup.noNAT.R, limmaup.R)
overlap(limmadown.noNAT.R.NP, limmadown)

limmadown.noNAT<-rownames(tableDEA.limma.noNAT[which(tableDEA.limma.noNAT$direction == "down"),])
limmaup.noNAT<-rownames(tableDEA.limma.noNAT[which(tableDEA.limma.noNAT$direction == "up"),])

limmadown.noNAT.NP<-rownames(tableDEA.limma.noNAT.NP[which(tableDEA.limma.noNAT.NP$direction == "down"),])
limmaup.noNAT.NP<-rownames(tableDEA.limma.noNAT.NP[which(tableDEA.limma.noNAT.NP$direction == "up"),])

limmadown.noNAT.R<-rownames(tableDEA.limma.noNAT.R[which(tableDEA.limma.noNAT.R$direction == "down"),])
limmaup.noNAT.R<-rownames(tableDEA.limma.noNAT.R[which(tableDEA.limma.noNAT.R$direction == "up"),])

limmadown.noNAT.R.NP<-rownames(tableDEA.limma.noNAT.NP.R[which(tableDEA.limma.noNAT.NP.R$direction == "down"),])
limmaup.noNAT.R.NP<-rownames(tableDEA.limma.noNAT.NP.R[which(tableDEA.limma.noNAT.NP.R$direction == "up"),])

########################################################################
########################################################################


#########correcting for library size#####

y<-calcNormFactors(DGEList(counts = as.matrix(dataPrep.all2.NF), group = cancer.dbgap))
eff.y<-y$samples$norm.factors*y$samples$lib.size
y.norm<-y$counts/eff.y
head(y.norm)
v.norm.dbgap<-voom(y, design.matrix, plot=TRUE)
mds.dbgap2<-myMDSplot(y,batch.dbgap,cancer.dbgap)

colmds.lusc<-c()
for(i in cancer.)

plotMDS(DGEList(counts = as.matrix(lusc.all.full)), labels=NULL, pch=19 ,top=1800, col=as.numeric(as.factor(purityvec.R)), gene.selection="common", prior.count=5)

legend("bottomleft", 
       legend = c("G", "NT", "P", "M"), 
       col = c("red", "cyan", "purple", "green"),
       pch = c(19,19))
title("MDS plot of combined GTEx and TCGA datasets: RNASEQDB")



plotMDS(DGEList(counts = as.matrix(dataPrep.all2.full)), labels=NULL, pch=19 ,top=1800, col=as.numeric(as.factor(purityvec.dbgap)), gene.selection="common", prior.count=5)

legend("bottomleft", 
       legend = c("G", "NT", "P", "M"), 
       col = c("red", "cyan", "purple", "green"),
       pch = c(19,19))
title("MDS plot of combined GTEx and TCGA datasets")

tableDEA.dbgap.limma.norm<-DE_limma(contr.matrix.limma, v.norm.dbgap, designDEA.limma, 1, 0.01)
limmadown.norm<-rownames(tableDEA.dbgap.limma.norm[which(tableDEA.dbgap.limma.norm$direction == "down"),])
limmaup.norm<-rownames(tableDEA.dbgap.limma.norm[which(tableDEA.dbgap.limma.norm$direction == "up"),])

#######EA barplots#########
EA<-TCGAanalyze_EAcomplete(TFname="DEA genes", RegulonList =
                             DEGs.limma.C$ID)

TCGAvisualize_EAbarplot(tf = rownames(EA$ResBP),
                        GOBPTab = EA$ResBP, 
                        GOCCTab = EA$ResCC,
                        GOMFTab = EA$ResMF, 
                        PathTab = EA$ResPat,  
                        nRGTab = DEGs.limma.C$ID,
                        nBar = 20,
                        filename = "EA_Barplot.pdf")

#####gPCA batch #####
metadata.pvca <- data.frame(labelDescription=c("Tumor status",
                                          "plate"),
                       row.names=c("cancer", "batch"))

Adf<-AnnotatedDataFrame(data.frame(row.names=colnames(lusc.noNAT.NF), cancer=as.factor(cancer.noNAT.R),
                                   batch=as.factor(plates.R.noNAT)),
                        varMetadata = metadata.pvca)

eset.noNAT.R<-ExpressionSet(as.matrix(lusc.noNAT), phenoData = Adf )

pvca.R.noNAT<-pvcaBatchAssess(eset.noNAT.R, batch.factors = c("batch"), threshold = 0.6)


bp <- barplot(pvcaObj$dat, xlab = "Effects",
              ylab = "Weighted average proportion variance",
              ylim= c(0,1.1),col = c("blue"),
              las=2,main="PVCA estimation bar chart")

gpca.tcga<-gPCA.batchdetect(t(dataPrep.tcga), as.numeric(as.factor(plates.tcga)), filt = 1500 )
gpca.all2.b<-gPCA.batchdetect(t(dataPrep.all2), as.numeric(as.factor(batch.dbgap)), filt = 1500 )


gpca.all2<-gPCA.batchdetect(t(dataPrep.all2), as.numeric(as.factor(plates)), filt = 1500 )
gpca.R<-gPCA.batchdetect(t(lusc.all.rnaseqdb), as.numeric(as.factor(plates.R)), filt = 1500 )

gpca.R.noNAT.NF<-gPCA.batchdetect(t(lusc.noNAT.NF), as.numeric(as.factor(plates.R.noNAT)), filt = 1500 )
gpca.all2.noNAT<-gPCA.batchdetect(t(dataPrep.noNAT), as.numeric(as.factor(plates.noNAT)), filt= 1500)

PCplot(gpca.R.noNAT.NF,ug="guided",type="1v2")
PCplot(gpca.R.noNAT.NF,ug="guided",type="comp",npcs=3)

####Hist for p values####
qplot(c(tableDEA.dbgap.limma$adj.P.Val, 
        tableDEA.limma.R$adj.P.Val,
        tableDEA.limma.noNAT$adj.P.Val,
        tableDEA.limma.noNAT.R$adj.P.Val,
        tableDEA.limma.noNAT.NP$adj.P.Val,
        tableDEA.limma.R.NP$adj.P.Val,
        tableDEA.limma.noNAT.NP.R$adj.P.Val),
      xlab = "FDR ajusted p-values",
      geom="histogram",
      fill=I("grey"), 
      col=I("red"),
      bins=50,
      main = "FDR adjusted p-values distribution")

qplot(c(tableDEA.limma.NAT$adj.P.Val),
      xlab = "FDR ajusted p-values",
      geom="histogram",
      fill=I("grey"), 
      col=I("red"),
      bins=50,
      main = "FDR adjusted p-values distribution") 

gpca.R.noNAT.NF$p.val
