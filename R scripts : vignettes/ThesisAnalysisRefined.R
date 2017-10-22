source("~/Desktop/Thesis/Scripts_Thesis/ThesisFunctions.R")
source("~/Desktop/Thesis/Scripts_Thesis/GOplot_function.R")
source("~/Desktop/Thesis/Scripts_Thesis/GeneOntology_functions.R")
source("~/Desktop/Thesis/Scripts_Thesis/GO-ALL.R")
load("~/Desktop/Thesis/Data_Thesis/ThesisWorlk.RData")
load("~/Desktop/Thesis/Data_Thesis/HUGO.RData")

####genes will be reported with gene names (HUGO symbol)


####DATA PREPARATION########################
gtex.annot<-read.delim("GTEx_Data_V6_Annotations_SampleAttributesDS.txt", sep='\t', header=TRUE)
discordant.lusc<-read.table("discordant_samples.txt",stringsAsFactors=FALSE)$V1
#####download lung data through GDC##########
#############################################


query.lung.htseq <- GDCquery(project = "TCGA-LUSC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification", 
                       workflow.type = "HTSeq - Counts")

samplesDown.lusc.htseq <- getResults(query.lung.htseq,cols=c("cases"))


dataSmTP.lusc.htseq <- TCGAquery_SampleTypes(barcode = samplesDown.lusc.htseq,
                                       typesample = "TP")

dataSmNT.lusc.htseq <- TCGAquery_SampleTypes(barcode = samplesDown.lusc.htseq,
                                       typesample = "NT")


######Tumor purity########
query.lung.htseq2 <- GDCquery(project = "TCGA-LUSC",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - Counts",
                             barcode = c(dataSmTP.lusc.htseq, dataSmNT.lusc.htseq))

GDCdownload(query=query.lung.htseq2)

dataPrep.tcga.full1<-GDCprepare(query=query.lung.htseq2,
                                  save=TRUE)
dataPrep.tcga.full<-TCGAanalyze_Preprocessing(object = dataPrep.tcga.full1, 
                                              cor.cut = 0.6)
rownames(dataPrep.tcga.full)<-rowData(dataPrep.tcga.full1)$external_gene_name

TP.pure<-TCGAtumor_purity(dataSmTP.lusc.htseq, 0, 0, 0, 0, 0.6)$pure_barcodes
TP.notpure<-setdiff(dataSmTP.lusc.htseq, TP.pure)

cancer.htseq.full<-get_IDs(dataPrep.tcga.full)$condition
cancer.htseq.full<-gsub("normal","N", cancer.htseq.full)
cancer.htseq.full<-gsub("cancer","C", cancer.htseq.full)

dataPrep.tcga.full.N.htseq<-TCGAanalyze_Normalization(tabDF = dataPrep.tcga.full,
                                                 geneInfo = geneInfo,
                                                 method = "gcContent")

dataPrep.tcga.full.NF.htseq<- TCGAanalyze_Filtering(tabDF = dataPrep.tcga.full.N.htseq,
                                               method = "quantile", 
                                               qnt.cut =  0.25)


mds.htseq.P<-myMDSplot(dataPrep.tcga.full, purityvec.tcga, cancer.htseq.full)
#########################
purityinfo.htseq<-TCGAtumor_purity(dataSmTP.lusc.htseq, 0, 0, 0, 0, 0.6)

###barcodes with 60% purity at least and no discordant samples####
dataSmTP.lusc.pure.htseq<-purityinfo.htseq$pure
dataSmTP.lusc.pure.htseq<-setdiff(dataSmTP.lusc.pure.htseq, discordant.lusc)


###HTseq Counts
queryDown.lung.htseq <- GDCquery(project = "TCGA-LUSC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification", 
                       workflow.type = "HTSeq - Counts",
                       barcode = c(dataSmTP.lusc.pure.htseq, dataSmNT.lusc.htseq))


GDCdownload(query = queryDown.lung.htseq)


dataPrep1.tcga.htseq <- GDCprepare(query = queryDown.lung.htseq, 
                              save = TRUE )



dataPrep.tcga.htseq <- TCGAanalyze_Preprocessing(object = dataPrep1.tcga.htseq, 
                                                cor.cut = 0.6
                                                #,datatype = "raw_count"
                                               )

rownames(dataPrep.tcga.htseq)<-rowData(dataPrep1.tcga.htseq)$external_gene_name


#####Not necessary#######
#rownames(dataPrep.tcga.htseq)<-gsub("\\|.*","",rownames(dataPrep.tcga.htseq))

q<-get_IDs(dataPrep.tcga.htseq)$condition
cancer.htseq<-gsub("normal","N", cancer.htseq)
cancer.htseq<-gsub("cancer","C", cancer.htseq)

plates.tcga<-c(as.character(get_IDs(dataPrep.tcga.NF)$plate))
plates.tcga<-make.names(plates.tcga)

batch.htseq<-rep("TCGA", ncol(dataPrep.tcga.htseq))
mds.htseq<-myMDSplot(dataPrep.tcga.htseq, batch.htseq, cancer.htseq)

dataPrep.tcga.N.htseq<-TCGAanalyze_Normalization(tabDF = dataPrep.tcga.htseq,
                                           geneInfo = geneInfo,
                                           method = "gcContent")

dataPrep.tcga.NF.htseq<- TCGAanalyze_Filtering(tabDF = dataPrep.tcga.N.htseq,
                                         method = "quantile", 
                                         qnt.cut =  0.25)


###########
dataPrep.tcga.NF<-dataPrep.tcga.NF.htseq

plates.tcga<-c(as.character(get_IDs(dataPrep.tcga.NF)$plate))
plates.tcga<-make.names(plates.tcga)

v.tcga<-voom(dataPrep.tcga.NF, design.matrix.tcga, plot=TRUE)

designDEA.limma.tcga <- model.matrix(~0+cancer.htseq+plates.tcga)
designDEA.limma.tcga.NP<-model.matrix(~0+cancer.htseq)

contr.matrix.limma.tcga<-makeContrasts(contrasts="cancer.htseqC-cancer.htseqN", levels=designDEA.limma.tcga)
contr.matrix.limma.tcga.NP<-makeContrasts(contrasts="cancer.htseqC-cancer.htseqN", levels=designDEA.limma.tcga.NP)

DEGs<-TCGAanalyze_DEA(dataPrep.tcga.NF[,dataSmNT.lusc.htseq],
                      dataPrep.tcga.NF[,dataSmTP.lusc.pure.htseq],"Normal", "Tumor", pipeline = "edgeR",
                      method = "glmLRT", fdr.cut = 0.01, logFC.cut = 1,
                      batch.factors = "Plate")



DEGs.limma<-DE_limma(contr.matrix.limma.tcga, v.tcga, designDEA.limma.tcga, 1, 0.01)
DEGs.limma.NP<-DE_limma(contr.matrix.limma.tcga.NP, v.tcga, designDEA.limma.tcga.NP, 1, 0.01)

tt<-DEGs
index.up <- which(tt$logFC >= 1 & tt$FDR < 0.01)
index.down <- which(tt$logFC <= -1 & tt$FDR < 0.01)
direction <- c()
direction[index.up] <- "up"
direction[index.down] <- "down"
direction[!(1:nrow(tt) %in% union(index.up,index.down))] <- "no DE"
tt <- cbind(tt,direction)

tt.limma<-DEGs.limma
index.up <- which(tt.limma$logFC >= 1 & tt.limma$adj.P.Val < 0.01)
index.down <- which(tt.limma$logFC <= -1 & tt.limma$adj.P.Val < 0.01)
direction <- c()
direction[index.up] <- "up"
direction[index.down] <- "down"
direction[!(1:nrow(tt.limma) %in% union(index.up,index.down))] <- "no DE"
tt.limma <- cbind(tt.limma,direction)

tt.limma.NP<-DEGs.limma.NP
index.up <- which(tt.limma.NP$logFC >= 1 & tt.limma.NP$adj.P.Val < 0.01)
index.down <- which(tt.limma.NP$logFC <= -1 & tt.limma.NP$adj.P.Val < 0.01)
direction <- c()
direction[index.up] <- "up"
direction[index.down] <- "down"
direction[!(1:nrow(tt.limma.NP) %in% union(index.up,index.down))] <- "no DE"
tt.limma.NP <- cbind(tt.limma.NP,direction)


edgeRdown.tcga<-rownames(tt[which(tt$direction == "down"),])
edgeRup.tcga<-rownames(tt[which(tt$direction == "up"),])
limmadown.tcga<-rownames(tt.limma[which(tt.limma$direction == "down"),])
limmaup.tcga<-rownames(tt.limma[which(tt.limma$direction == "up"),])
limmadown.tcga.NP<-rownames(tt.limma.NP[which(tt.limma.NP$direction == "down"),])
limmaup.tcga.NP<-rownames(tt.limma.NP[which(tt.limma.NP$direction == "up"),])

overlap(limmadown.tcga, limmadown.tcga.NP)

####Overlap ratio####
overlap(edgeRup.tcga, limmaup.tcga)

write.table(edgeRdown.tcga, file="edgeRdowntcga.txt")
write.table(edgeRup.tcga, file="edgeRuptcga.txt")
write.table(limmadown.tcga, file="limmadowntcga.txt")
write.table(limmaup.tcga, file="limmauptcga.txt")

#####DbGap GTEX lung Data, GDC TCGA#######
###########################################

####data from preprocessing function######

dataPrep.tcga<-dataPrep.tcga.htseq

load("GTex_Lung.Rdata")

###Remove genes with "." in the name
dataPrep.gtex<-dataGTEX_Tissue_only
dataPrep.gtex<-dataPrep.gtex[-which(grepl("\\.", rownames(dataPrep.gtex))),]


rownames(dataPrep.tcga)<-gsub("\\|.*","",rownames(dataPrep.tcga))
#rownames(dataPrep.gtex)<-gsub("\\..*","",rownames(dataPrep.gtex))
###common gene rows#####
dataPrep.gtex.filt<-dataPrep.gtex[intersect(rownames(dataPrep.gtex), rownames(dataPrep.tcga)),] ##GTEx dbgap healthy samples
dataPrep.tcga.filt<-dataPrep.tcga[intersect(rownames(dataPrep.gtex), rownames(dataPrep.tcga)),] ##TCGA GDC

dataPrep.tcga.tumor<-dataPrep.tcga.filt[,dataSmTP.lusc.pure.htseq] ##TCGA tumors
dataPrep.tcga.healthy<-dataPrep.tcga.filt[,dataSmNT.lusc.htseq] ##TCGA normal samples


dataPrep.all.healthy2<-merge(dataPrep.gtex.filt, dataPrep.tcga.healthy, by="row.names") ##All normal/healthy, tissues dbgap &GDC
rownames(dataPrep.all.healthy2)<-dataPrep.all.healthy2$Row.names
dataPrep.all.healthy2$Row.names<-NULL

dataPrep.all2<-merge(dataPrep.all.healthy2, dataPrep.tcga.tumor, by="row.names") ##All data: tumor and healthy
rownames(dataPrep.all2)<-dataPrep.all2$Row.names
dataPrep.all2$Row.names<-NULL
dim(dataPrep.all2)
####Condition and Project vector######

cancer.dbgap<-rep(c("N", "C"), c(ncol(dataPrep.all.healthy2), ncol(dataPrep.tcga.tumor)))
batch.dbgap<-rep(c("GTEX", "TCGA"), c(ncol(dataPrep.gtex.filt), ncol(dataPrep.tcga.filt)))
colorvec<-rep(c("GN", "TT", "TN"), c(ncol(dataPrep.gtex.filt), ncol(dataPrep.tcga.healthy), ncol(dataPrep.tcga.tumor) ))

##GTEx samples with batch info######
gtex.batch<-gtex.annot[which(gtex.annot$SAMPID %in% colnames(dataPrep.gtex)),]$SMNABTCH
gtex.batch<-make.names(as.character(gtex.batch))

####Keys for variables#############
#C: Combat corrected
#V: voom corrected
#N: Normalized
#F: filtered
#R: rnaseqdb
#L: Log transformed
#S: Standardised
#M: matrix
#E: eset object
#cqn: conditional quantile normalized 


####Normalization- GC content######
dataPrep.all2.N<-TCGAanalyze_Normalization(tabDF = dataPrep.all2,
                                                    geneInfo = geneInfo,
                                                    method = "gcContent")


dge.dbgap<-DGEList(dataPrep.all2.NF)
y.dbgap<-calcNormFactors(dge.dbgap)

###Filtering
dataPrep.all2.NF<- TCGAanalyze_Filtering(tabDF = dataPrep.all2.N,
                                         method = "quantile", 
                                         qnt.cut =  0.25)

dge.dbgap<-DGEList(dataPrep.all2.NF)
y.dbgap<-calcNormFactors(dge.dbgap)
plotMDS(y.dbgap, labels=NULL, pch=19 ,top=500, col=as.numeric(cancer.dbgap), gene.selection="common", prior.count=5)


###design matrices && batch info###
design.matrix<- model.matrix(~cancer.dbgap)
design.mod.combat<-model.matrix(~cancer.dbgap)
design.mod.combat.R<-model.matrix(~cancer.rnaseqdb)

plates<-c(gtex.batch, as.character(get_IDs(dataPrep.tcga.filt)$plate))
plates<-make.names(plates)
v<-voom(dataPrep.all2.NF, design.matrix, plot=TRUE)

###

ComBat(dataPrep.all2.NF, batch=batch.dbgap, mod=design.mod.combat, par.prior=TRUE, prior.plots = TRUE)
ComBat(lusc.all.NF, batch=plates.R, mod=design.mod.combat.R, par.prior=TRUE, prior.plots = TRUE)

###MDS unprocessed data
mds.dbgap<-myMDSplot(dataPrep.all2, batch.dbgap, cancer.dbgap)

###MDS before ComBat (NOT voom transformed)
mds.dbgap.N<-myMDSplot(dataPrep.all2.NF, batch.dbgap, cancer.dbgap)


####DEA######
plates<-make.names(plates)
length(unique(plates))

cancer.dbgap<-factor(x=cancer.dbgap, levels=unique(cancer.dbgap))
designDEA.limma <- model.matrix(~0+cancer.dbgap+plates)
designDEA.edgeR <- model.matrix(~cancer.dbgap+plates)

colnames(designDEA.edgeR)[1:length(levels(cancer.dbgap))]<-levels(cancer.dbgap)

aDGEList <- edgeR::DGEList(counts = as.matrix(dataPrep.all2.NF), group = cancer.dbgap)
aDGEList$samples$condition<-cancer.dbgap
aDGEList$samples$condition<-relevel(aDGEList$samples$condition, ref="N")

aDGEList <- edgeR::estimateGLMCommonDisp(aDGEList, designDEA.edgeR)
aDGEList <- edgeR::estimateGLMTagwiseDisp(aDGEList, designDEA.edgeR)
aGlmFit <- edgeR::glmFit(aDGEList, designDEA.edgeR, dispersion = aDGEList$tagwise.dispersion,
                         prior.count.total=0)
my.lrt <- edgeR::glmLRT(aGlmFit, coef = 2)

contr.matrix.limma<-makeContrasts(contrasts="cancer.dbgapC-cancer.dbgapN", levels=designDEA.limma)


tableDEAedgeR.test<-edgeR_DEA(Data=dataPrep.all2.NF,
                              design=designDEA.edgeR, 
                              groupvec=cancer.dbgap,
                              lfc=1, fdr=0.01, ref="N")

tableDEA.dbgap.limma<-DE_limma(contr.matrix.limma, v, designDEA.limma, 1, 0.01)
tableDEA.dbgap.edgeR<-DE_edgeR(my.lrt, aGlmFit, 1, 0.01)

############THIS!!##################

#####up and down regulated genes
edgeRdown<-rownames(tableDEA.dbgap.edgeR[which(tableDEA.dbgap.edgeR$direction == "down"),])
edgeRup<-rownames(tableDEA.dbgap.edgeR[which(tableDEA.dbgap.edgeR$direction == "up"),])
limmadown<-rownames(tableDEA.dbgap.limma[which(tableDEA.dbgap.limma$direction == "down"),])
limmaup<-rownames(tableDEA.dbgap.limma[which(tableDEA.dbgap.limma$direction == "up"),])

write.table(edgeRdown, file="edgeRdown.txt")
write.table(edgeRup, file="edgeRup.txt")
write.table(limmadown, file="limmadown.txt")
write.table(limmaup, file="limmaup.txt")


overlap(limmaup, limmaup.R)


martaedgeRup<-read.table("up_edgeR_LUSC_all_plate_TumorPurity.txt", header=FALSE)
martaedgeRdown<-read.table("down_edgeR_LUSC_all_plate_TumorPurity.txt", header=FALSE)
martalimmaup<-read.table("up_limma_LUSC_all_plate_TumorPurity.txt", header=FALSE)
martalimmadown<-read.table("down_limma_LUSC_all_plate_TumorPurity.txt", header=FALSE)

intersect(martalimmaup,limmaup)
length(edgeRdown)

head(martaedgeRup)
####RNASEQDB normalized data##########################
######################################################
lusc.gtex.rnaseqdb<-read.table("lung-rsem-count-gtex.txt", sep='\t', header=TRUE)
lusc.tcga.rnaseqdb.tumor<-read.table("lusc-rsem-count-tcga-t.txt", sep='\t', header=TRUE)
lusc.tcga.rnaseqdb.healthy<-read.table("lusc-rsem-count-tcga.txt", sep='\t', header=TRUE)

colnames(lusc.gtex.rnaseqdb)<-gsub("\\.", "-", colnames(lusc.gtex.rnaseqdb))
colnames(lusc.tcga.rnaseqdb.tumor)<-gsub("\\.", "-", colnames(lusc.tcga.rnaseqdb.tumor))
colnames(lusc.tcga.rnaseqdb.healthy)[3:length(colnames(lusc.tcga.rnaseqdb.healthy))]<-gsub("\\.", "-", colnames(lusc.tcga.rnaseqdb.healthy)[3:length(colnames(lusc.tcga.rnaseqdb.healthy))])


lusc.tcga.rnaseqdb.healthy$Hugo_Symbol<-NULL
lusc.tcga.rnaseqdb.healthy$Entrez_Gene_Id<-NULL

lusc.tcga.rnaseqdb.tumor$Hugo_Symbol<-NULL
lusc.tcga.rnaseqdb.tumor$Entrez_Gene_Id<-NULL
###barcodes with 60% purity at least and no discordant samples####
purityinfo.R<-TCGAtumor_purity(colnames(lusc.tcga.rnaseqdb.tumor), 0, 0, 0, 0, 0.6)
dataSmTP.lusc.pure.R<-purityinfo.R$pure
dataSmTP.lusc.pure.R<-setdiff(dataSmTP.lusc.pure.R, discordant.lusc)
lusc.tcga.rnaseqdb.tumor<-lusc.tcga.rnaseqdb.tumor[,dataSmTP.lusc.pure.R]

num_normal<-ncol(lusc.gtex.rnaseqdb)+ncol(lusc.tcga.rnaseqdb.healthy)-2
num_cancer<-ncol(lusc.tcga.rnaseqdb.tumor)
num_tcga<-ncol(lusc.tcga.rnaseqdb.tumor)+ncol(lusc.tcga.rnaseqdb.healthy)
num_gtex<-ncol(lusc.gtex.rnaseqdb)-2

cancer.rnaseqdb<-rep(c("N", "C"), c(num_normal,num_cancer))
batch.rnaseqdb<-rep(c("gtex", "tcga"), c(num_gtex,num_tcga)) 
colnames(lusc.gtex.rnaseqdb)

lusc.all.rnaseqdb<-cbind(lusc.gtex.rnaseqdb, lusc.tcga.rnaseqdb.healthy, lusc.tcga.rnaseqdb.tumor)
rownames(lusc.all.rnaseqdb)<-lusc.all.rnaseqdb$Hugo_Symbol
lusc.all.rnaseqdb$Hugo_Symbol<-NULL
lusc.all.rnaseqdb$Entrez_Gene_Id<-NULL



###MDS####

mds.rnaseqdb<-myMDSplot(lusc.all.rnaseqdb, batch.rnaseqdb, cancer.rnaseqdb)

##GTEx samples with batch info######
gtex.batch.rnaseqdb<-as.character(gtex.annot[which(gtex.annot$SAMPID %in% colnames(lusc.all.rnaseqdb)),]$SMNABTCH)

####Normalization- GC content######
lusc.all.rnaseqdb.N<-TCGAanalyze_Normalization(tabDF = lusc.all.rnaseqdb,
                                           geneInfo = geneInfo,
                                           method = "gcContent")

lusc.all.NF<- TCGAanalyze_Filtering(tabDF = lusc.all.rnaseqdb.N,
                                         method = "quantile", 
                                         qnt.cut =  0.25)


###RNASEQDB design matrices && batch info###
design.matrix.R<- model.matrix(~cancer.rnaseqdb)


plates.R<-c(gtex.batch.rnaseqdb, get_IDs(cbind(lusc.tcga.rnaseqdb.healthy, lusc.tcga.rnaseqdb.tumor))$plate)

v.R<-voom(lusc.all.NF, design.matrix.R, plot=TRUE)
####DEA RNASEQDB######
plates.R<-make.names(plates.R)
cancer.rnaseqdb<-factor(x=cancer.rnaseqdb, levels=unique(cancer.rnaseqdb))

designDEA.limma.R <- model.matrix(~0+cancer.rnaseqdb+plates.R)
designDEA.limma.R.NP <- model.matrix(~0+cancer.rnaseqdb)
designDEA.edgeR.R <- model.matrix(~cancer.rnaseqdb+plates.R)

colnames(designDEA.edgeR.R)[1:length(levels(cancer.rnaseqdb))]<-levels(cancer.rnaseqdb)

aDGEList.R <- edgeR::DGEList(counts = lusc.all.NF, group = cancer.rnaseqdb)
aDGEList.R$samples$condition<-cancer.rnaseqdb
aDGEList.R$samples$condition<-relevel(aDGEList.R$samples$condition, ref="N")

aDGEList.R <- edgeR::estimateGLMCommonDisp(aDGEList.R, designDEA.edgeR.R)
aDGEList.R <- edgeR::estimateGLMTagwiseDisp(aDGEList.R, designDEA.edgeR.R)
aGlmFit.R <- edgeR::glmFit(aDGEList.R, designDEA.edgeR.R, dispersion = aDGEList.R$tagwise.dispersion,
                         prior.count.total=0)
my.lrt.R <- edgeR::glmLRT(aGlmFit.R, coef = 2)

contr.matrix.limma.R<-makeContrasts(contrasts="cancer.rnaseqdbC-cancer.rnaseqdbN", levels=designDEA.limma.R)
contr.matrix.limma.R.NP<-makeContrasts(contrasts="cancer.rnaseqdbC-cancer.rnaseqdbN", levels=designDEA.limma.R.NP)


tableDEA.limma.R<-DE_limma(contr.matrix.limma.R, v.R, designDEA.limma.R, 1, 0.01)
tableDEA.limma.R.NP<-DE_limma(contr.matrix.limma.R.NP, v.R, designDEA.limma.R.NP, 1, 0.01)

tableDEA.edgeR.R<-DE_edgeR(my.lrt.R, aDGEList.R, 1, 0.01)

tableDEA.edgeR.R<-edgeR_DEA(Data=lusc.all.NF,
          design=designDEA.edgeR.R, 
          groupvec=cancer.rnaseqdb,
          lfc=1, fdr=0.01, ref="N")

#####up and down regulated genes
edgeRdown.R<-rownames(tableDEA.edgeR.R[which(tableDEA.edgeR.R$direction == "down"),])
edgeRup.R<-rownames(tableDEA.edgeR.R[which(tableDEA.edgeR.R$direction == "up"),])

limmadown.R<-rownames(tableDEA.limma.R[which(tableDEA.limma.R$direction == "down"),])
limmaup.R<-rownames(tableDEA.limma.R[which(tableDEA.limma.R$direction == "up"),])

limmadown.R.NP<-rownames(tableDEA.limma.R.NP[which(tableDEA.limma.R.NP$direction == "down"),])
limmaup.R.NP<-rownames(tableDEA.limma.R.NP[which(tableDEA.limma.R.NP$direction == "up"),])

overlap(limmadown.R, limmadown.R.NP)

write.table(edgeRdown.R, file="edgeRdownR.txt")
write.table(edgeRup.R, file="edgeRupR.txt")
write.table(limmadown.R, file="limmadownR.txt")
write.table(limmaup.R, file="limmaupR.txt")

######Vendiagrams########
VennPlot(list(limmaup, edgeRup))
VennPlot(list(limmaup.tcga, edgeRup.tcga))

VennPlot(list(limmaup, limmaup.tcga, limmaup.R))
VennPlot(list(limmadown, limmadown.tcga, limmadown.R))
VennPlot(list(edgeRdown, edgeRdown.tcga, edgeRdown.R))
#########################

####PCA##########################################

###tumor purity plotting####
purityvec.tcga<-c()
purityvec.R<-c()
purityvec.dbgap<-c()
cancer.dbgap.full<-c()
cancer.lusc.full<-c()

####TCGA only#####
for(i in colnames(dataPrep.tcga.full)){
  if(i %in% TP.pure){
    purityvec.tcga<-c(purityvec.tcga, "P")
    
  } 
  else if(i %in% TP.notpure){
    purityvec.tcga<-c(purityvec.tcga, "M")
  } 
  else{
    purityvec.tcga<-c(purityvec.tcga, "H")
  }
}

###GTEx+TCGA#####

dataPrep.gtex.filt2<-dataPrep.gtex[intersect(rownames(dataPrep.gtex), rownames(dataPrep.tcga.full)),] ##GTEx dbgap healthy samples
dataPrep.tcga.filt2<-dataPrep.tcga.full[intersect(rownames(dataPrep.gtex), rownames(dataPrep.tcga.full)),] ##TCGA GDC

dataPrep.tcga.tumor2<-dataPrep.tcga.filt2[,dataSmTP.lusc.htseq] ##TCGA tumors
dataPrep.tcga.healthy2<-dataPrep.tcga.filt2[,dataSmNT.lusc.htseq] ##TCGA normal samples

dataPrep.all.healthy2.P<-merge(dataPrep.gtex.filt2, dataPrep.tcga.healthy2, by="row.names") ##All normal/healthy, tissues dbgap &GDC
rownames(dataPrep.all.healthy2.P)<-dataPrep.all.healthy2.P$Row.names
dataPrep.all.healthy2.P$Row.names<-NULL

dataPrep.all2.full<-merge(dataPrep.all.healthy2.P, dataPrep.tcga.tumor2, by="row.names") ##All data: tumor and healthy
rownames(dataPrep.all2.full)<-dataPrep.all2$Row.names
dataPrep.all2.full$Row.names<-NULL


cancer.dbgap.full<-c()

for(i in colnames(dataPrep.all2.full)){
  if(i %in% TP.pure){
    purityvec.dbgap<-c(purityvec.dbgap,"P")
    cancer.dbgap.full<-c(cancer.dbgap.full,"C")
  } 
  else if(i %in% TP.notpure){
    purityvec.dbgap<-c(purityvec.dbgap, "M")
    cancer.dbgap.full<-c(cancer.dbgap.full,"C")
  } 
  else if(grepl("GTEX", i)){
    purityvec.dbgap<-c(purityvec.dbgap, "G")
    cancer.dbgap.full<-c(cancer.dbgap.full,"N")
  }
  else{
    purityvec.dbgap<-c(purityvec.dbgap, "NT")
    cancer.dbgap.full<-c(cancer.dbgap.full, "N")
  }
}

####RNASEQDB######

lusc.gtex.rnaseqdb<-read.table("lung-rsem-count-gtex.txt", sep='\t', header=TRUE)
lusc.tcga.rnaseqdb.tumor<-read.table("lusc-rsem-count-tcga-t.txt", sep='\t', header=TRUE)
lusc.tcga.rnaseqdb.healthy<-read.table("lusc-rsem-count-tcga.txt", sep='\t', header=TRUE)
nrow(lusc.all.rnaseqdb)

colnames(lusc.gtex.rnaseqdb)<-gsub("\\.", "-", colnames(lusc.gtex.rnaseqdb))
colnames(lusc.tcga.rnaseqdb.tumor)<-gsub("\\.", "-", colnames(lusc.tcga.rnaseqdb.tumor))
colnames(lusc.tcga.rnaseqdb.healthy)[3:length(colnames(lusc.tcga.rnaseqdb.healthy))]<-gsub("\\.", "-", colnames(lusc.tcga.rnaseqdb.healthy)[3:length(colnames(lusc.tcga.rnaseqdb.healthy))])


lusc.tcga.rnaseqdb.healthy$Hugo_Symbol<-NULL
lusc.tcga.rnaseqdb.healthy$Entrez_Gene_Id<-NULL

lusc.tcga.rnaseqdb.tumor$Hugo_Symbol<-NULL
lusc.tcga.rnaseqdb.tumor$Entrez_Gene_Id<-NULL


lusc.all.full<-cbind(lusc.gtex.rnaseqdb, lusc.tcga.rnaseqdb.healthy, lusc.tcga.rnaseqdb.tumor)
rownames(lusc.all.full)<-lusc.all.full$Hugo_Symbol
lusc.all.full$Hugo_Symbol<-NULL
lusc.all.full$Entrez_Gene_Id<-NULL

purityvec.R<-c()
for(i in colnames(lusc.all.full)){
  if(i %in% TP.pure){
    purityvec.R<-c(purityvec.R,"P")
    cancer.lusc.full<-c(cancer.lusc.full,"C")
  } 
  else if(i %in% TP.notpure){
    purityvec.R<-c(purityvec.R, "M")
    cancer.lusc.full<-c(cancer.lusc.full,"C")
  } 
  else if(grepl("GTEX", i)){
    purityvec.R<-c(purityvec.R, "G")
    cancer.lusc.full<-c(cancer.lusc.full,"G")
  }
  else{
    purityvec.R<-c(purityvec.R, "NT")
    cancer.lusc.full<-c(cancer.lusc.full, "NT")
  }
}
purityvec.R[705]


####PCA s################
Aov.corr(lusc.all.NF, plates.R)
title("Analysis of variance plot: RNASEQDB")
Aov.corr(dataPrep.all2.NF, plates)
title("Analysis of variance plot")
###########

Fstats(lusc.all.rnaseqdb, plates.R)
Fstats(dataPrep.all2, plates)

PC.1D(dataPrep.all2.NF, plates)
title("1D PCA plots for each plate/batch for all samples")
PC.1D(lusc.all.NF, plates.R)
title("1D PCA for each plate/batch for all samples")

pca.tcga.full<-myPCA(t(dataPrep.tcga.full), purityvec.tcga, title="1D PCA with tumor purity grouping")
pca.tcga.full.lib<-myPCA(t(libsize.corr(dataPrep.tcga.full)), purityvec.tcga, title="1D PCA with tumor purity grouping")
pca.tcga.full$melted

pca.dbgap.full<-myPCA(t(dataPrep.all2.full), purityvec.dbgap, title="PCA plot(s) with tumor purity grouping")
pca.lusc.full<-myPCA(t(lusc.all.full), purityvec.R, "PCA plot(s) with tumor purity grouping (RNASEQDB)")

pca.dbgap.full$melted
pca.lusc.full$melted

ggsave(filename="1DPCAfull.jpg", path="~/Desktop/Thesis/ku-template-latex-thesis/images", dpi=500, plot=pca.dbgap.full$melted)
ggsave(filename="1DPCAfullR.jpg", path="~/Desktop/Thesis/ku-template-latex-thesis/images", dpi=500, plot=pca.lusc.full$melted)


pca.lusc.full$PC2
pca.dbgap.full$PC2

pca.dbgap.NF<-myPCA(t(dataPrep.all2.NF), cancer.dbgap, title="1D PCA plots of normal vs cancer samples")
pca.dbgap.NF.lib<-myPCA(t(libsize.corr(dataPrep.all2.NF)), cancer.dbgap, title="PCA plot of normal vs cancer samples")

pca.dbgap.NF$melted

pca.dbgap.full$melte
pca.dbgap.full$PC2<-pca.dbgap.full$PC2+ggtitle("PCA plot of normal vs cancer samples")
pca.lusc.full$PC2


dev.off()
par(mfrow = c(3, 3)) 

# 3 rows and 2 columns
pcadata<-dataPrep.all2.full
pcatest<-prcomp(t(pcadata), center=TRUE)
my.data<-as.data.frame(t(pcadata))
my.data$grouping<-cancer.dbgap.full
PCi<-data.frame(pcatest$x[,1:9],group=cancer.dbgap.full)
PCi.melt<-melt(PCi, id.vars=c('group'))


VennPlot(list(limmadown.R, 
              limmadown.noNAT.R, 
              limmadown.noNAT.R.NP, 
              limmadown.tcga,
              limmadown.noNAT), 
         c("L.R", "L.noNAT.R", "L.noNAt.R.NP", "L.tcga", "L.noNAT"))

VennM1up<-VennPlot(list(limmaup.tcga,
              limmaup, limmaup.R,
              limmaup.noNAT,
              limmaup.noNAT.R), 
         c("Nt+", "Ntg+", "Rtg+", "Ntg-", "Rtg-" ))
grid.arrange(gTree(children=VennM1up), top="Venn diagram for up-regulated gene lists")

VennM1down<-VennPlot(list(limmadown.tcga,
              limmadown, limmadown.R,
              limmadown.noNAT,
              limmadown.noNAT.R), 
         c("Nt+", "Ntg+", "Rtg+", "Ntg-", "Rtg-" ))
grid.arrange(gTree(children=VennM1down), top="Venn diagram for down-regulated gene lists")
#######
VennM2up<-VennPlot(list(limmaup.tcga.NP,
                        limmaup.R.NP, limmaup.noNAT.NP,
                        limmaup.noNAT.R.NP), 
                   c("Nt+", "Rtg+", "Ntg-", "Rtg-" ))
grid.arrange(gTree(children=VennM2up), top="Venn diagram for up-regulated gene lists")

VennM2down<-VennPlot(list(limmadown.tcga.NP,
                          limmadown.R.NP, limmadown.noNAT.NP,
                          limmadown.noNAT.R.NP), 
                     c("Nt+", "Rtg+", "Ntg-", "Rtg-" ))
grid.arrange(gTree(children=VennM2down), top="Venn diagram for down-regulated gene lists")

