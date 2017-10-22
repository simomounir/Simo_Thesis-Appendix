###Pathway Enrichment Analysis#####


####M1
setwd("~/Desktop/Thesis/Data_Thesis/M1pathway")

path.tcgaUP<-pathwayRes(dataPrep.tcga.NF, 
                        limmaup.tcga,
                        "up-tcgaP", TRUE, 0.05)

path.tcgaDOWN<-pathwayRes(dataPrep.tcga.NF, 
                          limmadown.tcga,
                          "down-tcgaP", TRUE, 0.05)


path.dbgapUP<-pathwayRes(dataPrep.all2.NF, 
                          limmaup,
                          "up-dbgapP", TRUE, 0.05)

path.dbgapDOWN<-pathwayRes(dataPrep.all2.NF, 
                         limmadown,
                         "down-dbgapP", TRUE, 0.05)


path.RDOWN<-pathwayRes(lusc.all.NF, 
                          limmadown.R,
                           "down-RP", TRUE, 0.05)

path.RUP<-pathwayRes(lusc.all.NF, 
                         limmaup.R,
                         "up-RP", TRUE, 0.05)


path.R.noNAT.DOWN<-pathwayRes(lusc.noNAT.NF, 
                       limmadown.noNAT.R,
                       "down-noNAT-RP", TRUE, 0.05)


path.R.noNAT.UP<-pathwayRes(lusc.noNAT.NF, 
                     limmaup.noNAT.R,
                     "up-noNAT-RP", TRUE, 0.05)


path.noNAT.DOWN<-pathwayRes(dataPrep.noNAT.NF, 
                       limmadown.noNAT,
                       "down-noNATP", TRUE, 0.05)


path.noNAT.UP<-pathwayRes(dataPrep.noNAT.NF, 
                     limmaup.noNAT,
                     "up-noNATP", TRUE, 0.05)

###NATvsGTEx####
path.NATUP<-pathwayRes(NATvsGTEx.NF, 
                        limmaup.NAT,
                        "up-NAT", TRUE, 0.05)

path.NATDOWN<-pathwayRes(NATvsGTEx.NF, 
                          limmadown.NAT,
                          "down-NAT", TRUE, 0.05)


write.csv(path.tcgaUP, "up-tcga.csv")
write.csv(path.tcgaDOWN, "down-tcga.csv")
write.csv(path.dbgapUP, "up-dbgap.csv")
write.csv(path.dbgapDOWN, "down-dbgap.csv")
write.csv(path.RUP, "up-R.csv")
write.csv(path.RDOWN, "down-R.csv")
write.csv(path.R.noNAT.DOWN, "down-noNAT-R.csv")
write.csv(path.R.noNAT.UP, "up-noNAT-R.csv")
write.csv(path.noNAT.DOWN, "down-noNAT.csv")
write.csv(path.noNAT.UP, "up-noNAT.csv")

write.csv(path.NATDOWN, "down-NAT.csv")
write.csv(path.NATUP, "up-NAT.csv")

cons.data<-merge(lusc.all.NF, dataPrep.all2.NF, by="row.names") ##All normal/healthy, tissues dbgap &GDC
rownames(cons.data)<-cons.data$Row.names
cons.data$Row.names<-NULL


path.down.consensus<-pathwayRes(cons.data, 
                                unique(c(limmadown.tcga,limmadown,
                                                       limmadown.R, limmadown.noNAT,
                                                       limmadown.noNAT.R, limmadown.tcga.NP,
                                         limmadown.R.NP, limmadown.noNAT.R.NP,
                                         limmadown.noNAT.NP)),
                                "consensusdownALL2", TRUE, 0.05)

path.down.consensusNONAT<-pathwayRes(cons.data, 
                                unique(c(limmadown.noNAT,
                                         limmadown.noNAT.R,
                                         limmadown.noNAT.R.NP,
                                         limmadown.noNAT.NP)),
                                "consensusdownNONAT", TRUE, 0.05)

path.down.consensusNONAT.inter<-pathwayRes(cons.data, 
                                           Reduce(intersect, list(limmadown.noNAT,
                                                                  limmadown.noNAT.R,
                                                                  limmadown.noNAT.R.NP,
                                                                  limmadown.noNAT.NP)),
                                     "consensusdownNONATinter", TRUE, 0.05)


path.up.consensus<-pathwayRes(cons.data, 
                                unique(c(limmaup.tcga,limmaup,
                                         limmaup.R, limmaup.noNAT,
                                         limmaup.noNAT.R, limmaup.tcga.NP,
                                         limmaup.R.NP, limmaup.noNAT.R.NP,
                                         limmaup.noNAT.NP)),
                                "consensusupALL2", TRUE, 0.05)

View(path.up.consensus)

path.down.consensus<-pathwayRes(cons.data, 
                                unique(c(limmadown.tcga.NP,
                                         limmadown.R.NP, limmadown.noNAT.R.NP,
                                         limmadown.noNAT.NP)),
                                "consensusdownALL-NP", TRUE, 0.05)

path.down.consensus.NP<-pathwayRes(cons.data, 
                                   Reduce(intersect, list(limmadown.tcga.NP,
                                                         limmadown.R.NP,
                                                          limmadown.noNAT.NP, limmadown.noNAT.R.NP)),
                                   "consensusdownNP", TRUE, 0.05)
####M2
setwd("~/Desktop/Thesis/Data_Thesis/M2pathway")

path.tcga.DOWN.NP<-pathwayRes(dataPrep.tcga.NF, 
                                 limmadown.tcga.NP,
                                 "down-tcga-NP", TRUE, 0.05)

path.tcga.UP.NP<-pathwayRes(dataPrep.tcga.NF, 
                               limmaup.tcga.NP,
                               "up-tcga-NP", TRUE, 0.05)


path.R.noNAT.DOWN.NP<-pathwayRes(lusc.noNAT.NF, 
                              limmadown.noNAT.R.NP,
                              "down-noNAT-R-NP", TRUE, 0.05)

path.R.noNAT.UP.NP<-pathwayRes(lusc.noNAT.NF, 
                            limmaup.noNAT.R.NP,
                            "up-noNAT-R-NP", TRUE, 0.05)


path.noNAT.DOWN.NP<-pathwayRes(dataPrep.noNAT.NF, 
                            limmadown.noNAT.NP,
                            "down-noNAT-NP", TRUE, 0.05)

path.noNAT.UP.NP<-pathwayRes(dataPrep.noNAT.NF, 
                          limmaup.noNAT.NP,
                          "up-noNAT-NP", TRUE, 0.05)

path.R.DOWN.NP<-pathwayRes(lusc.all.NF, 
                                 limmadown.R.NP,
                                 "down-R-NP", TRUE, 0.05)

path.R.UP.NP<-pathwayRes(lusc.all.NF, 
                               limmaup.R.NP,
                               "up-R-NP", TRUE, 0.05)


write.csv(path.tcga.UP.NP, "up-tcgaNP.csv")
write.csv(path.tcga.DOWN.NP, "down-tcgaNP.csv")

write.csv(path.R.UP.NP, "up-RNP.csv")
write.csv(path.R.DOWN.NP, "down-RNP.csv")

write.csv(path.R.noNAT.DOWN.NP, "down-noNAT-RNP.csv")
write.csv(path.R.noNAT.UP.NP, "up-noNAT-RNP.csv")

write.csv(path.noNAT.DOWN.NP, "down-noNATNP.csv")
write.csv(path.noNAT.UP.NP, "up-noNATNP.csv")


#####Gene Ontology#########
load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/LUAD_PreprocessedData_paired_TumorPurity.rda")
my.gene.univers <- rownames(dataFilt)
up_paired_LUAD <- read.table("/data/user/marta/pipeline/DE/new_LUAD-LUSC/VennDiagram_three_methods/VennDiagram_LUAD-LUSC/final_genes_list/tumor_Purity/new_up_paired_LUAD.txt")
my.up.genes <- up_paired_LUAD$V1
logFC_file <- read.csv("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/limma/limma_LUAD_paired_tumorPurity.csv")
plot <- "circle_BP_up_paired_LUAD.png"
GOanalysis_plot("BP",my.gene.univers,my.up.genes,HUGO,20,
                "GO.up.paired.LUAD",logFC_file,plot)


###EXAMPLE######
setwd("~/Desktop/Thesis/Data_Thesis/GO10")

allgenes.tcga<-rownames(dataPrep.tcga.NF)
allgenes<-rownames(dataPrep.all2.NF)
allgenes.R<-rownames(lusc.all.NF)
allgenes.noNAT<-rownames(dataPrep.noNAT.NF)
allgenes.noNAT.R<-rownames(lusc.noNAT.NF)
allgenes.NAT<-rownames(NATvsGTEx)

####M1
GOanalysis_plot("BP", allgenes.tcga, limmadown.tcga, HUGO, 10, "GO.down.tcga", DEGs.limma, "GO.down.tcga.plot10.png")
dev.off()
GOanalysis_plot("BP", allgenes.tcga, limmaup.tcga, HUGO, 10, "GO.up.tcga", DEGs.limma, "GO.up.tcga.plot10.png")
dev.off()

GOanalysis_plot("BP", allgenes, limmadown, HUGO, 10, "GO.down.dbgap", tableDEA.dbgap.limma, "GO.down.dbgap.plot10.png")
dev.off()
GOanalysis_plot("BP", allgenes, limmaup, HUGO, 10, "GO.up.dbgap", tableDEA.dbgap.limma, "GO.up.dbgap.plot10.png")
dev.off()

GOanalysis_plot("BP", allgenes.R, limmadown.R, HUGO, 10, "GO.down.R", tableDEA.limma.R, "GO.down.R.plot10.png")
dev.off()
GOanalysis_plot("BP", allgenes.R, limmaup.R, HUGO, 10, "GO.up.R", tableDEA.limma.R, "GO.up.R.plot10.png")
dev.off()

GOanalysis_plot("BP", allgenes.noNAT, limmadown.noNAT, HUGO, 10, "GO.down.noNAT", tableDEA.limma.noNAT, "GO.down.noNAT.plot10.png")
dev.off()
GOanalysis_plot("BP", allgenes.noNAT, limmaup.noNAT, HUGO, 10, "GO.up.noNAT", tableDEA.limma.noNAT, "GO.up.noNAT.plot10.png")
dev.off()

GOanalysis_plot("BP", allgenes.noNAT.R, limmadown.noNAT.R, HUGO, 10, "GO.down.noNAT.R", tableDEA.limma.noNAT.R, "GO.down.R.noNAT.plot10.png")
dev.off()
GOanalysis_plot("BP", allgenes.noNAT.R, limmaup.noNAT.R, HUGO, 10, "GO.up.noNAT.R", tableDEA.limma.noNAT.R, "GO.up.R.noNAT.plot10.png")
dev.off()

###M2
GOanalysis_plot("BP", allgenes.noNAT, limmadown.noNAT.NP, HUGO, 10, "GO.down.noNAT.NP", tableDEA.limma.noNAT.NP, "GO.down.noNAT.NP.plot10.png")
dev.off()
GOanalysis_plot("BP", allgenes.noNAT, limmaup.noNAT.NP, HUGO, 10, "GO.up.noNAT.NP", tableDEA.limma.noNAT.NP, "GO.up.noNAT.NP.plot10.png")
dev.off()

GOanalysis_plot("BP", allgenes.noNAT.R, limmadown.noNAT.R.NP, HUGO, 10, "GO.down.noNAT.R.NP", tableDEA.limma.noNAT.NP.R, "GO.down.noNAT.R.NP.plot10.png")
dev.off()
GOanalysis_plot("BP", allgenes.noNAT.R, limmaup.noNAT.R.NP, HUGO, 10, "GO.up.noNAT.R.NP", tableDEA.limma.noNAT.NP.R, "GO.up.noNAT.R.NP.plot10.png")
dev.off()

GOanalysis_plot("BP", allgenes.R, limmadown.R.NP, HUGO, 10, "GO.down.R.NP", tableDEA.limma.R.NP, "GO.down.R.NP.plot10.png")
dev.off()
GOanalysis_plot("BP", allgenes.R, limmaup.R.NP, HUGO, 10, "GO.up.R.NP", tableDEA.limma.R.NP, "GO.up.R.NP.plot10.png")
dev.off()

GOanalysis_plot("BP", allgenes.NAT, limmadown.NAT, HUGO, 10, "GO.down.NAT", tableDEA.limma.NAT, "GO.down.NAT.plot10.png")
dev.off()
GOanalysis_plot("BP", allgenes.NAT, limmaup.NAT, HUGO, 10, "GO.up.NAT", tableDEA.limma.NAT, "GO.up.NAT.plot10.png")
dev.off()

###DONT FORGET#####

load("~/Desktop/Thesis/Data_Thesis/GO10/GO.up.R.RData")
View(res)
