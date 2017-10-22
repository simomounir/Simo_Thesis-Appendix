 #######FUNCTIONS and libraries#########################
 
 setwd("~/Desktop/Thesis/Data_Thesis")
 library(TCGAbiolinks)
 library(SummarizedExperiment)
 library(DT)
 library(edgeR)
 library(limma)
 library(devtools)
 library(roxygen2)
 library(sva)
 library(biomaRt)
 library(ggplot2)
 library(cqn)
 library(Biobase)
 library(e1071)
 library(dendextend)
 library(plyr)
 library(RColorBrewer)
 library(gplots)
 library(rafalib)
 library(VennDiagram)
 library(ggfortify)
 library(ReactomePA)
 library(org.Hs.eg.db)
 library(topGO)
 library(GOplot)
 library(grid)
 library(gridExtra)
 library(gPCA)
 library(reshape2)
 
 

####Thilde's function for MDS plots######
Heatmap.full<-function(my.data,my.group){
  ####Clustering########
  # By default calculates the distance between rows
  my.data.L = log2(my.data + 1)
  colnames(my.data.L)<-1:ncol(my.data.L)
  dist.L = dist(t(my.data.L), method='man')
  
  ## Look at distance matrix
  colramp = colorRampPalette(c(3,"white",2))(9)
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  cols <- palette(brewer.pal(8, "Dark2"))[as.numeric(as.factor(my.group))]
  #heatmap(as.matrix(dist.L),col=colramp,Colv=NA,Rowv=NA)
  hclust.L = hclust(dist.L, method = "complete")
  
  ###dendextend package
  dend <- as.dendrogram(hclust.L)
  dend<- rotate(dend, 1:ncol(my.data.L))
  dend <- color_branches(dend, k=3)
  
  # We shall add the flower type to the labels:
  labels(dend) <- my.group
  # We hang the dendrogram a bit:
  dend <- hang.dendrogram(dend,hang_height=0.1)
  dend <- set(dend, "labels_cex", 0.3)
  
  rv <- rowVars(as.matrix(my.data.L))
  idx <- order(-rv)[1:1000]
  
  heatmap.2(as.matrix(my.data.L[idx,]), labCol=as.numeric(as.factor(my.group)),
            trace="none", 
            ColSideColors=cols, 
            col=hmcol)
  
  legend("topright", legend = levels(as.factor(my.group)), col = cols[factor(levels(as.factor(my.group)))])
  
  some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
  
  gplots::heatmap.2(as.matrix(t(dataPrep.all2.L)), 
                    main = "Heatmap for GTEX&TCGA dataset",
                    srtCol = 20,
                    dendrogram = "row",
                    Rowv = dend,
                    Colv = "NA", # this to make sure the columns are not ordered
                    trace="none",          
                    margins =c(5,0.1),      
                    key.xlab = "Cm",
                    denscol = "grey",
                    density.info = "density",
                    RowSideColors = rev(labels_colors(dend)), # to add nice colored strips        
                    col = some_col_func)
  
  
  
}
 
 
PC.1D<-function(my.data, my.group){
  
  s <- svd(my.data)
  count.data<-as.data.frame(table(my.group))
  rownames(count.data)<-count.data$my.group
  
  unique.batches<-which(my.group%in%rownames(count.data[count.data$Freq==1,])==FALSE)
  
  s2<-s$v[unique.batches, unique.batches]
  variable <- as.factor(my.group)
  variable2 <- as.factor(my.group[unique.batches])
  
  mypar(2,1)
  d0<-s$d
  #plot(d0^2/sum(d0^2),ylim=c(0,.25))
  #mypar(2,1)
  
  boxplot(split(s2[,1],variable2),las=3,range=0, col=1, xaxt = "n")
  stripchart(split(s2[,1],variable2),add=TRUE,vertical=TRUE,pch=2,cex=.5,col=1)
  
  boxplot(split(s$v[,1],variable),las=3,range=0, col=1, xaxt = "n")
  stripchart(split(s$v[,1],variable),add=TRUE,vertical=TRUE,pch=2,cex=.5,col=1)
  
}

#Fstats 
Aov.corr<-function(my.data, my.group){
  s<-svd(my.data)
  corr <- sapply(1:ncol(s$v),function(i){
    fit <- lm(s$v[,i]~as.factor(my.group))
    return( summary(fit)$adj.r.squared  )
  })
  mypar()
  plot(seq(along=corr), corr, xlab="PC")
  
  
}



Fstats<-function(my.data, my.group){
  s <- svd(my.data)
  
  Fstats.res<- sapply(1:ncol(s$v),function(i){
    fit <- lm(s$v[,i]~as.factor(my.group))
    Fstat <- summary(aov(fit))[[1]][1,4]
    return(Fstat)
  })
  print(Fstats.res[1])
  mypar()
  plot(seq(along=Fstats.res),sqrt(Fstats.res))
  p <- length(unique(my.group))
  abline(h=sqrt(qf(0.995,p-1,ncol(s$v)-1))) 
}
 
UniqueDEGs<-function(l1,l2){
   return(c(length(setdiff(l1,l2)), length(setdiff(l2,l1))))
 }
 
myPCA<-function(my.data, my.group, title="PCs"){
   
   data.pca<-prcomp(my.data, center=TRUE)
   my.data<-as.data.frame(my.data)
   my.data$grouping<-my.group
   PCi<-data.frame(data.pca$x[,1:9],grouping=my.group)
   PCi.melt<-melt(PCi, id.vars=c('grouping'))
   
   plot.melted<- ggplot(PCi.melt, aes(x = grouping, y = value, fill = grouping)) +
     geom_boxplot() +
     scale_fill_manual(values = palette(rainbow(length(unique(my.group)))) ) +
     facet_wrap(~variable, scales="free_y") + ggtitle(title) 
   
   
   my.plot<-ggplot(PCi,aes(x=PC1,y=PC2,col=grouping))+
     geom_point(size=3,alpha=0.5)+
     #scale_color_manual(values = palette(rainbow(length(unique(my.group)))))+ #your colors here
     theme_classic()+ 
     xlab(paste0("PC1 ",summary(data.pca)$importance[2,1]*100, "%"))+
     ylab(paste0("PC2 ",summary(data.pca)$importance[2,2]*100, "%"))+
     ggtitle(title)
   
   P <- list()
   vars <- names(PCi)
   
   for (var in names(PCi[1:9])){
     PCi2<-PCi[,c(var,'grouping')]
     names(PCi2)[names(PCi2) == var] <- "PC"
     p <- ggplot(PCi2, aes(x = grouping, y = PC, fill = grouping)) +
       geom_boxplot() +
       scale_fill_manual(values = palette(rainbow(length(unique(my.group)))) ) +
       ggtitle(var)
     P<-c(P, list(p))
   }
   
   #return(ggplot(aes(y = PC1, x = grouping, fill=grouping), data = PCi) + geom_boxplot()+ ggtitle("PC1"))
   return(list(plots=P, PC2=my.plot, melted=plot.melted))
 } 
 
 

 
myMDSplot <- function(my.data, my.group, my.labels) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  p <- ggplot(data=res)
  p + geom_point(aes(x=M1,y=M2,color=my.group)) + geom_text(aes(x=M1,y=M2, label= my.labels, color=my.group)) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_gray() + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold"), axis.text.x=element_blank()) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) + theme(axis.text = element_text(colour = "black"))
}



###filter low counts genes with rowsums or avg
Filter_genes<-function(SE, avg, cpm.thresh, row_sums, groupvec, n=1){
  
  DGE.SE<-DGEList(SE, group=groupvec)
  #DGE.SE.norm<-calcNormFactors(DGE.SE)
  
  cpm<-cpm(DGE.SE, log=FALSE)
  logcpm.norm<-cpm(DGE.SE.norm, log=TRUE, scale=TRUE)
  
  keep <-rowSums(cpm>cpm.thresh) >= row_sums
  keep.exprs.avg<-aveLogCPM(DGE.SE) >avg
  ###keep.exprs.zero<-
  
  if(n==1){
    keep <-rowSums(cpm>cpm.thresh) >= row_sums
  }
  
  else{
    keep<-aveLogCPM(DGE.SE) >avg
  }
  
  return(DGE.SE[keep,,keep.lib.sizes=FALSE])
}

####DEA Limma and EdgeR#####
DE_limma <- function(my.contrast, my.data, my.design, coLFC, coFDR) {
  fit3 <- eBayes(contrasts.fit(lmFit(my.data, my.design), my.contrast))
  tt <- toptable(fit3, coef=1, adjust='fdr', number=nrow(my.data))
  index.up <- which(tt$logFC >= coLFC & tt$adj.P.Val < coFDR)
  index.down <- which(tt$logFC <= -coLFC & tt$adj.P.Val < coFDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  direction[!(1:nrow(tt) %in% union(index.up,index.down))] <- "no DE"
  tt <- cbind(tt,direction)
  return(tt)
}

DE_edgeR <- function(my.lrt, my.data, coLFC, coFDR) {
  my.tags <- topTags(my.lrt, n=nrow(my.data$counts))
  my.tags <- my.tags$table
  
  tt <- cbind(my.lrt$table, FDR = p.adjust(my.lrt$table$PValue, "fdr"))
  index.up <- which(tt$logFC >= coLFC & tt$FDR < coFDR)
  index.down <- which(tt$logFC <= -coLFC & tt$FDR < coFDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  direction[!(1:nrow(tt) %in% union(index.up,index.down))] <- "no DE"
  tt <- cbind(tt,direction)
  
  return(tt)
}

get_IDs <- function(data) {
  IDs <- strsplit(c(colnames(data)), "-")
  IDs <- ldply(IDs, rbind)
  colnames(IDs) <- c('project', 'tss','participant', 'sample', "portion", "plate", "center")
  cols <- c("project", "tss", "participant")
  IDs$patient <- apply(IDs[,cols],1,paste,collapse = "-" )
  barcode <- colnames(data)
  IDs <- cbind(IDs, barcode)
  condition <- gsub("11+[[:alpha:]]", "normal", as.character(IDs$sample))
  condition  <- gsub("01+[[:alpha:]]", "cancer", condition)
  IDs$condition <- condition
  IDs$myorder  <- 1:nrow(IDs)
  return(IDs)
}

overlap<-function(v1,v2){
  l1<-length(v1)
  l2<-length(v2)
  l.i<-length(intersect(v1,v2))
  return(l.i/(l1+l2-l.i))
}

####function for cqn normalization
condit_QN<-function(Data, groupvec){
  geneInfoDF<-as.data.frame(geneInfo)
  Data<-Data[(rownames(Data)%in%rownames(geneInfoDF)),]
  
  
  remove<-rownames(geneInfoDF[which(is.na(geneInfoDF$gcContent)),])
  Data<-Data[!rownames(Data) %in% remove,]
  remove<-rownames(geneInfoDF[which(is.na(geneInfoDF$geneLength)),])
  Data<-Data[!rownames(Data) %in% remove,]
  
  gc<-as.numeric(as.character(geneInfoDF[rownames(Data),]$gcContent))
  geneLength<-as.numeric(geneInfoDF[rownames(Data),]$geneLength)
  
  DGE<-DGEList(Data, group=groupvec)
  lib.size<-DGE$samples$lib.size
  print(dim(Data))
  Data.cqn<-cqn(counts = Data, 
                x = gc, 
                lengths = geneLength,
                sizeFactors = lib.size, 
                verbose = TRUE)
  return(Data.cqn)
}

#test
#condition.tcga<-get_IDs(dataPrep.tcga.rsem.copy)$condition
#dataPrep.tcga.rsem.cqn<-condit_QN(dataPrep.tcga.rsem.copy, condition.tcga)

edgeR_DEA<-function(Data, design, groupvec, lfc, fdr, ref){
  groupvec<-factor(x=groupvec, levels=unique(groupvec))
  colnames(design)[1:length(levels(groupvec))]<-levels(groupvec)
  
  aDGEList <- edgeR::DGEList(counts = Data, group = groupvec)
  aDGEList$samples$condition<-as.factor(groupvec)
  aDGEList$samples$condition<-relevel(aDGEList$samples$condition, ref=ref)
  
  aDGEList <- edgeR::estimateGLMCommonDisp(aDGEList, design)
  aDGEList <- edgeR::estimateGLMTagwiseDisp(aDGEList, design)
  aGlmFit <- edgeR::glmFit(aDGEList, design, dispersion = aDGEList$tagwise.dispersion,
                           prior.count.total=0)
  
  my.lrt <- edgeR::glmLRT(aGlmFit, coef = 2)
  return(DE_edgeR(my.lrt, Data, lfc, fdr))
  
}

VennPlot<-function(a, annot){
  grid.newpage()
  
  if(length(a)==2){
    venn.plot <- venn.diagram(list(a[[1]], a[[2]]), NULL,
                              fill=c("red", "green"), 
                              alpha=c(0.5,0.5), cex = 2, cat.fontface=4,
                              category.names=annot)
    
  }
  
  else if(length(a)==3){
    venn.plot <- venn.diagram(list(a[[1]], a[[2]], a[[3]]), NULL,
                              fill=c("red", "green", "blue"), 
                              alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, 
                              category.names=annot)
    
  }
  
  else if(length(a)==4){
    venn.plot <- venn.diagram(list(a[[1]], a[[2]], a[[3]], a[[4]]), NULL, 
                              fill=c("red", "green", "blue", "orange"), 
                              alpha=c(0.5,0.5,0.5,0.5), cex = 2, cat.fontface=4,
                              category.names=annot)
  
  }
  
  else if(length(a)==5){
    venn.plot <- venn.diagram(list(a[[1]], a[[2]], a[[3]], a[[4]], a[[5]]), NULL,
                              fill=c("red", "green", "blue", "orange", 'cyan'),
                              alpha=c(0.5,0.5,0.5,0.5, 0.5), cex = 2, cat.fontface=4, 
                              category.names=annot)
  }
  else if(length(a)==6){
    venn.plot <- venn.diagram(list(a[[1]], a[[2]], a[[3]], a[[4]], a[[5]], a[[6]]), NULL,
                              fill=c("red", "green", "blue", "orange", 'cyan', "purple"),
                              alpha=c(0.5,0.5,0.5,0.5,0.5,0.5), cex = 2, cat.fontface=4, 
                              category.names=annot)   
  }
  grid.draw(venn.plot)
  return(venn.plot)
}

enrich_pathway <- function(background, my.genes, my.name, my.plot, value=0.05) {
  entrez.data <- background[background$geneName %in% my.genes, ]
  my.pathway <- enrichPathway(as.character(entrez.data$entrez), organism = "human", pvalueCutoff = value, pAdjustMethod = "fdr", qvalueCutoff = value, as.character(background$entrez), minGSSize = 2, readable = T)
  if (my.plot == TRUE) {
    png(filename = paste0(my.name,".png"), height = 4500, width = 4500)
    cnetplot(my.pathway, categorySize="pvalue", showCategory = 8)
    dev.off()
  }
  return(my.pathway)
} 

conversion.background <- function(dataframe){
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  total_genes <- rownames(dataframe)
  
  #covert the gene names (universal genes) into entrez IDs
  background <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'hgnc_symbol',
                      values = total_genes , mart = ensembl)
  colnames(background) <- c("geneName","entrez")
  return (background)
}

pathwayRes<-function(my.data, my.genes, name, my.plot, cutoff){
  convert.genes<-conversion.background(my.data)
  return(enrich_pathway(convert.genes, my.genes, name, my.plot, cutoff))
}

libsize.corr<-function(my.data){
  y<-calcNormFactors(DGEList(counts = as.matrix(my.data)))
  eff.y<-y$samples$norm.factors*y$samples$lib.size
 return(y$counts/eff.y)
}

OverlapMatrix<-function(DEGslist, R){
  o.m<-matrix(nrow=length(DEGslist), ncol=length(DEGslist))
  u.m<-matrix(nrow=length(DEGslist), ncol=length(DEGslist))
  
  rownames(o.m)<-R
  colnames(o.m)<-R
  rownames(u.m)<-R
  colnames(u.m)<-R
    
  n<-1
  m<-1
  for(i in DEGslist){
    m<-1
    for(j in DEGslist){
      o.m[n,m]=overlap(i,j)
      temp<-as.character(UniqueDEGs(i,j))
      u.m[n,m]=paste0(temp[1], ",", temp[2])
      m<-m+1
    }
    n<-n+1
  }
  
  return(list(u=u.m, o=o.m))
}

