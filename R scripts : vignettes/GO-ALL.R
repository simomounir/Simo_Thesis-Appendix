#-----------------------------------------------------------------------------
#ont= "BP" or "CC", "MF"
#nameGO: topGO output name
#logFC_file: file that includes logFC (for example from limma analysis)
GOanalysis_plot <- function(ont,my.gene.univers,my.DE.genes,HUGO,nTerm,nameGO,logFC_file,plot){
  
  GO <- TOPGO(ont, my.gene.univers, my.DE.genes, HUGO, nTerm, nameGO)
  vector <- make_vector(my.DE.genes,GO,HUGO)
  genes_interest <- make_genes_interest(my.DE.genes,GO,HUGO)
  terms <- make_terms_dataframe(GO,vector,ont)
  genes <- make_genes_dataframe(logFC_file,genes_interest)
  circ <- circle_dat(terms, genes)
  process <- terms$term
  chord <- chord_dat(data = circ, genes, process)
  png(filename=plot,height = 1000, width = 1500)
  GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(genes$logFC), lfc.max = max(genes$logFC))
  #dev.off()
  
  #return(dev.off())
}

#-------------------------------------------------------------------------------------------
# LUAD paired
#----------------------------------------------------------------------------------------
#up

