run_KEGG <- function(genes, type = "SYMBOL", species, output, plot = FALSE, topTerm = 20){
  if (is.null(output)){
    output = paste(getwd(), "/output/run_KEGG", sep = "")
    dir.create(output)
  } else {
    dir.create(output)
  }

  if (species == "human"){
    organism = "hsa"
  } else if (species == "mouse"){
    organism = "mmu"
  } else {
    print("only human or mouse supported now!")
  }

  if (type == "SYMBOL"){
    gene.df <- clusterProfiler::bitr(genes, fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb=OrgDb)
  } else if (type == "ENSEMBL"){
    gene.df <- clusterProfiler::bitr(genes, fromType="ENSEMBL",toType=c("SYMBOL","ENTREZID"),OrgDb=OrgDb)
  } else if (type == "ENTREZID"){
    gene.df <- clusterProfiler::bitr(genes, fromType="ENTREZID",toType=c("SYMBOL","ENSEMBL"),OrgDb=OrgDb)
  }

  kegg <- clusterProfiler::enrichKEGG(gene.df$ENTREZID, organism = organism, keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = 'BH', use_internal_data = FALSE)
  KEGG_result <- kegg@result
  for (j in 1:nrow(KEGG_result)){
    geneID <- KEGG_result$geneID[j]
    KEGG_result$symbol[j] <- paste(up.df$SYMBOL[up.df$ENTREZID %in% strsplit(geneID,"/")[[1]]], collapse = "/")
  }

  write.table(KEGG_result,paste(ouput, "/KEGG.csv",sep=""),sep = ",",quote = F,row.names = F,col.names = T)

  if (plot){
    # default figure
    pdf(paste(output, "/KEGG_topterms_type1.pdf", sep = ""), width = 5, height = 8)
    clusterProfiler::dotplot(kegg, showCategory=topTerm)
    dev.off()

    # custom figure
    Terms <- KEGG_result[1:topTerm,]
    Terms$logpva <- -log10(Terms$p.adjust)
    Terms$Description <- factor(Terms$Description, levels = rev(Terms$Description))
    p7 <- ggplot2::ggplot(Terms, aes(x = Description, y = logpva)) + ggplot2::geom_bar(stat = 'identity', position = 'dodge', width = 0.8, fill="#2166AC") + ggplot2::coord_flip()
    p8 <- p7 + ggplot2::theme(axis.line = element_line(colour = "black")) + ggplot2::theme(panel.background=element_blank()) + ggplot2::scale_y_continuous(expand = c(0,0),limits = c(0, 5.5))
    p9 <- p8 + ggplot2::xlab(NULL) + ggplot2::ylab("-Log10(p.adjust)") + ggplot2::labs(title = "")
    p10 <- p9 + ggplot2::theme(axis.text.x = element_text(colour = "black", size = 20), axis.text.y = element_text(colour = "black", size = 20), axis.title = element_text(size = 20), plot.title = element_text(size = 20))
    pdf(paste(output, "/KEGG_topterms_type2.pdf", sep = ""), width = 10, height = 8)
    p10
    dev.off()
  }
}
