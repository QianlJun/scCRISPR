run_GSEA <- function(markers, species, output, plot = FALSE, topTerm = 5){
  if (is.null(output)){
    output = paste(getwd(), "/output/run_GSEA", sep = "")
    dir.create(output)
  } else {
    dir.create(output)
  }

  if (species == "human"){
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  } else if (species == "mouse"){
    OrgDb <- org.Mm.eg.db::org.Mm.eg.db
  } else {
    print("only human or mouse supported now!")
  }

  markers.genes <- rownames(markers)
  gene.df <- bitr(genes, fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb=OrgDb)
  df <- merge(x = gene.df, y = markers, by = "ENSEMBL")
  GSEA.input <- df$avg_log2FC
  names(GSEA.input) <- as.character(df$ENTREZID)
  GSEA.input <- sort(GSEA.input, decreasing = TRUE)
  gseGO.res <- gseGO(GSEA.input, OrgDb = OrgDb, ont = 'BP', pvalueCutoff = 0.05, keyType = 'ENTREZID')
  write.table(as.data.frame(gseGO.res), paste(output, "/gseGO.csv", sep = ""), col.names = T, row.names = T, quote = F,sep = ",")

  gseKEGG.res <- gseKEGG(GSEA.input, pvalueCutoff = 0.05, keyType = 'kegg', use_internal_data = FALSE, minGSSize = 5)
  write.table(as.data.frame(gseKEGG.res), paste(output, "/gseKEGG.csv", sep = ""), col.names = T, row.names = T, quote = F,sep = ",")

  if (plot){
    # default figure 1
    pdf(paste(output, "/gseGO_topterms.pdf", sep = ""), width = 10, height = 8)
    ridgeplot(gseGO.res, topTerm)
    dev.off()

    # default figure 2
    for (i in 1:topTerm){
      pdf(paste(output, "/gseGO_top", i, "terms.pdf", sep = ""), width = 10, height = 8)
      gseaplot2(gseGO.res, i, title = gseGO.res$Description[i])
      dev.off()
    }
  }
}
