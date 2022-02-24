sgRNA_funceval <- function(filepath, nontarget = NULL, output){
  if (is.null(output)){
    output = paste(getwd(), "/output/sgRNA_funceval", sep = "")
    dir.create(output)
  } else {
    output = paste(output, "/output/sgRNA_funceval", sep = "")
    dir.create(output)
  }

  if (is.null(nontarget)){
    print("Please allocate the name of non-targert gene use parameter:nontarget")
  }

  seurat.obj <- readRDS(filepath)

  TargetGeneRegulons <- readRDS("int/scoreTargetGeneregulon.rds")
  TargetGeneRegulons <- TargetGeneRegulons@assays@data$AUC

  nontarget.cells <- colnames(seurat.obj)[seurat.obj$target_gene == nontarget]

  for (i in rownames(TargetGeneRegulons)){
    sgRNAs <- unique(seurat.obj$sgRNA_name[seurat.obj$target_gene == i])
    nontarget.df <- data.frame(sgRNA = "non_target", regulon_score = TargetGeneRegulons[i,nontarget.cells])

    for (j in sgRNAs){
      sgRNA.cells <- colnames(seurat.obj)[seurat.obj$sgRNA_name == j]
      sgRNA.df <- data.frame(sgRNA = j, regulon_score = TargetGeneRegulons[i,sgRNA.cells])
      df <- rbind(sgRNA.df,nontarget.df)
      p1 <- ggplot(df, aes(x=regulon_score, color=sgRNA)) + geom_density() + theme_bw()
      ggsave(paste(output, "/",j,".png",sep = ""),p1,width = 6,height = 4)

      em <- normalmixEM(sgRNA.df$regulon_score, k = 2)
      em.df <- data.frame(lambda = em$lambda, mu = em$mu, sigma = em$sigma)
      write.table(em.df, paste(output, "/",j,"_mixmodel.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

      p2 <- plot(em, whichplots = 2)
      ggsave(paste(output, "/",j,"_DensityCurves.png",sep = ""),p1,width = 6,height = 4)
    }

    sgRNA.cells <- colnames(seurat.obj)[seurat.obj$sgRNA_name %in% sgRNAs]
    sgRNA.df <- data.frame(targetgene = i, regulon_score = TargetGeneRegulons[i,sgRNA.cells])
    df <- rbind(sgRNA.df,nontarget.df)
    p1 <- ggplot(df, aes(x=regulon_score, color=sgRNA)) + geom_density() + theme_bw()
    ggsave(paste(output, "/",i,"_all.png",sep = ""),p1,width = 6,height = 4)

    em <- normalmixEM(sgRNA.df$regulon_score, k = 2)
    em.df <- data.frame(lambda = em$lambda, mu = em$mu, sigma = em$sigma)
    write.table(em.df, paste(output, "/",i,"_all_mixmodel.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

    p2 <- plot(em, whichplots = 2)
    ggsave(paste(output, "/",i,"_all_DensityCurves.png",sep = ""),p1,width = 6,height = 4)
  }

}
