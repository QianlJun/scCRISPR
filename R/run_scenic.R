run_scenic <- function(filepath, species, refgenome, dbDir, myDatasetTitle = "Mydata",
                       output = NULL, nontarget, method = "GENIE3", seed = NULL, interest.genes = NULL, cores = NULL){
  if (is.null(output)){
    output = paste(getwd(), "/output/run_scenic", sep = "")
    dir.create(output)
  } else {
    output = paste(output, "/output/run_scenic", sep = "")
    dir.create(output)
  }

  prepath <- getwd()
  setwd(output)

  if (species == "human"){
    org <- "hgnc"
    if (refgenome == "hg19" || refgenome == "GRCh37"){
      dbs <- c('500bp'= 'hg19-500bp-upstream-7species.mc9nr.feather', '10kb' = 'hg19-tss-centered-10kb-7species.mc9nr.feather')
    } else if (refgenome == "hg38" || refgenome == "GRCh38"){
      dbs <- c('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
    }
  } else if (species == "mouse"){
    org <- "mgi"
    if (refgenome == "mm9"){
      dbs <- c('500bp'= 'mm9-500bp-upstream-7species.mc9nr.feather', '10kb' = 'mm9-tss-centered-10kb-7species.mc9nr.feather')
    } else if (refgenome == "mm10"){
      dbs <- c('500bp'= 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', '10kb' = 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
    }
  }

  if (is.null(cores)){
    ncores = 1
  } else {
    ncores = cores
  }

  # prepare to run scenic
  seurat.obj <- readRDS(filepath)
  exprMat <- as.matrix(seurat.obj@assays$RNA@counts)
  cellInfo <- seurat.obj@meta.data
  scenicOptions <- SCENIC::initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=ncores)
  saveRDS(cellInfo, file="int/cellInfo.rds")
  scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

  if (!is.null(seed)){
    scenicOptions@settings$seed <- seed
  } else {
    seed <- 201
  }

  genesKept <- SCENIC::geneFiltering(exprMat, scenicOptions=scenicOptions,
                             minCountsPerGene=3*.01*ncol(exprMat),
                             minSamples=ncol(exprMat)*.01)

  # check whether target genes and interesting genes filtered
  target.genes <- unique(seurat.obj$target_gene)
  target.genes <- target.genes[-which(target.genes == nontarget)]

  if (is.null(interest.genes)){
    interestingGenes <- target.genes
  } else {
    interestingGenes <- c(target.genes, interest.genes)
  }

  filtered.genes <- interestingGenes[!interestingGenes %in% genesKept]
  write.table(filtered.genes, paste(output,"/int/interestingGenes filtered by geneFiltering.txt", sep = ""),quote = F,col.names = F,row.names = F)

  genesKept <- unique(c(interestingGenes, genesKept))
  saveRDS(genesKept, "int/1.1_genesKept.Rds")

  exprMat_filtered <- exprMat[genesKept, ]
  SCENIC::runCorrelation(exprMat_filtered, scenicOptions)
  exprMat_filtered <- log2(exprMat_filtered+1)


  # except TFs in database, we add some interesting genes, such as target genes as candidate TFs to build GRN
  candidate.TFs <- unique(c(SCENIC::getDbTfs(scenicOptions), interestingGenes))

  if (method == "GENIE3"){
    SCENIC::runGenie3(exprMat_filtered, scenicOptions, allTFs = candidate.TFs)
  } else if (method == "GRNBoost2"){
    SCENIC::exportsForArboreto(exprMat_filtered, scenicOptions, dir = "int")
    run_grnboost2(output, cores, seed)
    reticulate::py_run_file(paste(getwd(),"/.run_grnboost2.py",sep=""))
    grnboost.op <- read.table(paste(output,"/int/output_grnboost.tsv", sep=""),header = F,sep = "\t")
    colnames(grnboost.op) <- c("TF", "Target", "weight")
    saveRDS(grnboost.op, paste(output,"/int/1.4_GENIE3_linkList.Rds", sep=""))
  } else {
    print("Only GENIE3 or GRNBoost2 support now.")
  }

  exprMat_log <- log2(exprMat+1)

  scenicOptions <- SCENIC::runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- SCENIC::runSCENIC_2_createRegulons(scenicOptions)


  scenicOptions <- SCENIC::initializeScenic(org=org, dbDir=dbDir, dbs=mm10_dbs, datasetTitle=myDatasetTitle, nCores=1)


  # check if interestingGenes exist or not in regulons
  regulons <- readRDS("int/2.6_regulons_asGeneSet.Rds")
  regulons <- regulons[lengths(regulons)>=10]
  TFs <- unique(c(names(regulons)[grep("_extended", names(regulons), invert = TRUE)],
                  gsub("_extended","",names(regulons)[grep("_extended", names(regulons))])))
  Tfs.add <- interestingGenes[!interestingGenes %in% TFs]
  write.table(Tfs.add, paste(output,"/int/interestingGenes not in regulons_asGeneSet.txt", sep = ""),quote = F,col.names = F,row.names = F)

  cormat <- readRDS("int/1.2_corrMat.Rds")
  weightMatrix <- readRDS("int/1.4_GENIE3_linkList.Rds")
  weightMatrix <- weightMatrix[weightMatrix$weight > 0.001,]
  weightMatrix <- weightMatrix[weightMatrix$TF %in% Tfs.add,]
  Tfs.list <- vector("list", length(Tfs.add))
  names(Tfs.list) <- Tfs.add
  for (i in 1:length(Tfs.add)){
    tf <- Tfs.add[i]
    targets <- as.vector(weightMatrix[weightMatrix$TF == tf,]$Target)
    posGenes <- names(cormat[tf,])[cormat[tf,] > 0.03]
    postargets <- intersect(targets, posGenes)
    Tfs.list[[i]] <- postargets
  }
  regulons <- c(regulons, Tfs.list)
  saveRDS(regulons, "int/2.6_regulons_asGeneSet.Rds")

  # pick target gene regulons for score separately
  tag.list <- vector("list", length = length(target.genes))
  names(tag.list) <- target.genes
  for (j in 1:length(target.genes)){
    tf <- target.genes[j]
    if (tf %in% names(regulons)[grep("_extended", names(regulons), invert = TRUE)]){
      tag.list[[j]] <- regulons[[tf]]
    } else if (paste(tf,"_extended",sep = "") %in% names(regulons)[grep("_extended", names(regulons), invert = TRUE)]){
      tag.list[[j]] <- regulons[[paste(tf,"_extended",sep = "")]]
    }
  }
  cells_rankings <- AUCell::AUCell_buildRankings(exprMat_log, nCores=1, plotStats=TRUE)
  cells_AUC.regulon <- AUCell::AUCell_calcAUC(tag.list, cells_rankings, nCores = 1)
  saveRDS(cells_AUC.regulon, "int/scoreTargetGeneregulon.rds")


  # filter regulon < 10 genes first, then score regulon
  scenicOptions <- SCENIC::runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
  scenicOptions@fileNames$output["loomFile",] <- "output/SCENIC.loom"
  #saveRDS(scenicOptions, file="int/scenicOptions.Rds")
  scenicOptions <- SCENIC::runSCENIC_4_aucell_binarize(scenicOptions)
  saveRDS(scenicOptions, file="int/scenicOptions.Rds")
  SCENIC::export2loom(scenicOptions, exprMat)

  # score regulons(targets genes) separately

  # delete tmp file
  if (file.exists(paste(getwd(),"/.run_grnboost2.py",sep=""))){
    file.remove(paste(getwd(),"/.run_grnboost2.py",sep=""))
  }

  setwd(prepath)
}
