# convert to grange file
convert_to_gr <- function(data.frame){
  # unify as chr format
  chr_name <- function(data.frame){
    if (TRUE %in% grepl("chr", data.frame$match_chr, ignore.case = TRUE)){
      seqnames = data.frame$match_chr
    } else {
      seqnames = paste("chr", data.frame$match_chr,sep = "")
    }
    return(seqnames)
  }

  data.frame$match_pos <- as.numeric(data.frame$match_pos)
  data.frame$match_width <- as.numeric(data.frame$match_width)
  gr <- GenomicRanges::GRanges(seqnames = chr_name(data.frame),
                ranges = IRanges::IRanges(start = data.frame$match_pos, end = data.frame$match_pos-1+data.frame$match_width),
                strand = data.frame$match_strand,
                pos = data.frame$match_pos,
                match_width = data.frame$match_width,
                match_seq = data.frame$match_seq,
                sgRNA_name = data.frame$sgRNA_name, sgRNA_sequence = data.frame$sgRNA_sequence, target_gene = data.frame$target_gene
  )
  return(gr)
}

sgRNA_search <- function(sgRNA.data, fasta){
  if (nrow(sgRNA.data) > 0){
    pb <- utils::txtProgressBar(style=3)

    result_list <- vector("list", length = nrow(sgRNA.data))
    for (i in 1:nrow(sgRNA.data)){

      sgname <- sgRNA.data$oligo[i]
      sgRNA <- sgRNA.data$sequence[i]
      target_gene <- sgRNA.data$gene[i]
      chr_pos <- grep("chr",names(fasta))
      chr_name <- names(fasta)[grep("chr",names(fasta))]
      chr_list <- vector("list", length = length(chr_pos)*2)
      # search in 5'-> 3' direction
      for (j in 1:length(chr_pos)){
        chr <- chr_name[j]
        sg_pos <- unlist(gregexpr(sgRNA, fasta[[chr_pos[j]]][1]))
        strand <- "+"
        if (sg_pos != -1){
          pos_list <- vector("list", length = length(sg_pos))
          for (k in 1:length(sg_pos)){
            pos_list[[k]] <- data.frame(sgRNA_name = sgname, sgRNA_sequence = sgRNA, target_gene = target_gene,
                                        match_chr = chr, match_pos = sg_pos[k], match_width = nchar(sgRNA), match_strand = strand,
                                        match_seq = substr(fasta[[chr_pos[j]]][1], sg_pos[k], sg_pos[k]+nchar(sgRNA)+2))
          }
          chr_list[[j]] <- do.call(rbind.data.frame, pos_list)
        }

      }
      # search in 3'->5' direction
      for (j in 1:length(chr_pos)){
        chr <- chr_name[j]
        sg_pos <- unlist(gregexpr(stri_reverse(mgsub(sgRNA,c("A","T","G","C","a","t","g","c","N","n"),
                                                     c("T","A","C","G","t","a","c","g","N","n"))), fasta[[chr_pos[j]]][1]))
        strand <- "-"
        if (sg_pos != -1){
          pos_list <- vector("list", length = length(sg_pos))
          for (k in 1:length(sg_pos)){
            pos_list[[k]] <- data.frame(sgRNA_name = sgname, sgRNA_sequence = sgRNA, target_gene = target_gene,
                                        match_chr = chr, match_pos = sg_pos[k], match_width = nchar(sgRNA), match_strand = strand,
                                        match_seq = substr(fasta[[chr_pos[j]]][1], sg_pos[k], sg_pos[k]+nchar(sgRNA)+2))
          }
          chr_list[[j+length(chr_pos)]] <- do.call(rbind.data.frame, pos_list)
        }
      }
      result_list[[i]] <- do.call(rbind.data.frame, chr_list)

      utils::setTxtProgressBar(pb, i/nrow(sgRNA.data))
    }
    result <- do.call(rbind.data.frame, result_list)
  } else {
    result_list <- list()
    result_list[[1]] <- data.frame(sgRNA_name = NA, sgRNA_sequence = NA, target_gene = NA,match_chr = NA, match_pos = NA, match_width = NA, match_strand = NA, match_seq = NA)
  }
  return(result_list)
}

sgRNA_annotation_org <- function(sgRNA.data, fasta, species, refgenome, txdb){

  sglist <- vector("list", length = data.table::getDTthreads())
  cutnum <- ceiling(nrow(sgRNA.data)/28)
  for (m in 1:data.table::getDTthreads()){
    sglist[[m]] <- sgRNA.data[((m-1)*cutnum+1):(m*cutnum),][!is.na(sgRNA.data[((m-1)*cutnum+1):(m*cutnum),][,1]),]
  }
  result <- do.call(rbind.data.frame,
                    parallel::pvec(seq_along(sglist), function(i){sgRNA_search(sglist[[i]], fasta)}, mc.cores=data.table::getDTthreads()))
  result <- result[!is.na(result$sgRNA_name),]
  #result <- sgRNA_search(sgRNA.data, fasta)

  result.gr <- convert_to_gr(result)
  if (species == "mouse" || species == "human"){
    if (species == "mouse"){
      anno.result <- annotatePeak(result.gr, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Mm.eg.db")
    } else if (species == "human"){
      anno.result <- annotatePeak(result.gr, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
    }
    anno.result <- as.data.frame(anno.result)
    anno.result <- data.frame(sgRNA_name = anno.result$sgRNA_name, sgRNA_sequence = anno.result$sgRNA_sequence, target_gene = anno.result$target_gene,
                              match_chr = anno.result$seqnames, match_pos = anno.result$pos, match_width = anno.result$match_width,
                              match_strand = anno.result$strand, match_seq = anno.result$match_seq,annotation = anno.result$annotation,
                              geneChr = anno.result$geneChr, geneStart = anno.result$geneStart, geneEnd = anno.result$geneEnd,
                              geneStrand = anno.result$geneStrand, geneId = anno.result$geneId, ENSEMBL = anno.result$ENSEMBL, SYMBOL = anno.result$SYMBOL)

    anno.result$offtarget = c(anno.result$target_gene != anno.result$SYMBOL)

    sgname <- sgRNA.data$oligo
    if (FALSE %in% (sgname %in% anno.result$sgRNA_name)){
      sgname.unfind <- sgRNA.data[-which(sgname %in% anno.result$sgRNA_name),]
      unfind.df <- data.frame(sgRNA_name = sgname.unfind$oligo, sgRNA_sequence = sgname.unfind$sequence, target_gene = sgname.unfind$gene,
                              match_chr = NA, match_pos = NA, match_width = NA,
                              match_strand = NA, match_seq = NA, annotation = NA,
                              geneChr = NA, geneStart = NA, geneEnd = NA,
                              geneStrand = NA, geneId = NA, ENSEMBL = NA, SYMBOL = NA, offtarget = NA)
      anno.result <- rbind(anno.result, unfind.df)
    }
  } else {
    print("only human/mouse support gene annotation")
    anno.result <- c("no results")
  }
  anno.result$no_match <- is.na(anno.result$SYMBOL)
  anno.result$match_seq <- paste(anno.result$match_seq, "(5'->3')", sep = "")
  anno.result$genome <- refgenome

  if (species == "human"){
    if (refgenome == "hg19" || refgenome == "GRCh37"){
      anno.result$source <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
    } else if (refgenome == "hg38" || refgenome == "GRCh38"){
      anno.result$source <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    }
  } else if (species == "mouse"){
    if (refgenome == "mm9"){
      anno.result$source <- "TxDb.Mmusculus.UCSC.mm9.knownGene"
    } else if (refgenome == "mm10"){
      anno.result$source <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
    }
  }

  return(anno.result)
}

sgRNA_annotate <- function(sgPath, fasta, species, refgenome, output = NULL){
  if (is.null(output)){
    output = paste(getwd(), "/output/sgRNA_annotate", sep = "")
    dir.create(output)
  } else {
    output = paste(output, "/output/sgRNA_annotate", sep = "")
    dir.create(output)
  }

  if (species == "human"){
    if (refgenome == "hg19" || refgenome == "GRCh37"){
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if (refgenome == "hg38" || refgenome == "GRCh38"){
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }
  } else if (species == "mouse"){
    if (refgenome == "mm9"){
      txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
    } else if (refgenome == "mm10"){
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    }
  }

  sgRNA.data <- read.table(sgPath,header = T)
  fasta <- read.fasta(fasta, as.string = TRUE, forceDNAtolower = FALSE)
  sgRNA.anno <- sgRNA_annotation_org(sgRNA.data, fasta, species, refgenome, txdb)
  write.table(sgRNA.anno, paste(output,"/sgRNA_annotation.txt", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
}


