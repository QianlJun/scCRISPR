# load bamfile
load_bam <- function(bampath, params, outpath, cores){
  start_time <- Sys.time()
  print("loading bamfile")
  print(paste("start time: ", start_time, sep = ""))

  if (is.null(cores)){
    mc.cores = data.table::getDTthreads()
  } else {
    mc.cores = cores
  }
  yieldSize(bampath) <- yieldsize
  open(bampath)
  end = 0
  while(end == 0){
    bam <- Rsamtools::scanBam(bampath, param = params)
    bam[[1]]$CB <- bam[[1]]$tag$CB
    bam[[1]]$UB <- bam[[1]]$tag$UB
    bam <- bam[[1]][-which(names(bam[[1]]) == "tag")]
    bam$seq <- as.vector(bam$seq)
    if (length(bam$qname) < yieldsize){
      end = 1
    }

    if (length(bam$qname) > 0){
      qnames <- parallel::pvec(seq_along(bam$seq), function(i){bam$qname[i][grep(bc, bam$seq[i], ignore.case = FALSE)]}, mc.cores=mc.cores)
      num <- which(bam$qname %in% qnames)
      bam <- lapply(bam, function(x){x[num]})

      if (length(bam$qname) > 0){
        data.table::fwrite(data.table::setDT(bam), outpath, sep = "\t", quote = F, col.names = F, row.names = F, append = TRUE)
      }
    }
  }

  rm(bam)
  gc()

  close(bampath)
  yieldSize(bampath) <- NA

  end_time <- Sys.time()
  print(paste("end time: ", end_time, sep = ""))
}

# call line contain sgRNA
call_matchlist <- function(data.frame){
  match_list <- vector("list", nrow(sgRNA.data))
  for (i in 1:nrow(sgRNA.data)){
    match_line <- grep(sgRNA.data$sequence[i], data.frame$seq)
    if (length(match_line) > 0){
      match_list[[i]] <- cbind(sgRNA.data[rep(i,nrow(data.frame[match_line,])),],data.frame[match_line,])
    }
  }
  return(match_list)
}

# convert to grange file
convert_to_gr <- function(bam){
  # 统一为chr格式
  chr_name <- function(bam){
    if (TRUE %in% grepl("chr", bam$rname, ignore.case = TRUE)){
      seqnames = bam$rname
    } else {
      seqnames = paste("chr",bam$rname,sep = "")
    }
    return(seqnames)
  }
  bam$pos <- as.numeric(bam$pos)
  bam$qwidth <- as.numeric(bam$qwidth)
  gr <- GenomicRanges::GRanges(seqnames = chr_name(bam),
                ranges = IRanges::IRanges(start = bam$pos, end = bam$pos-1+bam$qwidth),
                strand = bam$strand,
                qwidth = bam$qwidth,
                qname = bam$qname,
                flag = bam$flag,
                seq = bam$seq,
                CB = bam$CB,
                UB = bam$UB
  )
  return(gr)
}

# call results from unmapped reads
call_result_unmap <- function(data.frame){
  result <- data.frame(sgRNA_name = data.frame$oligo, sgRNA_sequence = data.frame$sequence, target_gene = data.frame$gene,
                       matched_name = data.frame$qname, matched_chr = NA, matched_seq = data.frame$seq,
                       matched_CB = data.frame$CB, matched_UMI = data.frame$UB, matched_annotation = NA,
                       matched_transcriptId = NA, matched_SYMBOL = NA)
}

# call results from mapped reads
call_result_map <- function(data.frame){
  result <- data.frame(sgRNA_name = data.frame$oligo, sgRNA_sequence = data.frame$sequence, target_gene = data.frame$gene,
                       matched_name = data.frame$qname,
                       matched_chr = paste(data.frame$seqnames, paste(data.frame$start, data.frame$end, sep = "-"),sep = ":"),
                       matched_seq = data.frame$seq,
                       matched_CB = data.frame$CB, matched_UMI = data.frame$UB, matched_annotation = data.frame$annotation,
                       matched_transcriptId = data.frame$transcriptId, matched_SYMBOL = data.frame$SYMBOL)
}

# decide sgRNA in each cell
sgRNA_count <- function(unmap.result,map.result){
  UMI_rm <- c()
  for (i in 1:nrow(unmap.result)){
    CB <- unmap.result[i,]$matched_CB
    UB <- unmap.result[i,]$matched_UMI
    if (UB %in% map.result$matched_UMI[map.result$matched_CB == CB]){
      UMI_rm <- c(UMI_rm, i)
    }
  }

  if (length(UMI_rm) > 0){
    unmap.result <- unmap.result[-UMI_rm,]
    sgRNA_counts <- filter_result(unmap.result)
  } else {
    sgRNA_counts <- filter_result(unmap.result)
  }

  cells <- unique(sgRNA_counts$matched_CB)
  cell_sgRNA_list <- vector("list", length(cells))
  for (i in 1:length(cells)){
    cell.df <- sgRNA_counts[sgRNA_counts$matched_CB == cells[i], ]
    sgRNAs <- unique(cell.df$sgRNA_sequence)
    sgRNA_list <- vector("list", length(sgRNAs))
    for (j in 1:length(sgRNAs)){
      sgRNA.df <- cell.df[cell.df$sgRNA_sequence == sgRNAs[j], ]
      read_counts <- length(sgRNA.df$matched_UMI)
      umi_counts <- length(unique(sgRNA.df$matched_UMI))
      sgRNA_list[[j]] <- data.frame(cell_barcode = cells[i], sgRNA_name = unique(sgRNA.df$sgRNA_name),
                                    sgRNA_sequence = unique(sgRNA.df$sgRNA_sequence), target_gene = unique(sgRNA.df$target_gene),
                                    sgRNA_readcounts = read_counts, sgRNA_umicounts = umi_counts)
    }
    cell_sgRNA_list[[i]] <- do.call(rbind.data.frame, sgRNA_list)
  }
  results <- do.call(rbind.data.frame, cell_sgRNA_list)

  return(results)
}

# decide sgRNA in each cell, if no mapped reads got
sgRNA_count_unmap <- function(unmap.result){

  unmap.result <- filter_result(unmap.result)

  cells <- unique(unmap.result$matched_CB)
  cell_sgRNA_list <- vector("list", length(cells))
  for (i in 1:length(cells)){
    cell.df <- unmap.result[unmap.result$matched_CB == cells[i], ]
    sgRNAs <- unique(cell.df$sgRNA_sequence)
    sgRNA_list <- vector("list", length(sgRNAs))
    for (j in 1:length(sgRNAs)){
      sgRNA.df <- cell.df[cell.df$sgRNA_sequence == sgRNAs[j], ]
      read_counts <- length(sgRNA.df$matched_UMI)
      umi_counts <- length(unique(sgRNA.df$matched_UMI))
      sgRNA_list[[j]] <- data.frame(cell_barcode = cells[i], sgRNA_name = unique(sgRNA.df$sgRNA_name),
                                    sgRNA_sequence = unique(sgRNA.df$sgRNA_sequence), target_gene = unique(sgRNA.df$target_gene),
                                    sgRNA_readcounts = read_counts, sgRNA_umicounts = umi_counts)
    }
    cell_sgRNA_list[[i]] <- do.call(rbind.data.frame, sgRNA_list)
  }
  results <- do.call(rbind.data.frame, cell_sgRNA_list)

  return(results)
}

# remove sequence in ummapped reads which can match sgRNA but without anchor sequence from vector
filter_result <- function(unmap.result){

  preseq <- c()
  behseq <- c()
  if (nrow(unmap.result) > 1000){
    for (i in 1:1000){
      sg <- unmap.result$sgRNA_sequence[i]
      seq <- unmap.result$matched_seq[i]
      pos <- regexpr(sg, seq)
      len <- nchar(seq)
      if (pos-4 >= 0 & pos+4 <= len){
        preseq <- c(preseq, substr(seq, pos-4, pos-1))
        behseq <- c(behseq, substr(seq, pos+nchar(sg), pos+nchar(sg)+3))
      } else if (pos-4 >= 0 & pos+4 > len) {
        preseq <- c(preseq, substr(seq, pos-4, pos-1))
      } else if (pos-4 < 0 & pos+4 <= len) {
        behseq <- c(behseq, substr(seq, pos+nchar(sg), pos+nchar(sg)+3))
      }
    }
  } else {
    for (i in 1:(nrow(unmap.result))){
      sg <- unmap.result$sgRNA_sequence[i]
      seq <- unmap.result$matched_seq[i]
      pos <- regexpr(sg, seq)
      len <- nchar(seq)
      if (pos-4 >= 0 & pos+4 <= len){
        preseq <- c(preseq, substr(seq, pos-4, pos-1))
        behseq <- c(behseq, substr(seq, pos+nchar(sg), pos+nchar(sg)+3))
      } else if (pos-4 >= 0 & pos+4 > len) {
        preseq <- c(preseq, substr(seq, pos-4, pos-1))
      } else if (pos-4 < 0 & pos+4 <= len) {
        behseq <- c(behseq, substr(seq, pos+nchar(sg), pos+nchar(sg)+3))
      }
    }
  }

  pre.anchor <- names(table(preseq))[which.max(table(preseq))]
  beh.anchor <- names(table(behseq))[which.max(table(behseq))]

  keep <- c()
  for (j in 1:nrow(unmap.result)){
    sg <- unmap.result$sgRNA_sequence[j]
    seq <- unmap.result$matched_seq[j]
    pos <- regexpr(sg, seq)
    len <- nchar(seq)
    if (pos-4 >= 0 & pos+4 <= len){
      preseq <- substr(seq, pos-4, pos-1)
      behseq <- substr(seq, pos+nchar(sg), pos+nchar(sg)+3)
      if (preseq == pre.anchor | behseq == beh.anchor){
        keep <- c(keep, j)
      }
    } else if (pos-4 >= 0 & pos+4 > len) {
      preseq <- substr(seq, pos-4, pos-1)
      if (preseq == pre.anchor){
        keep <- c(keep, j)
      }
    } else if (pos-4 < 0 & pos+4 <= len) {
      behseq <- substr(seq, pos+nchar(sg), pos+nchar(sg)+3)
      if (behseq == beh.anchor){
        keep <- c(keep, j)
      }
    }
  }
  unmap.result <- unmap.result[keep,]
  return(unmap.result)
}

sgRNA_detect <- function(bamPath, sgPath, output = NULL, species, refgenome,
                         cores = NULL, search_mapped = FALSE, yieldsize = 10000000){
  if (is.null(output)){
    output = paste(getwd(), "/output/sgRNA_detect", sep = "")
    dir.create(output)
  } else {
    output = paste(output, "/output/sgRNA_detect", sep = "")
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
  sgRNA <- paste(sgRNA.data$sequence, collapse = "|")
  bc <- toupper(sgRNA)
  bamFile <- Rsamtools::BamFile(bamPath)

  params.unmap <- Rsamtools::ScanBamParam(what = c("qname","flag","rname","strand","pos","qwidth","seq"), tag = c("CB","UB"), flag = scanBamFlag(isUnmappedQuery = TRUE))
  tmp1 <- paste(output,"/tmp1.txt", sep = "")
  load_bam(bamFile, params.unmap, tmp1, cores)

  if (file.exists(tmp1)){
    bam.unmap <- data.table::fread(tmp1, sep = "\t", data.table = FALSE, showProgress = FALSE, na.strings = NULL, skip = "-")
    colnames(bam.unmap) <- c("qname","flag","rname","strand","pos","qwidth","seq","CB","UB")
    num.filt <- which(bam.unmap$CB == ""|bam.unmap$UB == "")
    if (length(num.filt) > 0){
      bam.unmap <- bam.unmap[-num.filt,]
    }
    CB <- unique(bam.unmap$CB)

    unmap.match <- do.call(rbind.data.frame, call_matchlist(bam.unmap))
    unmap.result <- call_result_unmap(unmap.match)

    rm(bam.unmap, unmap.match)
    gc()

  } else {
    print("no unmapped reads found")
  }

  if (search_mapped){
    params.map <- Rsamtools::ScanBamParam(what = c("qname","flag","rname","strand","pos","qwidth","seq"), tag = c("CB","UB"), flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE), tagFilter = list(CB = CB))
    tmp2 <- paste(output,"/tmp2.txt", sep = "")
    load_bam(bamFile, params.map, tmp2)
  } else {
    tmp2 <- paste(output,"/tmp2.txt", sep = "")
  }

  if (file.exists(tmp2)){
    bam.map <- data.table::fread(tmp2, sep = "\t", data.table = FALSE, showProgress = FALSE, na.strings = NULL)
    colnames(bam.map) <- c("qname","flag","rname","strand","pos","qwidth","seq","CB","UB")
    # check if CB or UB has NA
    num.filt <- which(bam.map$CB == ""|bam.map$UB == "")
    if (length(num.filt) > 0){
      bam.map <- bam.map[-num.filt,]
    }

    # gene Annotation
    gr.map.keep <- convert_to_gr(bam.map)
    if (species == "mouse"){
      anno.map.keep <- ChIPseeker::annotatePeak(gr.map.keep, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Mm.eg.db")
      anno.map.keep <- as.data.frame(anno.map.keep)
      anno.map.keep.match <- do.call(rbind.data.frame, call_matchlist(anno.map.keep))
      map.result <- call_result_map(anno.map.keep.match)
    } else if (species == "human"){
      anno.map.keep <- ChIPseeker::annotatePeak(gr.map.keep, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
      anno.map.keep <- as.data.frame(anno.map.keep)
      anno.map.keep.match <- do.call(rbind.data.frame, call_matchlist(anno.map.keep))
      map.result <- call_result_map(anno.map.keep.match)
    } else {
      print("only human/mouse genome support gene annotation now")
      map.result <- bam.map
    }
    write.table(map.result, paste(output, "/sgRNA_in_mapped_reads.txt", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)

    rm(bam.map, anno.map.keep, anno.map.keep.match)
    gc()

  } else {
    print("no mapped reads found")
  }

  if (file.exists(tmp1) & file.exists(tmp2)){
    sgRNA_counts <- sgRNA_count(unmap.result,map.result)
    print("running successfully!")
  } else if (file.exists(tmp1) & file.exists(tmp2) == FALSE){
    print("running successfully!")
    print("no sgRNA detected in mapped reads!")
    sgRNA_counts <- sgRNA_count_unmap(unmap.result)
  } else if (file.exists(tmp1) == FALSE & file.exists(tmp2)){
    print("no cells detected reliable sgRNA!")
    sgRNA_counts <- data.frame()
  } else if (file.exists(tmp1) == FALSE & file.exists(tmp2) == FALSE){
    print("no cells detected sgRNA!")
    sgRNA_counts <- data.frame()
  }
  write.table(sgRNA_counts,paste(output, "/sgRNA_counts.txt", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
}




