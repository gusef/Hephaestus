
#' @title extract_GRanges
#' 
#' @description Generates a GenomicRanges object from a .gff file
#'
#' @param gff_file Location of the .gff file
#' 
#' @return GRanges object
#' @docType methods
#' @export
#' 
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges values<-
#' 
#' @rdname extract_GRanges
extract_GRanges <- function(gff_file){
  
  #grab the gff file and remove all comment lines
  annot <- readLines(gff_file)
  annot <- annot[-grep('##',annot)]
  annot <- do.call(rbind,strsplit(annot,'\t'))
  
  #use only the genes
  annot <- annot[annot[,3] == 'gene',]
  
  #only use genes that have an actual chromosome assigned to them
  annot <- annot[grep('chr',annot[,1]),]
  
  #split the rest
  gene_names <- annot[,9]
  annot <- annot[,-9]
  colnames(annot) <- c('chr','source','feature','start','end','score','strand','frame')
  
  gene_names <- lapply(gene_names,function(x){
    split <- strsplit(x,';')[[1]]
    split <- strsplit(split,'=')
    if(!all(sapply(split,length) == 2)){
      stop('There were more than one "=" in the string so splitting will lead to errors.')
    }
    entries <- sapply(split,function(x)x[2])
    names(entries) <- sapply(split,function(x)x[1])
    
    #split off the IDs as well
    ids <- entries[2]
    entries <- entries[-2]
    ids <- strsplit(ids,',')[[1]]
    names(ids) <- sub(':.+$','',ids)
    ids <- sub('^[^:]+:','',ids)
    names(ids)[1] <- 'GeneID'
    entries <- c(ids,entries)
    return(entries)
  })
  
  #the entries vary from entry to entry so we need to get all possible identifiers to put them into a matrix
  all_identifiers <- unique(unlist(sapply(gene_names,names)))
  gene_names <- t(sapply(gene_names,function(x)x[all_identifiers]))
  colnames(gene_names) <- all_identifiers
  
  #there are some entries like exception, start range, end range that are mostly NA with a few exceptions
  #removing those
  gene_names <- gene_names[,colSums(is.na(gene_names)) < 10000]
  
  #combine all of them
  annot <- cbind(annot,gene_names)
  annot <- data.frame(annot,stringsAsFactors = F)
  
  #removing the duplicates
  dups <- annot$GeneID[duplicated(annot$GeneID)]
  dups <- annot[annot$GeneID %in% dups,]
  dups <- dups[order(dups$GeneID),]
  
  #remove the duplicate genes on the Y chromosome that are also on the X chr
  annot <- annot[!(annot$GeneID %in% dups$GeneID & annot$chr == 'chrY'),]
  
  #now we can use the gene names as row identifiers since there are no duplicates
  rownames(annot) <- annot$Name
  
  #Build a genomic ranges object
  gr <- makeGRangesFromDataFrame(annot) 
  stopifnot(all(names(gr) == rownames(annot)))
  values(gr) <- annot[,9:ncol(annot)]
  
  return(gr)
}


#' Extract gtf attribute
#' @description Extract specific feature of gtf attribute by attribute name.
#' 
#' @param gtf_attributes the gtf attribute separated by ";"
#' @param att_of_interest name of the attributes (e.g. gene_name, gene_biotype)
#' @return attribute of interest in vector format
#' 
#' @export
#' 
#' @rdname extract_attributes
#' 
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

#' @title extract_GRanges_from_gtf
#' 
#' @description Generates a GenomicRanges object from a .gtf file
#'
#' @param gtf_file Location of the .gtf file
#' 
#' @return GRanges object
#' @docType methods
#' @export
#' 
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges values<-
#' 
#' @rdname extract_GRanges_from_gtf
#' 
extract_GRanges_from_gtf <- function(gtf_file){
  
  #grab the gff file and remove all comment lines
  annot <- readLines(gtf_file)
  annot <- annot[-grep('##',annot)]
  annot <- do.call(rbind,strsplit(annot,'\t'))
  
  #use only the genes
  annot <- annot[annot[,3] == 'gene',]
  
  #only use genes that have an actual chromosome assigned to them
  annot <- annot[grep('chr',annot[,1]),]
  
  #split the rest
  attributes <- annot[,9]
  annot <- annot[,-9]
  colnames(annot) <- c('chr','source','feature','start','end','score','strand','frame')
  
  gene_names <- unlist(lapply(attributes, extract_attributes, "gene_name"))
  gene_type <- unlist(lapply(attributes, extract_attributes, "gene_type"))
  ensembl_id <- unlist(lapply(attributes, extract_attributes, "gene_id"))
  
  #the entries vary from entry to entry so we need to get all possible identifiers to put them into a matrix
  # all_identifiers <- unique(unlist(sapply(gene_names,names)))
  # gene_names <- t(sapply(gene_names,function(x)x[all_identifiers]))
  # colnames(gene_names) <- all_identifiers
  
  #there are some entries like exception, start range, end range that are mostly NA with a few exceptions
  #removing those
  # gene_names <- gene_names[,colSums(is.na(gene_names)) < 10000]
  
  #combine all of them
  annot <- cbind(annot,gene_names, gene_type, ensembl_id)
  annot <- data.frame(annot,stringsAsFactors = F)
  
  #removing the duplicates
  dups <- annot$gene_names[duplicated(annot$gene_names)]
  dups <- annot[annot$gene_names %in% dups,]
  dups <- dups[order(dups$gene_names),]
  
  #remove the duplicate genes on the Y chromosome that are also on the X chr
  annot <- annot[!(annot$gene_names %in% dups$gene_names & annot$chr == 'chrY'),]
  
  #removing additional duplicates
  dups <- annot$gene_names[duplicated(annot$gene_names)]
  dups <- annot[annot$gene_names %in% dups,]
  dups <- dups[order(dups$gene_names),]
  annot <- annot[!(annot$gene_names %in% dups$gene_names),]
  
  #now we can use the gene names as row identifiers since there are no duplicates
  rownames(annot) <- annot$gene_names
  
  #Build a genomic ranges object
  gr <- makeGRangesFromDataFrame(annot) 
  stopifnot(all(names(gr) == rownames(annot)))
  values(gr) <- annot[,9:ncol(annot)]
  
  return(gr)
}