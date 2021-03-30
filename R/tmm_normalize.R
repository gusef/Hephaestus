
#' @title tmm_normalize
#' 
#' @description Standard function to normalize a RNAseq cube
#'
#' @param se SummarizedExperiment 
#' 
#' @return SummarizedExperiment that contains the normalized data
#' 
#' @docType methods
#' @export
#' 
#' @importFrom edgeR DGEList
#' @importFrom edgeR calcNormFactors
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assays<-
#' 
#' @rdname tmm_normalize
tmm_normalize <- function(se){
  
  #Using edgeR to extract the normalization factors
  dge <- DGEList(counts = assays(se)$raw_counts)
  dge <- calcNormFactors(dge)
  se$sampleNormFactors <- dge$samples$norm.factors
  
  #adjust using the normalization factors
  assays(se)$log2_tmm_cpm <- sweep(assays(se)$raw_counts,2,se$sampleNormFactors,'/')
  
  #counts per million
  assays(se)$log2_tmm_cpm <- sweep(assays(se)$log2_tmm_cpm,2,1e6,'*')
  assays(se)$log2_tmm_cpm <- sweep(assays(se)$log2_tmm_cpm,2,dge$samples$lib.size,'/')
  
  #bring into log2 space (using 1 because log2 results in 0)
  assays(se)$log2_tmm_cpm <- log2(assays(se)$log2_tmm_cpm + 1)
  
  #there should be one assay name that points to the data after the last processing step no matter if we use Liver Confounding
  assays(se)$final <- assays(se)$log2_tmm_cpm   
  
  return(se)
  
}