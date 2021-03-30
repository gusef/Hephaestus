#' @title deriveLiverScore
#' 
#' @description
#' Derive a liver score for each sample based on the 22 liver gene signature
#'
#' @param dat Gene expression matrix - rows are genes and columns are samples
#' @param method Method used to estimate the liver score ('mean' or 'gsva') 
#' @return Vector of liver scores indicating the levels of liver confounding
#' @docType methods
#' @export
#' @importFrom GSVA gsva
#' @importFrom stats var
#' @importFrom utils data
#' 
#' @rdname deriveLiverScore
deriveLiverScore <- function(dat, method=c('mean','gsva')){
  
  #read the signature file stored in the package
  liver_signature <- readRDS(system.file("extdata", 'liver_signature.RDS', package = "Swiffer"))
  
  if(!all(liver_signature %in% rownames(dat))){
    stop('Not all genes in the liver signature are present in the rownames of the provided dataset.')
  }
  
  if (method[1] == 'mean'){
    liver_score <- colMeans(dat[liver_signature,])
    
    
  } else if (method[1] == 'gsva'){
    dat <- dat[apply(dat,1,var) > 0,]
    liver_score <- gsva(dat, 
                        list(liver=liver_signature),
                        method="gsva")
    
    #in GSVA < 1.26 the return value was a lis
    if (class(liver_score) == 'list'){
      liver_score <- liver_score$es.obs
    }
    
    #reshape the data
    liver_score <- liver_score[1,]
  }else {
    stop("Method can only be 'mean' or 'gsva'.")  
  }   
  
  return(liver_score)
}

#' @title Liver score PCA
#' 
#' @description
#' Plot a PCA and superimpose the liver score levels on the samples
#' @param dat Gene expression matrix - rows are genes and columns are samples
#' @param method Method used to estimate the liver score ('mean' or 'gsva')
#' @param se_matrix indicates which matrix to use in case a SE is passed. By default this is final, but can usually also be log2_tmm_cpm
#' @return Vector of gene scores indicating the levels of liver confounding
#' @docType methods
#' @export
#' @importFrom graphics plot
#' @importFrom stats prcomp
#' @rdname plotLiverScore
plotLiverScore <- function(dat, method=c('mean','gsva'), se_matrix = 'final'){
  
  #there are two modes - either just using a matrix or a SE
  if (class(dat) == 'RangedSummarizedExperiment'){
    dat <- assays(dat)[[se_matrix]]
  }
  
  #derive score
  score <- deriveLiverScore(dat, method)
  
  #break it into color bins
  col <- RColorBrewer::brewer.pal(11,"RdBu")[11:1]
  breaks <- seq(min(score),
                max(score),
                length = length(col) + 1 )
  grps <- cut(score, breaks = breaks, include.lowest = TRUE)
  col <- col[grps]
  
  #run a PCA
  pca <- prcomp(t(dat))
  var_explained <- round(pca$sdev ^ 2 / sum(pca$sdev ^ 2) * 100,digits = 1)
  
  #plot the score on the PCA
  plot(pca$x[,1], 
       pca$x[,2], 
       main='Liver confounding',
       xlab=paste0("Principal Component 1 (",var_explained[1],"%)"), 
       ylab=paste0("Principal Component 2 (",var_explained[2],"%)"), 
       pch=16, las=1, cex=1.5, cex.axis=1.2, cex.lab=1.2, 
       col=col)
  
  return(score)
}


#' @title Liver confounder removal
#' 
#' @description
#' Removes a liver confounder in a gene expression dataset
#'
#' @param dat Gene expression matrix - rows are genes and columns are samples.
#' @param plot Flag indicating whether a PCA should be plotted.
#' @param method Method used to estimate the liver score ('mean' or 'gsva')
#' @param se_matrix indicates which matrix to use in case a SE is passed. By default this is final, but can usually also be log2_tmm_cpm
#' 
#' @return Vector of gene scores indicating the levels of liver confounding
#' @docType methods
#' @export
#' 
#' @importFrom stats lm
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assays<-
#' 
#' @rdname removeLiverConfounder
removeLiverConfounder <- function(dat, plot = T, method=c('mean','gsva'), se_matrix = 'final'){
  
  using_se <- F
  
  #there are two modes - either just using a matrix or a SE
  if (class(dat) == 'RangedSummarizedExperiment'){
    using_se <- T
    se <- dat
    dat <- assays(dat)[[se_matrix]]
  }
  
  if (plot){
    liver_score <- plotLiverScore(dat, method = method)
  } else {
    liver_score <- deriveLiverScore(dat, method = method)
  }
  
  #remove the confounder
  dat <- t(apply(dat,
                 1,
                 function(y,x){
                   lm1 <- lm(y ~ x)
                   return(lm1$residuals + lm1$coefficients['(Intercept)'])
                 },
                 liver_score))
  
  #if a SE was passed it, we pass a SE back
  if (using_se){
    se$Liver.scores <- liver_score
    assays(se)$final <- dat   
    dat <- se
  }
  
  return(dat)
}