

#####################\##################################################################################
################################################################################################################################
############################################################################################
## PCA functions

#' @title plotPCA
#' 
#' @description Plots a PCA and colors the dots according to the classes of a covariate
#'
#' @param cov Covariate that should be plotted
#' @param pca Pre-computed pca using the prcomp functions
#' @param var_explained Variance explained by the single principle components
#' @param covs Data frame (or SummarizedExperiment) containing the covariate that should be plotted
#' @param pos Position of the legend
#' @param col_set_factor RColorBrewer color scheme for factors
#' @param col_set_numeric RColorBrewer color scheme for numbers
#' @docType methods
#' @export
#' @importFrom graphics legend
#' @importFrom graphics plot
#' @importFrom stats prcomp
#' @importFrom stats var
#' @rdname plotPCA
plotPCA <- function(cov, pca, var_explained, covs, pos="bottomleft", col_set_factor="Set1", col_set_numeric="Blues"){
  
  #get the right number of colors
  if (is.numeric(covs[[cov]])){
    score<-covs[[cov]]
    col <- RColorBrewer::brewer.pal(9,col_set_numeric)[2:9]
    breaks <- seq(min(score),
                  max(score),
                  length = length(col) + 1 )
    breaks <- signif(breaks, 3)
    grps <- cut(score, breaks = breaks, include.lowest = TRUE)
    cols <- col[grps]
    legends <- paste0(breaks[-length(breaks)],"-",breaks[-1])
  }
  else{
    #make sure the covariate is a factor
    covs[[cov]] <- as.factor(as.character(covs[[cov]]))
    
    if (length(levels(covs[[cov]])) == 2){
      col <- RColorBrewer::brewer.pal(3,col_set_factor)[1:2]
    } else{
      col <- RColorBrewer::brewer.pal(length(levels(covs[[cov]])),col_set_factor)
    }
    cols <- col[covs[[cov]]]
    legends <- levels(covs[[cov]])
  }
  
  plot(pca$x[,1], 
       pca$x[,2], 
       main=paste("PCA -",cov), 
       xlab=paste0("Principal Component 1 (",var_explained[1],"%)"), 
       ylab=paste0("Principal Component 2 (",var_explained[2],"%)"), 
       pch=16, 
       las=1, cex=1.5, cex.axis=1.2, cex.lab=1.2, 
       col=cols)
  legend(pos, fill=col, legend=legends)
}   



#' @title autoPCAplotter
#' 
#' @description Autoplotter that plots a PCA for each selected covariate, which are those that have more than 1 class and less
#' than MAX_CLASSES_FOR_PCA. Will also plot numeric data.
#'
#'
#' @param se SummarizedExperiment
#' @param pca pca results as outputted by prcomp, by default this is NA and the function calculates the principle components based on assays(se)$final
#' @param var_explained variance explained by each principle component, by default this is NA and calculated based on the final data
#' @param MAX_CLASSES_FOR_PCA Maximum number of classes to consider
#' @param exclude An array of covariate names that should be excluded
#' @param col_set_factor RColorBrewer color scheme for factors
#' @param col_set_numeric RColorBrewer color scheme for numbers
#' @param plot_numeric Flag to indicate whether numberic covariates should be plotted as well
#' 
#' @docType methods
#' @export
#' @importFrom SummarizedExperiment colData
#' 
#' @rdname autoPCAplotter
autoPCAplotter <- function(se, 
                           pca = NA, 
                           var_explained = NA, 
                           MAX_CLASSES_FOR_PCA = 8, 
                           exclude = c(), 
                           col_set_factor = "Set1", 
                           col_set_numeric = "Blues",
                           plot_numeric = F){
  
  # Derive principle components based on the final data matrix
  if (class(class(pca)) != 'prcomp'){
    stopifnot('final' %in% names(assays(se)))
    pca <- prcomp(t(assays(se)$final))
  }
  
  # Derive the variance explained based on the PCA
  if (length(var_explained) == 1 && is.na(var_explained)){
    var_explained <- round(pca$sdev ^ 2 / sum(pca$sdev ^ 2) * 100, digits = 1)
  }
  
  #figure out the numver of classes
  num_of_classes <- apply(colData(se),
                          2,
                          function(x)
                            length(levels(as.factor(as.character(x)))))
  col_is_numeric <- sapply(colData(se),
                           is.numeric)
  
  #use only covs with 2-MAX_CLASSES_FOR_PCA levels
  covs_to_use <-  num_of_classes > 1 & num_of_classes <= MAX_CLASSES_FOR_PCA
  if (plot_numeric){
    covs_to_use <- covs_to_use | col_is_numeric
  } else {
    covs_to_use <- covs_to_use & !col_is_numeric
  }
  covs <- colData(se)[, covs_to_use]
  
  #exlude numeric columns that don't vary
  col_var <- sapply(colData(se)[,col_is_numeric],var)
  zero_variance <- names(which(col_var==0))
  covs <- covs[,!colnames(covs) %in% zero_variance]
  
  #exclude the covs that were pre-specified
  covs <- covs[,!colnames(covs) %in% exclude]
  
  #figure out the best place to put the legend
  pc <- pca$x
  prange <- apply(pc[,1:2],2,range)
  span <- prange[2,] - prange[1,]
  prange[1,] <- prange[1,] + span * 0.4
  prange[2,] <- prange[2,] - span * 0.4
  rownames(prange) <- c('1st Qu.','3rd Qu.')
  quad <- c(sum(pc[,1] < prange['1st Qu.','PC1'] & pc[,2] < prange['1st Qu.','PC2']),
            sum(pc[,1] < prange['1st Qu.','PC1'] & pc[,2] > prange['3rd Qu.','PC2']),
            sum(pc[,1] > prange['3rd Qu.','PC1'] & pc[,2] < prange['1st Qu.','PC2']),
            sum(pc[,1] > prange['3rd Qu.','PC1'] & pc[,2] > prange['3rd Qu.','PC2']))
  names(quad) <- c('bottomleft','topleft','bottomright','topright')
  pos <- names(which.min(quad))
  
  #plot all PCAs if there are any
  if (ncol(covs) > 0){
    dummy <- sapply(colnames(covs), plotPCA, pca, var_explained, covs, pos, col_set_factor, col_set_numeric)
  }
  
}
