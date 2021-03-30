#complex  heatmap_wrapper
#' @title hcopt
#' 
#' @description Ward clustering with the optimal leaf ordering algorithm
#' 
#' @param d a dissimilarity structure as produced by dist.
#' @param HC a hierarchical clustering object (default: NULL)
#' @param method the agglomeration method to be used.
#' @param members NULL or a vector with length size of d. See the 'Details' section fo the hclust function.
#' 
#' @export
#' 
hcopt <- function(d, HC=NULL, method = "ward.D", members = NULL)
{
  if ( is.null(HC) ) {
    HC <- stats::hclust(d,method=method,members=members)
  }
  if (ncol(d) > 2){
    ORD <- cba::order.optimal(d,merge=HC$merge)
    HC$merge <- ORD$merge
    HC$order <- ORD$order
  }
  HC
}


#' @title complex_heatmap
#' 
#' @description Complex heatmap function that shows a barplot on top, can show several covariates and
#' can display 3 panels of data - IHC, gene expression and gene sets
#' 
#' @param se SummarizedExperiment cube that contains all the data
#' @param selector logic vector with length of samples that should be included. This is mainly used in conjunction with GSVA as gene set projection method, because it works better with a larger set of samples
#' @param genes Genes to be displayed, which need to be contained in the SE cube
#' @param signatures List of gene vectors indicating different signatures / pathways / gene sets
#' @param ihc Array of names in the colData data.frame of the SE that include IHC
#' @param cluster_rows logic flag indicating whether the rows of the three data panels should be clustered. This is done using hierarchical clustering by Ward's method and using the optimal leaf ordering algorithm as implemented in the cba package.
#' @param cluster_cols logic flag indicating whether the samples should be clustered. Is ignored if a barplot is specificed and is automatically set to TRUE if not specified. If you want to impose a custom order order the matrix beforehand and set it to FALSE.
#' @param covs Array of names in the colData data.frame of the SE indicating additional covariates that should be displayed on top of the three panels
#' @param annotation_colors (default: NA) List of colors for the covariates that are shown on top. 
#' @param column_order Indicating the covariate (name of the colData column) that the samples should be ordered by. This numberic value is used for the top barplot.
#' @param barplot_color This can be either NA - which colors all bars the color specified in barplot_legend, or a factorial where the levels need to correspond to the names of the barplot legend.
#' @param barplot_legend This can be a single color indicating the colors of all bars or a color scheme, which is an array of colors with names that correspond to the levels indicated in barplot_color
#' @param pathway_projection Pathway projection method. 'means' just averages the gene expression values of the passed gene sets, 'gsva' uses the R/Bioconductor package GSVA (using a Gaussian kernel in log2 space)
#' @param rowscale Logic flag indicating whether each row should be row scaled
#' @param rowscale_cap Caps the values in all three panels to +/- the cap value (default: 3), which avoids that the color schemes are driven by outliers. A values of 3 indicates that only values within 3 standard deviations are shown.
#' @param barplot_decreasing_order Logic flag that indicates the sorting order of the barplot. Can be NA - no sorting
#' @param barplot_specific_order Instead of sorting a specific ordering can be passed in case there are multiple indications or datasets that should be grouped (This overrides the barplot_decreasing_order flag)
#' @param barplot_censor Censor vector if the barplot is used for survival data
#' @param barplot_height Height of the barplot in grub small px (default: 30)
#' @param labels_col Custom sample labels. By default these are the column labels of the passed SE
#' @param gaps_col array of indices that indicates where a gap in the columns should be placed
#' @param custom_legend Changes the labels for the color scales
#' @param fontsize Used fontsize in the heatmap
#' @param pathways_name_cap Pathway labels are be capped at 40  
#' @param main Title for the complex heatmap
#' 
#' @export
#' 
#' @importFrom grDevices gray.colors
#' @importFrom stats as.dendrogram
#' @importFrom stats as.dist
#' @importFrom stats cor
#' 
complex_heatmap <- function(se, 
                            selector = NULL,
                            genes, 
                            signatures = NULL, 
                            ihc = NULL, 
                            cluster_rows = TRUE,
                            cluster_cols = NULL,
                            covs = NULL, 
                            annotation_colors = NA,
                            column_order = NULL, 
                            barplot_color = NA,
                            barplot_legend = 'grey',
                            pathway_projection = 'means',
                            rowscale = TRUE,
                            rowscale_cap = 3,
                            barplot_decreasing_order = TRUE,
                            barplot_specific_order = NULL,
                            barplot_censor = NA,
                            barplot_height = 30,
                            labels_col = NULL,
                            gaps_col = NULL,
                            custom_legend = NULL,
                            fontsize = 8,
                            pathways_name_cap = 40, 
                            main){
  
  #select only the single genes of interest
  mat <- assays(se[genes,])$final
  if (cluster_rows && nrow(mat) > 2){
    mat <- mat[order.dendrogram(as.dendrogram(hcopt(as.dist(1-cor(t(mat))), method = "ward.D"))),]
  }
  
  if (is.null(labels_col)){
    labels_col <- colnames(mat)
  }
  
  #ihc values      
  if (!is.null(ihc)){
    ihc <- colData(se)[, ihc]
    ihc <- t(as.matrix(ihc))
    mat <- rbind(log2(ihc + 1), mat)
    
  }
  
  #row scale
  if (rowscale){
    mat <- rowScale(mat)
    #cap the rowscale values
    if (!is.na(rowscale_cap)){
      mat[mat > rowscale_cap] <- rowscale_cap
      mat[mat < -rowscale_cap] <- -rowscale_cap
    }
  }
  
  #pathway projection
  if (!is.null(signatures)){
    
    #pathway projection based on GSVA
    if (pathway_projection == 'GSVA'){
      pathways <- gsva(assays(se)$final, signatures)
      pathways <- pathways * 3
      
      #pathway projection based on means of log2 values
    } else if (pathway_projection == 'means'){
      if (!all(sapply(signatures,function(x,y)all(x%in%rownames(y)),assays(se)$final))){
        warning("Not all genes in the signatures are represented within the RNAseq data. Filtering them now.")
        signatures <- lapply(signatures, function(x,y) x[x%in%y],rownames(assays(se)$final))
      }
      
      pathways <- t(sapply(signatures,
                           function(gs,exp)
                             colMeans(exp[gs,]),
                           assays(se)$final))
      #rowscaling
      if (rowscale){
        pathways <- rowScale(pathways)
        if (!is.na(rowscale_cap)){
          pathways[pathways > rowscale_cap] <- rowscale_cap
          pathways[pathways < -rowscale_cap] <- -rowscale_cap
        }
      }         
      
    } else {
      stop('Pathway projection has to be "means" or "GSVA"')
    }      
    
    #if rows should be clustered
    if (cluster_rows && nrow(pathways) > 2){
      pathways <- pathways[order.dendrogram(as.dendrogram(hcopt(as.dist(1-cor(t(pathways))), method = "ward.D"))),]
    }
    
    if (!is.na(pathways_name_cap )){
      rownames(pathways) <- strtrim(rownames(pathways), pathways_name_cap)
    }
    
    #combine the sets
    mat <- rbind(mat, pathways)
    
  }
  
  #filter the samples
  if (!is.null(selector)){
    if (!is.null(barplot_specific_order)){
      barplot_specific_order <- match(barplot_specific_order, (1:ncol(mat))[selector])
      barplot_specific_order <- barplot_specific_order[!is.na(barplot_specific_order)]
    }
    
    se <- se[, selector]
    mat <- mat[, selector]
    labels_col <- labels_col[selector]
  }
  
  #drop all samples that have NA values in the column order covariate, unless a specific ordering was specified
  if (!is.null(column_order)){
    if (is.null(barplot_specific_order)){
      filter <- !is.na(se[[column_order]])
      mat <- mat[,filter]
      se <- se[,filter]
      labels_col <- labels_col[filter]
    }
    cluster_cols <- as.numeric(se[[column_order]])
  } else if (is.null(cluster_cols)){
    cluster_cols <- TRUE
  }
  
  #gaps and colors
  color <- list()
  gaps <- c()
  offset <- 0
  
  #IHC
  if (!is.null(ihc)){
    offset <- nrow(ihc)
    gaps <- c(gaps, offset)
    color$IHC <- list(index = 1:nrow(ihc),
                      color = gray.colors(11, start = 0.2, end = 1)[11:1])
  }
  
  #genes
  color$`Gene Expression` <- list(index = (offset + 1) : (offset + length(genes)),
                                  color = brewer.pal(11, "RdBu")[11:1])
  offset <- offset + length(genes)
  
  
  #pathways
  if (!is.null(signatures)){
    gaps <- c(gaps, offset)
    color$Signatures <- list(index = (offset + 1) : nrow(mat),
                             color = brewer.pal(11, "RdYlBu")[11:1])
  }
  
  #covariates   
  if (is.null(covs)){
    col_annot <- NA
  } else if (class(covs) == 'list'){
    col_annot <- covs
  } else if (length(covs) == 1){
    col_annot <- data.frame(colData(se)[,covs])
    rownames(col_annot) <- colnames(mat)
    colnames(col_annot) <- covs
  } else {
    col_annot <- data.frame(colData(se)[,covs])
  }
  
  #set the colors for the barplots
  if (is.na(barplot_color[1])){
    barplot_color <- rep(barplot_legend, ncol(mat))
    barplot_legend <- NULL
    
  } else{
    barplot_color <- barplot_legend[match(se[[barplot_color]],names(barplot_legend))]
  }
  
  #Change the color scale labels to the custom labels
  if (!is.null(custom_legend)){
    if(length(custom_legend) != length(color)){
      stop("Error: custom_legend needs to have the same length as the number of panels. 3 - IHC/GExp/GS, 2 - GExp/GS, ...")
    }
    names(color) <- custom_legend
  }
  
  #heatmap
  g_heatmap(mat,
            na_col = 'white',
            main = main,
            scale = 'none',
            color = color,
            annotation_col = col_annot,
            annotation_colors = annotation_colors,
            cluster_rows = FALSE,
            cluster_cols = cluster_cols,
            barplot_decreasing_order = barplot_decreasing_order,
            barplot_specific_order = barplot_specific_order,
            barplot_label = column_order,
            barplot_color = barplot_color,
            barplot_legend = barplot_legend, 
            barplot_height = 80,
            barplot_censor = barplot_censor,
            na_as_X = TRUE,
            gaps_row = gaps,
            gaps_col = gaps_col,
            labels_col = labels_col,
            fontsize = fontsize)
}


#' @title foldchange_heatmap
#' 
#' @description Complex heatmap function that shows a barplot on top, can show several covariates and
#' can display 3 panels of data - IHC, gene expression and gene sets - specifically tailored towards time series
#' data that clusters the samples of each subject and adds gaps between subjects.
#' 
#' @param se SummarizedExperiment cube that contains all the data
#' @param subject_ids array indicating the subject ids for each sample, there should be multiple samples per id 
#' @param time_points timepoints for each sample
#' @param timepoint_order ordering of time points within each sample indicated as array or strings
#' @param labels_col Custom sample labels. By default these are the column labels of the passed SE
#' @param genes Genes to be displayed, which need to be contained in the SE cube
#' @param signatures List of gene vectors indicating different signatures / pathways / gene sets
#' @param ihc Array of names in the colData data.frame of the SE that include IHC
#' @param cluster_rows logic flag indicating whether the rows of the three data panels should be clustered. This is done using hierarchical clustering by Ward's method and using the optimal leaf ordering algorithm as implemented in the cba package.
#' @param covs Array of names in the colData data.frame of the SE indicating additional covariates that should be displayed on top of the three panels
#' @param annotation_colors (default: NA) List of colors for the covariates that are shown on top. 
#' @param column_order Indicating the covariate (name of the colData column) that the subjects should be ordered by. Notice that only the first timepoint is used for that, which in case of best pecent change should be consistent.
#' @param barplot_color This can be either NA - which colors all bars the color specified in barplot_legend, or a factorial where the levels need to correspond to the names of the barplot legend.
#' @param barplot_legend This can be a single color indicating the colors of all bars or a color scheme, which is an array of colors with names that correspond to the levels indicated in barplot_color
#' @param pathway_projection Pathway projection method. 'means' just averages the gene expression values of the passed gene sets, 'gsva' uses the R/Bioconductor package GSVA (using a Gaussian kernel in log2 space)
#' @param rowscale Logic flag indicating whether each row should be row scaled
#' @param rowscale_cap Caps the values in all three panels to +/- the cap value (default: 3), which avoids that the color schemes are driven by outliers. A values of 3 indicates that only values within 3 standard deviations are shown.
#' @param barplot_decreasing_order Logic flag that indicates the sorting order of the barplot. Can be NA - no sorting
#' @param barplot_censor Censor vector if the barplot is used for survival data
#' @param specific_sample_order Instead of sorting a specific ordering can be passed in case there are multiple indications or datasets that should be grouped 
#' @param fontsize Used fontsize in the heatmap
#' @param pathways_name_cap Pathway labels are be capped at 40  
#' @param main Title for the complex heatmap
#' 
#' @export
#' 
#' @importFrom grDevices gray.colors
#' @importFrom stats as.dendrogram
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @importFrom stats order.dendrogram
#' 
foldchange_heatmap <- function(se, 
                               subject_ids,
                               time_points,
                               timepoint_order = NULL,
                               labels_col = NULL,
                               genes, 
                               signatures = NULL, 
                               ihc = NULL, 
                               cluster_rows = TRUE,
                               covs = NULL, 
                               annotation_colors = NA,
                               column_order, 
                               barplot_color = NA,
                               barplot_legend = 'grey',
                               pathway_projection = 'means',
                               rowscale = TRUE,
                               rowscale_cap = 3,
                               barplot_decreasing_order = TRUE,
                               barplot_censor = NA,
                               specific_sample_order = NA,
                               fontsize = 8,
                               pathways_name_cap = 40, 
                               main){
  
  #select only the single genes of interest
  mat <- assays(se[genes,])$final
  if (cluster_rows && nrow(mat) > 2){
    mat <- mat[order.dendrogram(as.dendrogram(hcopt(as.dist(1-cor(t(mat))), method = "ward.D"))),]
  }
  
  if (is.null(labels_col)){
    labels_col <- paste(se[[subject_ids]], se[[time_points]])
  }
  
  #ihc values      
  if (!is.null(ihc)){
    ihc <- colData(se)[, ihc]
    ihc <- t(as.matrix(ihc))
    mat <- rbind(log2(ihc + 1), mat)
  }
  
  #row scale
  if (rowscale){
    mat <- rowScale(mat)
    #cap the rowscale values
    if (!is.na(rowscale_cap)){
      mat[mat > rowscale_cap] <- rowscale_cap
      mat[mat < -rowscale_cap] <- -rowscale_cap
    }
  }
  
  #pathway projection
  if (!is.null(signatures)){
    if (!all(sapply(signatures,function(x,y)all(x%in%rownames(y)),assays(se)$final))){
      warning("Not all genes in the signatures are represented within the RNAseq data. Filtering them now.")
      signatures <- lapply(signatures, function(x,y) x[x%in%y],rownames(assays(se)$final))
    }
    
    #pathway projection based on GSVA
    if (pathway_projection == 'GSVA'){
      pathways <- gsva(assays(se)$final, signatures)
      pathways <- pathways * 3
      
      #pathway projection based on means of log2 values
    } else if (pathway_projection == 'means'){
      pathways <- t(sapply(signatures,
                           function(gs,exp)
                             colMeans(exp[gs,]),
                           assays(se)$final))
      #rowscaling
      if (rowscale){
        pathways <- rowScale(pathways)
        if (!is.na(rowscale_cap)){
          pathways[pathways > rowscale_cap] <- rowscale_cap
          pathways[pathways < -rowscale_cap] <- -rowscale_cap
        }
      }         
      
    } else {
      stop('Pathway projection has to be "means" or "GSVA"')
    }      
    
    #if rows should be clustered
    if (cluster_rows && nrow(pathways) > 2){
      pathways <- pathways[order.dendrogram(as.dendrogram(hcopt(as.dist(1-cor(t(pathways))), method = "ward.D"))),]
    }
    
    if (!is.na(pathways_name_cap )){
      rownames(pathways) <- strtrim(rownames(pathways), pathways_name_cap)
    }
    
    #combine the sets
    mat <- rbind(mat, pathways)
    
  }
  
  if (is.null(timepoint_order)){
    timepoint_order <- unique(se[[time_points]])
  }
  
  #drop all samples that have NA values in the column order covariate
  filter <- !is.na(se[[column_order]]) & !grepl("NA", colnames(se)) & !is.na(colnames(se))
  mat <- mat[,filter]
  se <- se[,filter]
  labels_col <- labels_col[filter]
  
  
  cluster_cols <- as.numeric(se[[column_order]])
  
  
  #order the samples based on the column order using the first timepoint_order sample
  if (all(is.na(specific_sample_order))){
    temp_ord <- cluster_cols
    names(temp_ord) <- se[[subject_ids]]
    temp_ord <- temp_ord[se[[time_points]] == timepoint_order[1]]
    temp_ord <- sort(temp_ord, decreasing = barplot_decreasing_order)
    samp_order <- names(temp_ord)
    
  } else {
    stopifnot(all(specific_sample_order %in% se[[subject_ids]]))
    stopifnot(all(se[[subject_ids]] %in% specific_sample_order))
    samp_order <- specific_sample_order
  }
  
  #figure out the sample order
  sample_order <- lapply(samp_order, 
                         function(x, SUBJ, TP, TO){
                           current_subj <- (1:length(SUBJ))[SUBJ == x]
                           ret <- current_subj[match(TO, TP[current_subj])]
                           return(ret)       
                         },
                         se[[subject_ids]],
                         se[[time_points]],
                         timepoint_order)
  sample_order <- unlist(sample_order)
  ## Remove NA samples
  sample_order <- sample_order[!is.na(sample_order)]
  
  #order all data items
  mat <- mat[, sample_order]
  se <- se[, sample_order]
  labels_col <- labels_col[sample_order]
  cluster_cols <- cluster_cols[sample_order]
  
  if (!all(is.na(barplot_censor))){
    barplot_censor <- barplot_censor[sample_order]
  }
  
  #gap between the samples
  gaps_col <- c()
  pre_subj <- se[[subject_ids]][1]
  ## Iterate to identify gap position when detect a new subject
  for(i in 2:length(labels_col)) {
    subj <- se[[subject_ids]][i]
    if(subj != pre_subj) {
      gaps_col <- c(gaps_col, i - 1)
      pre_subj <- subj
    }
  }
  
  #gaps_col <- grep(timepoint_order[length(timepoint_order)], se[[time_points]])
  #gaps_col <- gaps_col[-length(gaps_col)]
  
  #gaps and colors
  color <- list()
  gaps <- c()
  offset <- 0
  
  #IHC
  if (!is.null(ihc)){
    offset <- nrow(ihc)
    gaps <- c(gaps, offset)
    color$IHC <- list(index = 1:nrow(ihc),
                      color = gray.colors(11, start = 0.2, end = 1)[11:1])
  }
  
  #genes
  color$`Gene Expression` <- list(index = (offset + 1) : (offset + length(genes)),
                                  color = brewer.pal(11, "RdBu")[11:1])
  offset <- offset + length(genes)
  
  #pathways
  if (!is.null(signatures)){
    gaps <- c(gaps, offset)
    color$Signatures <- list(index = (offset + 1) : nrow(mat),
                             color = brewer.pal(11, "RdYlBu")[11:1])
  }
  
  #covariates   
  if (is.null(covs)){
    col_annot <- NA
  } else if (class(covs) == 'list'){
    col_annot <- covs
  } else if (length(covs) == 1){
    col_annot <- data.frame(colData(se)[,covs])
    rownames(col_annot) <- colnames(mat)
  } else {
    col_annot <- data.frame(colData(se)[,covs])
  }
  
  #set the colors for the barplots
  if (is.na(barplot_color[1])){
    barplot_color <- rep(barplot_legend, ncol(mat))
    barplot_legend <- NA
    
  } else{
    barplot_color <- barplot_legend[match(se[[barplot_color]],names(barplot_legend))]
  }
  
  #heatmap
  g_heatmap(mat,
            na_col = 'white',
            main = main,
            scale = 'none',
            color = color,
            annotation_col = col_annot,
            annotation_colors = annotation_colors,
            cluster_rows = FALSE,
            cluster_cols = cluster_cols,
            barplot_label = column_order,
            barplot_color = barplot_color,
            barplot_censor = barplot_censor,
            barplot_legend = barplot_legend, 
            barplot_height = 80,
            na_as_X = TRUE,
            gaps_row = gaps,
            gaps_col = gaps_col,
            labels_col = labels_col,
            fontsize = fontsize)
}