
##TODO: add Wilcoxon

#' @title run_limma
#' 
#' @description Function that streamlines differntial expression using limma
#'
#' @param se SummarizedExperiment
#' @param group group to do the differential expression with 
#' @param main title of the heatmap
#' @param col_annot Heatmap annotation. Data frame that specifies the annotations shown on top of the heatmap. Each row defines the features for a specific column The rows in the data and in the annotation are matched using corresponding column names. Note that color schemes takes into account if variable is continuous or discrete.
#' @param annot_colors Heatmap annotation colors. List for specifying annotation_row and annotation_col track colors manually. It is possible to define the colors for only some of the features.
#' @param num_hm_genes Number of genes to be used in the heatmap
#' @param design_mat custom design matrix if a more complicated model is needed
#' @param coef coefficient to use for the likelihood ratio test when specifying a custom design matrix
#' @param cluster_cols (Default: T) boolean value determining whether columns should be clustered or hclust object.
#' @param cluster_rows (Default: T) boolean value determining whether rows should be clustered or hclust object.
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param assay_name_counts (Default: "raw_counts") the assay used for gene expression counts (i.e. \code{assays(se)[[assay_name_counts]]})
#' @param assay_name_final (Default: "final") the assay used for final plotting (i.e. \code{assays(se)[[assay_name_final]]})
#' 
#' @return table of differentially expressed genes
#' 
#' @export
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom stats model.matrix
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#'  
run_limma <- function(se, 
                      group=NULL, 
                      main = '', 
                      col_annot = NA, 
                      annot_colors = NA, 
                      num_hm_genes = 30,
                      design_mat = NULL,
                      assay_name_counts = "log2_tmm_cpm",
                      assay_name_final = "final",
                      coef = NULL,
                      cluster_cols = T,
                      cluster_rows = T,
                      show_colnames = T,
                      show_rownames = T){
  
  #gene expression counts
  mat <- assays(se)[[assay_name_counts]]
  
  #if there was no design matrix provided we generate one
  if (is.null(design_mat)){
    if (is.null(group)){
      stop('Either a design matrix or a binary group assignment vector needs to be specified to run the differential analysis.')
    }
    #in the design matrix we implicitly model the liver score if there is one
    if ('Liver.scores' %in% colnames(colData(se))){
      print('Correcting for liver score found in the SummarizedExperiment.')
      design_mat <- data.frame(liver=se$Liver.scores, group=group)
      design_mat <- model.matrix(~ liver + group, data = design_mat)
      
    } else {
      print('No liver score found in the SummarizedExperiment - not correcting.')
      design_mat <- model.matrix(~ group)
    }   
    
    #coefficient to use
    coef <- ncol(design_mat)
    
  } else if(is.null(coef)){
    coef <- ncol(design_mat)
  }
  
  fit <- lmFit(mat, design_mat)
  fit <- eBayes(fit)
  diff_genes <- topTable(fit, number = nrow(mat),  coef=coef)
  
  #Gene description
  if ('Gene_description' %in% colnames(rowData(se))){
    diff_genes$Gene_description <- rowData(se[rownames(diff_genes),])$description
  }
  
  #plot a heatmap of the top 30 genes
  plot_hm(mat = assays(se)[[assay_name_final]], 
          genes = rownames(diff_genes)[1:num_hm_genes], 
          main, 
          col_annot, 
          annot_colors,
          cluster_cols = cluster_cols,
          cluster_rows = cluster_rows,
          show_colnames = show_colnames,
          show_rownames = show_rownames)
  
  return(diff_genes)
}



#' @title run_diff_pathways
#' 
#' @description Function that streamlines differntial expression of pathways using limma
#'
#' @param pathways matrix of gene sets x samples
#' @param se SummarizedExperiment
#' @param group group to do the differential expression with 
#' @param main title of the heatmap
#' @param design_mat custom design matrix if a more complicated model is needed
#' @param coef coeffient to use for the likelihood ratio test when specifying a custom design matrix
#' @param plot_hm Should a heatmap be plotted? (default: T)
#' @param col_annot Heatmap annotation. Data frame that specifies the annotations shown on top of the heatmap. Each row defines the features for a specific column The rows in the data and in the annotation are matched using corresponding column names. Note that color schemes takes into account if variable is continuous or discrete.
#' @param annot_colors Heatmap annotation colors. List for specifying annotation_row and annotation_col track colors manually. It is possible to define the colors for only some of the features.
#' @param num_hm_genes Number of genes to be used in the heatmap
#' @param cluster_cols (Default: T) boolean value determining whether columns should be clustered or hclust object.
#' @param cluster_rows (Default: T) boolean value determining whether rows should be clustered or hclust object.
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' 
#' @return data frame of differentially expressed pathways
#' 
#' @export
#' 
#' @importFrom limma topTable
#' @importFrom limma eBayes
#' @importFrom limma lmFit
#' @importFrom stats model.matrix
#' @importFrom stats p.adjust
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' 
run_diff_pathways <- function(pathways, 
                              se, 
                              group = NULL, 
                              main = NULL,
                              design_mat = NULL,
                              coef = NULL,
                              plot_hm = T,
                              col_annot = NA, 
                              annot_colors = NA,
                              num_hm_genes = 30,
                              show_colnames = T,
                              show_rownames = T,
                              cluster_cols = T,
                              cluster_rows = T){
  
  #if a design matrix was provided   
  if (is.null(design_mat)){
    if (is.null(group)){
      stop('You need to specify either a design matrix or a group vector.')
    }
    design_mat <- model.matrix(~ group)
  } 
  
  #make sure we are not trying to show more genes / pathways than we have 
  if (num_hm_genes > nrow(pathways)){
    num_hm_genes <- nrow(pathways)
  }
  
  #coefficient to use
  if(is.null(coef)){
    coef <- ncol(design_mat)
  }
  
  #fit a linear model using limma
  fit <- lmFit(pathways, design_mat)
  fit <- eBayes(fit)
  
  #correcting for multiple testing using BH - FDR
  diff_pathways <- data.frame(topTable(fit, coef, nrow(pathways)))
  
  if (plot_hm){
    #plot a heatmap of the top 30 genes
    plot_hm(pathways, 
            rownames(diff_pathways)[1:num_hm_genes], 
            main, 
            col_annot, 
            annot_colors,
            cluster_cols = cluster_cols,
            cluster_rows = cluster_rows,
            show_colnames = show_colnames,
            show_rownames = show_rownames)
  }
  
  return(diff_pathways)
}


#' @title extract_correlations
#' 
#' @description Function that streamlines differntial expression of pathways using limma
#'
#' @param mat matrix of genes x samples
#' @param num_cov numeric covariate to correlate against
#' @param method used for correlation ('pearson' or 'spearman')
#' @param filer_genes Should lowly expressed genes (rowMeans(mat) < 1) be filtered? (default: F)
#' 
#' @return table with correlations
#' 
#' @export
#' 
#' @importFrom stats p.adjust
#' @importFrom stats sd
#' @importFrom stats cor.test
#'  
extract_correlations <- function(mat, num_cov, method = 'pearson', filer_genes = F){
  
  stopifnot(method %in% c('pearson', 'spearman'))
  
  #drop genes that have no variance
  mat <- mat[apply(mat,1,sd) > 0, ]
  
  if (filer_genes){
    # drop genes that have less than 1 log2 CPM avg expression
    mat <- mat[rowMeans(mat) > 1, ]
  }
  
  #Calculate correlations beased on Spearman
  if (method == 'pearson'){
    all_cor <- apply(mat, 1, cor.test, num_cov)
    all_cor <- data.frame(t(sapply(all_cor,
                                   function(x)
                                     c(x$statistic, 
                                       x$estimate, 
                                       x$conf.int, 
                                       x$p.value))))
    names(all_cor) <- c('Statistic','cor','95.CI.lower','95.CI.upper','p.value')
    
    #Calculate correlations beased on Spearman
  } else {
    all_cor <- suppressWarnings(apply(mat, 1, cor.test, num_cov, method='spearman'))
    all_cor <- data.frame(t(sapply(all_cor,
                                   function(x)
                                     c(x$statistic, 
                                       x$estimate, 
                                       x$p.value))))
    names(all_cor) <- c('Statistic','cor','p.value')
  }
  
  all_cor$mean.log2.CPM <- rowMeans(mat)
  all_cor$FDR <- p.adjust(all_cor$p.value,method = 'BH')
  all_cor <- all_cor[order(all_cor$p.value,decreasing = F),]
  all_cor <- data.frame(all_cor, stringsAsFactors = F)
  
  return(all_cor)   
}


#' @title plot_volcano
#' 
#' @description Volcano plot to show top differential expressed genes
#'
#' @param diff_genes data.frame after differential gene analysis: with logFC, gene.name and p.val columns at least
#' @param p.val.cutoff cutoff of p.val to show the gene labels. 
#' @param p.val.name column name for the p value in diff_genes
#' @param p.val.lower.limit the lowest p.val limit, if p.val is smaller than it, it will be set to this value to avoid volcano becoming too high.
#' @param logfc.name column name for the log fold change in diff_genes
#' @param gene.name by default, it uses rownames as gene name, but if there is particular column for gene.name, specifity it here.
#' @param n.top (default: 20) number of top differential genes to show the gene labels.
#' @param n.max (default: 50) maximum number of genes to show on the volcano to avoid over-populating the volcano.
#' A gene will be labeled if it either meets the p.val.cutoff or meets the n.top criteria.
#' @param ylab label of y-axis
#' @param title title of the plot
#' 
#' @return a ggplot object.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' plot_volcano(diff_genes, p.val.cutoff = 1e-4, logfc.name = "log2FoldChange", p.val.name = "pvalue")
#' }
plot_volcano <- function(diff_genes, 
                         p.val.cutoff = 1e-3, 
                         p.val.name = "pvalue", 
                         p.val.lower.limit = 1e-8,
                         logfc.name = "log2FoldChange",
                         gene.name = "rownames",
                         n.top = 20,
                         n.max = 50, 
                         ylab = "-log10(p.val)",
                         title = ""){
  
  y.max <- -1 * log10(p.val.lower.limit)
  
  diff_genes <- diff_genes[order(diff_genes[,p.val.name], decreasing = T),]
  diff_genes$log.pval <- -log10(diff_genes[,p.val.name])
  diff_genes$log.pval[diff_genes$log.pval > y.max] <- y.max
  diff_genes$sig <- diff_genes[,p.val.name] < p.val.cutoff
  diff_genes$log2FC <- diff_genes[,logfc.name]
  
  if(gene.name == "rownames") {
    diff_genes$gene <- rownames(diff_genes)
  } else {
    diff_genes$gene <- diff_genes[,gene.name]
  }
  
  if(n.top > sum(diff_genes$sig, na.rm = T)) {
    highlight <- diff_genes[1:n.top,]
  } else {
    highlight <- subset(diff_genes, sig)
  }
  highlight <- highlight[!duplicated(highlight),]
  
  ggplot2::ggplot(diff_genes, aes(x = log2FC, y = log.pval)) + geom_point(aes(col=sig), alpha = 0.5) +
    scale_color_manual(values=c("black", "red")) + 
    ggrepel::geom_text_repel(data=highlight, aes(label=gene, col=sig), size = 3) +
    ylab(ylab) + xlab(logfc.name) + 
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position="none") + ggtitle(title)
  
}
