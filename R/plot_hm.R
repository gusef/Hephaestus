## simple heatmap wrapper


#' @title rowScale
#' 
#' @description Simple function to perform row scaling
#'
#' @param x gene expression matrix
#' @return row scaled gene expression matrix
#'
#' @export
#' 
rowScale <- function(x){
  rm <- rowMeans(x, na.rm = T)
  x_norm <- sweep(x, 1, rm)
  sx <- apply(x_norm, 1, sd, na.rm = T)
  x_norm <- sweep(x_norm, 1, sx, "/")
  return(x_norm)
}


#' @title plot_hm
#' 
#' @description Simple heatmap function that is used to output differential expression results
#'
#' @param mat gene expression matrix
#' @param genes genes to include
#' @param main title of the heatmap
#' @param col_annot data frame that specifies the annotations shown on top of the heatmap. Each row defines the features for a specific column The rows in the data and in the annotation are matched using corresponding column names. Note that color schemes takes into account if variable is continuous or discrete.
#' @param annot_colors list for specifying annotation_row and annotation_col track colors manually. It is possible to define the colors for only some of the features.
#' @param value_cap (Default: 3) Each expression value will be capped at +/- of this value after row scaling. This corresponds to a cap of 3 standard devivations.
#' @param string_cap (default: 40) Row name labels are capped at 40, which is particularly useful for pathways 
#' @param col (Default: Blue-Red) Color scheme
#' @param font.size (Default: 8) Size of the label fonts
#' @param cluster_cols (Default: T) boolean value determining whether columns should be clustered or hclust object.
#' @param cluster_rows (Default: T) boolean value determining whether rows should be clustered or hclust object.
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' 
#' @return row scaled gene expression matrix
#' @export
#' 
plot_hm <- function(mat, 
                    genes, 
                    main = '', 
                    col_annot = NA, 
                    annot_colors = NA, 
                    value_cap = 3, 
                    string_cap = 40, 
                    col = brewer.pal(11, "RdBu")[11:1],
                    font.size=8,
                    cluster_cols = T,
                    cluster_rows = T,
                    show_colnames = T,
                    show_rownames = T){
  #using only provided genes
  mat <- mat[genes,]
  
  #rowscale and cap the values
  mat <- rowScale(mat)
  mat[mat > value_cap] <- value_cap
  mat[mat < -value_cap] <- -value_cap
  
  #limit the label length
  rownames(mat) <- strtrim(rownames(mat), string_cap)
  
  if (!all(is.na(col_annot))){
    col_annot <- data.frame(col_annot,stringsAsFactors = F)
  }
  
  g_heatmap(mat,
            main = main,
            scale = 'none',
            annotation_col = col_annot,
            annotation_colors = annot_colors,
            cluster_cols = cluster_cols,
            cluster_rows = cluster_rows,
            color = col,
            fontsize = font.size,
            show_colnames = show_colnames,
            show_rownames = show_rownames)
}