#' Collecting cell information from the original experimental data
#'
#' This function collects information from the original experimental data. The library size, biological and batch information are extracted from the original experimental data, which will be prepared for fitting each gene into a generalized linear model for the next step.
#'
#' @param real_count real data count matrix. Rows: genes. Columns: cells.
#' @param biology_vec a factor vector of biology groups.
#' @param batch_vec a factor vector of batches.

#' @return A data frame to reorganize the information for biology and unwanted variations.
#' @export

GLMsim_cell_info <- function(real_count, biology_vec = NULL, batch_vec = NULL){
  real_count <- as.matrix(real_count)
  real_count <- real_count[which(rowSums(real_count > 0) >= 4),]
  if(is.null(rownames(real_count))){
    rownames(real_count) <- paste0("Gene", 1:nrow(real_count))
  }
  if(is.null(colnames(real_count))){
    colnames(real_count) <- pate0("Cell", 1:ncol(real_count))
  }
  lib.size <- log(colSums(real_count))
  cell.info <- data.frame("l" = lib.size)

  # Define biology data frame for glm fitting
  if(!is.null(biology_vec) || length(unique(biology_vec)) > 1){
    bio.group <- list()
    if(!is.factor(biology_vec)){
      biology_vec <- factor(biology_vec)
    }
    for(i in 2:length(levels(biology_vec))){
      bio.grp.vec <- rep(0, ncol(real_count))
      bio.grp.vec[which(biology_vec %in% levels(biology_vec)[i])] <- 1
      bio.grp.vec <- factor(bio.grp.vec, levels = c(0, 1))
      bio.group[[levels(biology_vec)[i]]] <- bio.grp.vec
    }
    for(i in 1:length(bio.group)){
      bio.df <- data.frame(bio.group[[i]])
      colnames(bio.df) <- paste0("x", i)
      cell.info <- cbind(cell.info, bio.df)
    }
    cell.info <- cell.info[,c(2:ncol(cell.info),1)]
  }

  # Define batch data frame for glm fitting
  if(!is.null(batch_vec) || length(unique(batch_vec)) > 1){
    batch.group <- list()
    if(!is.factor(batch.group)){
      batch_vec <- factor(batch_vec)
    }
    for(i in 2:length(levels(batch_vec))){
      batch.grp.vec <- rep(0, ncol(real_count))
      batch.grp.vec[which(batch_vec %in% levels(batch_vec)[i])] <- 1
      batch.grp.vec <- factor(batch.grp.vec, levels = c(0, 1))
      batch.group[[levels(batch_vec)[i]]] <- batch.grp.vec
    }
    for(i in 1:length(batch.group)){
      batch.df <- data.frame(batch.group[[i]])
      colnames(batch.df) <- paste0("w", i+1)
      cell.info <- cbind(cell.info, batch.df)
    }
  }

  return(cell.info)
}
