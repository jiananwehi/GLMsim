#' Fitting each gene count into a generalized linear model.
#'
#' This function fits UMI count for each gene into a generalized linear model, which aims to obtain the estimated coefficients for biological effect and unwanted variation.
#'
#' @param real_count real data count matrix. Rows: genes. Columns: cells.
#' @param cell_info a data frame of cell information that outputs from data_prepare function.
#' @param parallel whether to fit all genes into GLM in parallel. Default: TRUE.
#' @param ncores number of cores that can be used for parallel computing. Default: "all", which means using all cores. If that does not work, try 6 or 8.

#' @return A matrix of fitted coefficients from the real experimental data. In this matrix, each row represents a gene, and each column corresponds to a coefficient.
#'
#' @importFrom MASS glm.nb
#' @importFrom MASS glm
#' @import parallel
#' @export

GLMsim_fit <- function(real_count, cell_info, parallel = TRUE, ncores = "all"){
  real_count <- as.matrix(real_count)
  real_count <- real_count[which(rowSums(real_count > 0) >= 4),]
  if(is.null(rownames(real_count))){
    rownames(real_count) <- paste0("Gene", 1:nrow(real_count))
  }
  if(is.null(colnames(real_count))){
    colnames(real_count) <- pate0("Cell", 1:ncol(real_count))
  }

  # Fit each gene into NB GLM
  glm.form <- paste("y", paste(colnames(cell_info), collapse=" + "), sep=" ~ ")

  if(parallel){
    # Parallel fit GLM
    if(ncores == "all"){
      cl <- makeCluster(detectCores())
    }else{
      cl <- makeCluster(ncores)
    }
    clusterExport(cl, varlist = c("real_count", "cell_info", "glm.form"), envir=environment())
    clusterEvalQ(cl, library(MASS))
    gene.glm.coe <- parApply(cl, real_count, 1, FUN = function(gene.umi){
      df <- cbind(gene.umi, cell_info)
      colnames(df)[1] <- "y"
      tryCatch({
        m <- MASS::glm.nb(as.formula(glm.form), data = df)
        coe.list <- as.list(m$coefficients)
        names(coe.list) <- c("b0", colnames(df)[2:ncol(df)])
        coe.list[["theta"]] <- m$theta
        return(coe.list)
      },
      error = function(e){
        tryCatch({
          m <- MASS::glm.nb(as.formula(glm.form), data = df, maxit = 1000)
          coe.list <- as.list(m$coefficients)
          names(coe.list) <- c("b0", colnames(df)[2:ncol(df)])
          coe.list[["theta"]] <- m$theta
          return(coe.list)
        },
        error = function(e){
          m <- glm(as.formula(glm.form), data = df, family = poisson())
          coe.list <- as.list(m$coefficients)
          names(coe.list) <- c("b0", colnames(df)[2:ncol(df)])
          coe.list[["theta"]] <- NA
          return(coe.list)
        })
      })
    })
    stopCluster(cl)
  }else{
    gene.glm.coe <- apply(real_count, 1, FUN = function(gene.umi){
      df <- cbind(gene.umi, cell_info)
      colnames(df)[1] <- "y"
      tryCatch({
        m <- MASS::glm.nb(as.formula(glm.form), data = df)
        coe.list <- as.list(m$coefficients)
        names(coe.list) <- c("b0", colnames(df)[2:ncol(df)])
        coe.list[["theta"]] <- m$theta
        return(coe.list)
      },
      error = function(e){
        tryCatch({
          m <- MASS::glm.nb(as.formula(glm.form), data = df, maxit = 1000)
          coe.list <- as.list(m$coefficients)
          names(coe.list) <- c("b0", colnames(df)[2:ncol(df)])
          coe.list[["theta"]] <- m$theta
          return(coe.list)
        },
        error = function(e){
          m <- glm(as.formula(glm.form), data = df, family = poisson())
          coe.list <- as.list(m$coefficients)
          names(coe.list) <- c("b0", colnames(df)[2:ncol(df)])
          coe.list[["theta"]] <- NA
          return(coe.list)
        })
      })
    })
  }

  # Organize the fitted GLM coefficient from list to matrix
  coe.mat <- matrix(NA, nrow = length(gene.glm.coe), ncol = length(gene.glm.coe[[1]]))
  for(i in 1:length(gene.glm.coe[[1]])){
    coe.mat[,i] <- unlist(lapply(gene.glm.coe, FUN = function(x){x[[i]]}))
  }
  colnames(coe.mat) <- names(gene.glm.coe[[1]])
  rownames(coe.mat) <- names(gene.glm.coe)

  return(coe.mat)
}











