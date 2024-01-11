check_outlier <- function(sim_count, real_count, nmad){
  sim.mean <- log(rowMeans(sim_count+1))
  real.mean <- log(rowMeans(real_count+1))
  common_genes <- intersect(rownames(sim_count), rownames(real_count))
  mean_dif <- abs(sim.mean[common_genes] - real.mean[common_genes])
  outlier_cutoff <- nmad * mad(mean_dif)
  outlier_genes <- common_genes[which(mean_dif > outlier_cutoff)]
  return(outlier_genes)
}

fit_robust_glm <- function(gene_umi, design_mat, cell.info, glm_form, start_coe){
  df <- cbind(gene_umi, cell.info)
  colnames(df)[1] <- "y"
  tryCatch({
    robnbm <- glmrob.nb(y = gene_umi, X = design_mat)
    new.glm.fit <- MASS::glm.nb(as.formula(glm_form), data = df, start = robnbm$coef[2:length(robnbm$coef)], init.theta = 1/robnbm$coef[1])
    coe.list <- as.list(new.glm.fit$coefficients)
    names(coe.list) <- c("b0", colnames(df)[2:ncol(df)])
    coe.list[["theta"]] <- new.glm.fit$theta
    return(coe.list)
  },
  error = function(e){
    tryCatch({
      m <- MASS::glm.nb(as.formula(glm_form), data = df, start = start_coe[1:(length(start_coe)-1)], init.theta = start_coe[length(start_coe)])
      coe.list <- as.list(m$coefficients)
      names(coe.list) <- c("b0", colnames(df)[2:ncol(df)])
      coe.list[["theta"]] <- m$theta
      return(coe.list)
    },
    error = function(e){
      tryCatch({
        m <- MASS::glm.nb(as.formula(glm_form), data = df, start = start_coe[1:(length(start_coe)-1)], init.theta = start_coe[length(start_coe)], maxit = 1000)
        coe.list <- as.list(m$coefficients)
        names(coe.list) <- c("b0", colnames(df)[2:ncol(df)])
        coe.list[["theta"]] <- m$theta
        return(coe.list)
      },
      error = function(e){
        m <- glm(as.formula(glm_form), data = df, family = poisson())
        coe.list <- as.list(m$coefficients)
        names(coe.list) <- c("b0", colnames(df)[2:ncol(df)])
        coe.list[["theta"]] <- NA
        return(coe.list)
      })
    })
  })
}

#' Simulating counts with outlier genes.
#'
#' To mitigate the effects of extreme counts from outlier genes, this function will check if outlier genes exist at first. If outlier genes exist in the initial simulated data, this function allows users to choose one of three methods: robust negative binomial GLM, Winsorizing and trended coefficient to refit those outlier genes. Then simulate counts for those outlier genes using refitted coefficients.
#'
#' @param real_count real data count matrix. Rows: genes. Columns: cells.
#' @param cell_info output data frame from data_prepare function.
#' @param coe.mat output of glm fitted coefficients from GLMsim_fit function.
#' @param sim.glm.count the initial simulated count matrix from GLMsim_simulation function.
#' @param fit.method the method used to fit outlier genes. Users should choose one of the methods from “robust_glm” or “clipping” or “trend”. Default: “trend” method.
#' @param nmad number of median absolute deviation (MAD) to determine outlier genes. Default: 30.
#' @param parallel whether to fit outlier genes into GLM in parallel or not. Default: FALSE.
#' @param trunc_q a vector consisting of 2 elements to indicate the quantile range for the cut-off of the Winsorizing method. The default is c(0.05, 0.95).

#' @return A simulated single cell count matrix after correcting outlier genes, with rows corresponding to genes and columns corresponding to cells. If no outlier genes exist, the output simulated count will be same as that from the GLMsim_simulation function.
#' @export

GLMsim_sim_outlier <- function(real_count, cell_info, coe.mat, sim.glm.count, fit.method = "trend", nmad = 30,  parallel = FALSE, trunc_q = c(0.05, 0.95)){
  real_count <- as.matrix(real_count)
  real_count <- real_count[which(rowSums(real_count > 0) >= 4),]
  if(is.null(rownames(real_count))){
    rownames(real_count) <- paste0("Gene", 1:nrow(real_count))
  }
  if(is.null(rownames(sim.glm.count))){
    rownames(sim.glm.count) <- paste0("Gene", 1:nrow(sim.glm.count))
  }
  outlier.genes <- check_outlier(sim.glm.count, real_count, nmad)

  if(length(outlier.genes) > 0){
    if(fit.method == "robust_glm"){
      cell_design_mat <- apply(cell_info, 2, as.character)
      cell_design_mat <- apply(cell_design_mat, 2, as.numeric)
      rownames(cell_design_mat) <- rownames(cell_info)
      glm.form <- paste("y", paste(colnames(cell_info), collapse=" + "), sep=" ~ ")
      start_coe_mat <- coe.mat[outlier.genes,]
      if(length(outlier.genes) > 1){
        if(parallel){
          cl <- makeCluster(detectCores())
          clusterExport(cl, varlist = c("real_count", "cell_design_mat", "cell_info", "glm.form", "outlier.genes", "start_coe_mat", "fit_robust_glm"), envir=environment())
          clusterEvalQ(cl, library(MASS))
          outlier.coe.list <- parLapply(cl, 1:length(outlier.genes), fun = function(x){
            outlier.gene.fit <- fit_robust_glm(real_count[outlier.genes[x],], cell_design_mat, cell_info, glm.form, start_coe_mat[outlier.genes[x],])
            return(outlier.gene.fit)
          })
          stopCluster(cl)
        }else{
          outlier.coe.list <- lapply(1:length(outlier.genes), FUN = function(x){
            fit_robust_glm(real_count[outlier.genes[x],], cell_design_mat, cell_info, glm.form, coe.mat[outlier.genes,])
          })
        }
      }else{
        outlier.coe.list <- list()
        outlier.coe.list[[outlier.genes]] <- fit_robust_glm(real_count[outlier.genes,], cell_design_mat, cell_info, glm.form, start_coe_mat[outlier.genes,])
      }
      coe.outlier.mat <- matrix(NA, nrow = length(outlier.genes), ncol = ncol(coe.mat))
      for(i in 1:length(outlier.coe.list[[1]])){
        coe.outlier.mat[,i] <- unlist(lapply(outlier.coe.list, FUN = function(x){x[[i]]}))
      }
      colnames(coe.outlier.mat) <- names(outlier.coe.list[[1]])
      rownames(coe.outlier.mat) <- outlier.genes

    }else if(fit.method == "clipping"){
      coe.mat2 <- coe.mat
      coe.mat2[,"theta"] <- log(coe.mat2[,"theta"])
      truncate.cutoff <- apply(coe.mat2, 2, FUN = function(x){
        # Default quantile: 0.05, 0.95
        quantile(x, trunc_q, na.rm = TRUE)
      })

      coe.outlier.mat <- matrix(coe.mat2[outlier.genes,], ncol = ncol(coe.mat2), nrow = length(outlier.genes))
      colnames(coe.outlier.mat) <- colnames(coe.mat2)
      rownames(coe.outlier.mat) <- outlier.genes
      for(i in 1:ncol(coe.outlier.mat)){
        coe.value <- coe.outlier.mat[,i]
        coe.value[which(coe.value < truncate.cutoff[1,i])] <- truncate.cutoff[1,i]
        coe.value[which(coe.value > truncate.cutoff[2,i])] <- truncate.cutoff[2,i]
        coe.outlier.mat[,i] <- coe.value
      }
      coe.outlier.mat[,"theta"] <- exp(coe.outlier.mat[,"theta"])

    }else if(fit.method == "trend"){
      coe.mat2 <- coe.mat
      coe.mat2[,"theta"] <- log(coe.mat2[,"theta"])
      mean_real <- log(rowMeans(real_count+1))
      common.genes <- intersect(names(mean_real), rownames(coe.mat2))
      mean_real <- mean_real[common.genes]
      coe.outlier.mat <- apply(coe.mat2[common.genes,], 2, FUN = function(x){
        data.df <- data.frame("x" = mean_real, "y" = x)
        loess_m <- loess(y ~ x, data = data.df)
        outlier.df <- data.frame("x" = data.df[outlier.genes,1])
        new.coe <- predict(loess_m, newdata = outlier.df)
        return(new.coe)
      })
      coe.outlier.mat <- matrix(coe.outlier.mat, ncol = ncol(coe.mat2), nrow = length(outlier.genes))
      colnames(coe.outlier.mat) <- colnames(coe.mat2)
      rownames(coe.outlier.mat) <- outlier.genes
      coe.outlier.mat[,"theta"] <- exp(coe.outlier.mat[,"theta"])
    }
    outlier.sim.count <- GLMsim_simulation(coe.outlier.mat, cell_info)
    if(!is.matrix(outlier.sim.count)){
      outlier.sim.count <- matrix(outlier.sim.count, ncol = ncol(sim.glm.count), nrow = nrow(coe.outlier.mat))
      colnames(outlier.sim.count) <- colnames(sim.glm.count)
      rownames(outlier.sim.count) <- rownames(coe.outlier.mat)
    }
    sim.glm.count[rownames(outlier.sim.count),] <- outlier.sim.count
  }
  return(sim.glm.count)
}


