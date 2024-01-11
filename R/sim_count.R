#' Simulating single cell counts using estimated coefficients.
#'
#' This function performs an initial simulation of single cell counts using fitted coefficients from the previous step.
#'
#' @param coe.mat output of glm fitted coefficients from GLMsim_fit function.
#' @param cell_info output data frame from data_prepare function.
#' @param seed  used for set.seed() to generate random numbers. Default: 2022.

#' @return An initial simulated single cell count matrix, with rows corresponding to genes and columns corresponding to cells.
#' @export

GLMsim_simulation <- function(coe.mat, cell_info, seed = 2022){
  # Convert data frame for cell info df to numeric matrix
  cell_info_dat <- apply(cell_info, 2, as.character)
  cell_info_dat <- apply(cell_info_dat, 2, as.numeric)
  rownames(cell_info_dat) <- rownames(cell_info)
  cell_info_dat <- cbind(1, cell_info_dat)
  colnames(cell_info_dat)[1] <- "inter"

  # Calculate mu for entry of the simulated count
  sim.mu.mat <- coe.mat[,-ncol(coe.mat)] %*% t(cell_info_dat)
  sim.theta <- as.vector(coe.mat[,ncol(coe.mat)])

  # Simulate count
  sim.glm.count <- NULL
  set.seed(seed)
  theta.pois.genes <- which(is.na(sim.theta))
  theta.nb.genes <- which(!is.na(sim.theta))
  sim.pois.count <- NULL
  if(length(theta.pois.genes) > 0){
    sim.pois.count <- matrix(rpois(length(theta.pois.genes)*ncol(sim.mu.mat), lambda = exp(sim.mu.mat[theta.pois.genes,])), nrow = length(theta.pois.genes))
    rownames(sim.pois.count) <- rownames(coe.mat)[theta.pois.genes]
    colnames(sim.pois.count) <- rownames(cell_info)
  }
  sim.nb.count <- NULL
  if(length(theta.nb.genes) > 0){
    sim.nb.count <- matrix(rnbinom(length(theta.nb.genes)*ncol(sim.mu.mat), mu = exp(sim.mu.mat[theta.nb.genes,]), size = sim.theta[theta.nb.genes]), nrow = length(theta.nb.genes))
    rownames(sim.nb.count) <- rownames(coe.mat)[theta.nb.genes]
    colnames(sim.nb.count) <- rownames(cell_info)
  }
  sim.glm.count <- rbind(sim.nb.count, sim.pois.count)
  sim.glm.count <- sim.glm.count[rownames(coe.mat),]
  return(sim.glm.count)
}


