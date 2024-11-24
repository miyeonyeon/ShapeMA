#' Prepare multivariate varying coefficient model formats to read a dataset
#'
#' @param data a matrix including id, W (clinical confounders), X (genetic exposure), Psi (shape mediators), Y (neurocognitive outcome), and Phi (FPC scores).
#' @param m the number of features.
#' @return a L0 x 3 matrix of tractdata, a n x p matrix of designdata, a m list with L0 x n matrix of diffusionFiles, and the number of features.
#' @references Hongtu Zhu. Runze Li. Linglong Kong. "Multivariate varying coefficient model for functional responses." Ann. Statist. 40 (5) 2634 - 2666, October 2012.
#' @export

SMA_MVCMprep = function(data, m=1) {
  # n is the sample size.
  # p is the number of covariates X and W.
  # L0 is the number of locations.

  # tractdata
  if (L0 %% 2 == 0){
    t1 <- seq(-L0/2, L0/2-1)
  } else {
    t1 <- seq(-(L0-1)/2, (L0-1)/2)
  }
  t2 <- matrix(0, L0, 2)
  tractdata <- cbind(t1, t2)

  # diffusionFiles
  diffusionFiles <- vector("list", m)
  diffusionFiles[[1]] <- aperm(as.matrix(data[, grep("Psi",colnames(data))]), c(2, 1))

  nofeatures <- length(diffusionFiles)
  designdata <- as.matrix(cbind(data[, c("X")], data[, grep("W",colnames(data))]))

  return(list(tractdata=tractdata, designdata=designdata, diffusionFiles=diffusionFiles, nofeatures=nofeatures))
}
