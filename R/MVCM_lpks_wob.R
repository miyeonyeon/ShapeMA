#' MVCM_lpks_wob is to implement Zhu's (2010) method of local polynomial kernel smoothing (order = 1) without preselected bandwidth in MVCM
#'
#' Input:
#'     NoSetup    - col vector of [n, L0, p, m]
#'                  n = sample size
#'                  L0 = number of location of fiber tract
#'                  p = number of covariates
#'                  m = number of features
#'     arclength  - col vector of the arclength from one end to the other end
#'     Xdesign    - a n x p normalized design matrix.
#'     Ydesign    - a n x L0 x m matrix of fiber bundle diffusion properties
#'     kstr       - Kernel function
#' Output:
#'     mh         - a 1 x m vector of optimal bandwidth
#'     GCVs       - a nh x m matrix of GCVs
#'     vh         - a 1 x nh vector of bandwidth
#' ###################################################################################################################
#' Please run
#'    function [NoSetup, arclength, Xdesign,Ydesign]=MVCM_read(tractdata, designdata,  diffusionFiles, nofeatures, featuresname)
#' before you use MVCM_lpks_wob
#' ###################################################################################################################
#' March 27, 2010 @ AA
#'
#' @references Hongtu Zhu. Runze Li. Linglong Kong. "Multivariate varying coefficient model for functional responses." Ann. Statist. 40 (5) 2634 - 2666, October 2012.
#' @import pracma
#' @export

MVCM_lpks_wob <- function(NoSetup, arclength, Xdesign, Ydesign, kstr='exp(-.5*t^2)') {

  n <- NoSetup[1]
  L0 <- NoSetup[2]
  p <- NoSetup[3]
  m <- NoSetup[4]

  xrange <- base::max(arclength) - base::min(arclength)
  nh <- base::max(30, floor(L0 / 2))
  hmin <- 1.0 * xrange / L0
  hmax <- xrange / 8
  efitYdesign <- array(0, c(n, L0, m))
  GCVs <- matrix(0, nh, m)
  vh <- pracma::logspace(log10(hmin), log10(hmax), nh)
  Tmat0 <- arclength %*% matrix(1, 1, L0) - matrix(1, L0, 1) %*% t(arclength)

  for (nhii in 1:nh){
    h <- vh[nhii]
    Tmat <- Tmat0 / h
    t <- Tmat
    Kmat <- eval(parse(text = kstr)) / h

    Tempmat0 <- array(0, c(2*p, 2*p, L0, n))
    Tempmat <- array(0, c(n, L0, 2*p, L0))
    for (nii in 1:n){
      for (L0ii in 1:L0){
        TempA <- pracma::repmat(Xdesign[nii, ], 1, L0)
        TempB <- kronecker(t(Tmat[, L0ii]), matrix(1, 1, p))
        TempC <- rbind(TempA, TempA * TempB)
        TempD <- matrix(TempC, nrow = 2 * p, ncol = L0)
        TempDD <- t(TempD)
        TempEE <- TempDD * matrix(rep(Kmat[,L0ii], each = 2 * p), nrow = L0, byrow = TRUE)
        Tempmat[nii, , , L0ii] <- TempEE
        Tempmat0[, , L0ii, nii] <- t(TempEE) %*% TempDD
      }
    }

    # helper function to repeat copies of a 2D matrix in the 3rd and 4th dimensions
    repmat2to4D <- function(mat, dimcopy){
      array <- array(0, c(dim(mat)[1], dim(mat)[2], dimcopy[1], dimcopy[2]))
      for (i in 1:dimcopy[1]){
        for (j in 1:dimcopy[2]){
          array[, , i, j] <- mat
        }
      }
      return(array)
    }

    Tempmat1 <- array(0, c(n, 2*p, L0, m))
    for (mii in 1:m){
      AA <- repmat2to4D(Ydesign[, , mii], c(2*p, L0))
      Tempmat10 <- matrix(Tempmat, nrow=dim(Tempmat)[1]) * matrix(AA, nrow=dim(AA)[1])
      Tempmat11 <- array(Tempmat10, c(n, L0, 2*p, L0))
      Tempmat1[, , , mii] <- apply(Tempmat11, c(1,3,4), sum)
    }

    for (mii in 1:m){
      for (L0ii in 1:L0){
        Tempmat2 <- drop(rowSums(Tempmat0[, , L0ii, , drop = FALSE], dims = 3))
        Tempmat3 <- drop(colSums(Tempmat1[, , L0ii, mii, drop = FALSE]))
        for (nii in 1:n){
          Tempmat22 <- Tempmat2 - Tempmat0[, , L0ii, nii]
          Tempmat33 <- Tempmat3 - Tempmat1[nii, , L0ii, mii, drop = FALSE]
          AA <- solve(Tempmat22) %*% Tempmat33
          efitYdesign[nii, L0ii, mii] <- sum(Xdesign[nii, ] * kronecker(diag(p), matrix(c(1, 0), nrow = 1)) %*% AA)
        }
      }
      GCVs[nhii, mii] <- sum((Ydesign[, , mii] - efitYdesign[, , mii])^2) / (n*L0)
    }
  }

  flag <- apply(GCVs, 2, which.min)
  mh <- vh[flag]

  return(list(mh=mh, GCVs=GCVs, vh=vh))
}
