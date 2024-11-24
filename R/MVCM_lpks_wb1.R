#' MVCM_lpks_wb1 is to implement Zhu's (2010) method of local polynomial kernel smoothing (order = 1) with preselected bandwidth in MVCM
#'
#' Input:
#'     NoSetup      - col vector of [n, L0, p, m]
#'                    n = sample size
#'                    L0 = number of location of fiber tract
#'                    p = number of covariates
#'                    m = number of features
#'     arclength    - col vector of the arclength from one end to the other end
#'     Xdesign      - a n x p normalized design matrix.
#'     Ydesign      - a n x L0 x m matrix of fiber bundle diffusion properties
#'     mh           - a 1 x m vector of optimal bandwidth
#'     kstr         - Kernel function
#' Output:
#'     efitBetas    - a p x L0 x m matrix of estimated coefficients
#'     efitBetas1   - a 2*p x L0 x m matrix of estimated coefficients
#'     InvSigmats   - a p*2 x p*2 x L0 x m matrix
#'     efitYdesign  - a n x L0 x m matrix of estimated curves
#' ###################################################################################################################
#' Please run
#'    function [NoSetup, arclength, Xdesign,Ydesign]=MVCM_read(tractdata, designdata,  diffusionFiles, nofeatures, featuresname)
#' before you use MVCM_lpks_wb1
#' ###################################################################################################################
#' March 27, 2010 @ AA
#'
#' @references Hongtu Zhu. Runze Li. Linglong Kong. "Multivariate varying coefficient model for functional responses." Ann. Statist. 40 (5) 2634 - 2666, October 2012.
#' @import pracma
#' @export

MVCM_lpks_wb1 <- function(NoSetup, arclength, Xdesign, Ydesign, mh, kstr='exp(-.5*t^2)') {

  n <- NoSetup[1]
  L0 <- NoSetup[2]
  p <- NoSetup[3]
  m <- NoSetup[4]
  Tmat0 <- arclength %*% matrix(1, 1, L0) - matrix(1, L0, 1) %*% t(arclength)

  InvSigmats <- array(0, c(2*p, 2*p, L0, m))
  efitBetas <- array(0, c(p, L0, m))
  efitYdesign <- array(0, c(n, L0, m))
  efitBetas1 <- array(0, c(2*p, L0, m))

  for (mii in 1:m){
    h <- mh[mii]
    Tmat <- Tmat0 / h
    t <- Tmat
    Kmat <- eval(parse(text = kstr)) / h

    Sigmat <- array(0, c(2*p, 2*p, L0))
    Tempmat0 <- array(0, c(2*p, 2*p, L0))
    Tempmat <- array(0, c(n, L0, 2*p, L0))
    for (nii in 1:n){
      for (L0ii in 1:L0){
        TempA <- pracma::repmat(Xdesign[nii, ], 1, L0)
        TempB <- kronecker(t(Tmat[, L0ii, drop = FALSE]), matrix(1, 1, p))
        TempC <- rbind(TempA, TempA * TempB)
        TempD <- matrix(TempC, 2*p, L0)
        TempDD <- aperm(TempD, c(2, 1))
        TempEE <- TempDD * pracma::repmat(Kmat[, L0ii, drop = FALSE], 1, 2*p)
        Tempmat[nii, , , L0ii] <- TempEE
        Tempmat0[, , L0ii] <- t(TempEE) %*% TempDD
      }
      Sigmat <- Sigmat + Tempmat0
    }

    for (L0ii in 1:L0){
      InvSigmats[, , L0ii, mii] <- solve(Sigmat[, , L0ii])
    }

    AA <- array(replicate(m, Ydesign[, , mii]), c(n, L0, 2*p))
    for (L0ii in 1:L0){
      BB <- Tempmat[, , , L0ii]
      CC <- matrix(AA, nrow=dim(AA)[1]) * matrix(BB, nrow=dim(BB)[1])
      CCC <- matrix(apply(CC, 2, sum), nrow=1)
      KXYZmat <- apply(matrix(CCC, L0, 2*p), 2, sum)
      efitBetas1[, L0ii, mii] <- InvSigmats[, , L0ii, mii] %*% KXYZmat
    }

    efitBeta <- kronecker(diag(p), matrix(c(1, 0), nrow = 1)) %*% efitBetas1[, , mii]
    efitBetas[, , mii] <- efitBeta
    efitYdesign[, , mii] <- Xdesign %*% efitBeta
  }

  return(list(efitBetas=efitBetas, efitBetas1=efitBetas1, InvSigmats=InvSigmats, efitYdesign=efitYdesign))
}
