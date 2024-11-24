#' Transform a scalar-on-shape partial linear single-index model
#'
#' @param Ydesign a n x L0 x m matrix of mediators.
#' @param Xdesign a n × p normalized design matrix.
#' @param efitBetas a p × L0 × m matrix of estimated functional coefficients.
#' @param npcval prespecified value for the number of principal components.
#' @param trunc truncation proportion of each side of g function.
#' @param Y a n x 1 matrix of outcomes.
#' @return a L0 × npc matrix of estimated eigenfunctions, a tn x 1 matrix of index, a tn x p matrix of linear covariates, a tn x npc matrix of nonlinear covariates, and a tn x 1 matrix of response variable.
#' @references https://cran.r-project.org/web/packages/refund/index.html.
#' @import refund
#' @export

SMA_PLSIMread = function(Ydesign, Xdesign, efitBetas, npcval, trunc, Y) {

  n <- dim(Ydesign)[1]
  L0 <- dim(Ydesign)[2]
  m <- dim(Ydesign)[3]
  p <- dim(Xdesign)[2]

  eg1 = refund::fpca.sc(Y=Ydesign[,,1], center=TRUE, npc=npcval)
  basis = eg1$efunctions
  z1 = eg1$scores
  alpha0 = matrix(1/sqrt(npcval), npcval,1)
  uu1 = z1%*%alpha0

  uu1X1z1 = cbind(uu1, Xdesign, z1, Y)
  quan <- quantile(uu1, c(trunc, 1-trunc))
  uuXz <- uu1X1z1[uu1X1z1[,1]>quan[1] & uu1X1z1[,1]<quan[2], ]
  tn <- dim(uuXz)[1]
  uu = matrix(uuXz[,1], tn,1)
  xmat = uuXz[,2:(1+p)]
  zmat = uuXz[,(2+p):(1+p+npcval)]
  ymat = matrix(uuXz[,ncol(uuXz)], tn,1)

  colnames(xmat) <- NULL
  colnames(zmat) <- NULL
  colnames(ymat) <- NULL

  return(list(basis=basis, uu=uu, xmat=xmat, zmat=zmat, ymat=ymat))
}
