#' Get profile least squares estimator from selected bandwidth in a scalar-on-shape partial linear single-index model
#'
#' @param ymat a tn x 1 matrix of response variable.
#' @param xmat a tn x p matrix of normalized linear covariates.
#' @param zmat a tn x npc matrix of nonlinear covariates from centered mediators.
#' @param basis a L0 × npc matrix of estimated eigenfunctions.
#' @return a p x 1 matrix and a npc x 1 matrix of estimated coefficients, a tn x 1 matrix of estimated index and g function, and a L0 x 1 matrix of estimated beta function.
#' @references Liang H, Liu X, Li R, Tsai CL. Estimation and testing for partially linear single-index models. Annals of Statistics 2010; 38(6): 3811–3836.
#' @import PLSiMCpp
#' @export

SMA_PLSIMest = function(ymat, xmat, zmat, basis) {

  p <- dim(xmat)[2]
  tn <- dim(zmat)[1]
  npcval <- dim(zmat)[2]

  zeta_i = PLSiMCpp::plsim.ini(ymat~xmat|zmat)
  res_plsimest_cross = PLSiMCpp::plsim.bw(ymat~xmat|zmat, bandwidthList=c(0.02,0.04,0.06,0.08,0.10))
  fit = PLSiMCpp::plsim.est(ymat~xmat|zmat, zetaini=zeta_i, h=res_plsimest_cross$bandwidthBest)

  hat.v = matrix(fit$zeta[(npcval+1):(npcval+p)], p,1)
  hat.b = matrix(fit$zeta[1:npcval], npcval,1)
  hat.uu = zmat%*%hat.b
  EstG = matrix(fit$eta, tn,1)
  EstBeta = basis%*%hat.b

  return(list(hat.v=hat.v, hat.b=hat.b, hat.uu=hat.uu, EstG=EstG, EstBeta=EstBeta))
}
