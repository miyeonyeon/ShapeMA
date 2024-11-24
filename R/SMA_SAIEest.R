#' Estimate the spatial average indirect effects
#'
#' @param data a matrix including id, W (clinical confounders), X (genetic exposure), Psi (shape mediators), Y (neurocognitive outcome), and Phi (FPC scores).
#' @param m the number of features.
#' @param EstC the estimated constant of the average indirect effects.
#' @param trunc truncation proportion of each side of g function.
#' @return a L0 x 1 matrix of the estimated spatial average indirect effects.
#' @import refund PLSiMCpp
#' @export

SMA_SAIEest = function(data, m=1, EstC, trunc) {

  # alpha(s)
  SMA1a <- SMA_MVCMprep(data, m=1)
    tractdata <- SMA1a$tractdata
    designdata <- SMA1a$designdata
    diffusionFiles <- SMA1a$diffusionFiles
    nofeatures <- SMA1a$nofeatures

  SMA1b <- MVCM_read(tractdata, designdata, diffusionFiles, nofeatures)
    NoSetup <- SMA1b$NoSetup
    arclength <- SMA1b$arclength
    Xdesign0 <- SMA1b$Xdesign
    Ydesign <- SMA1b$Ydesign
    Xdesign <- cbind(designdata[, 1], Xdesign0[, 2:p])
  SMA1c <- MVCM_lpks_wob(NoSetup, arclength, Xdesign, Ydesign)
    mh <- SMA1c$mh
  SMA1d <- MVCM_lpks_wb1(NoSetup, arclength, Xdesign, Ydesign, mh)
    efitBetas <- SMA1d$efitBetas

  # beta(s)
  npcval = length(grep("Phi",colnames(data)))
  SMA2a = SMA_PLSIMread(Ydesign, Xdesign, efitBetas, npcval, trunc, Y=data[,c("Y")])
    basis <- SMA2a$basis
    xmat <- SMA2a$xmat
    zmat <- SMA2a$zmat
    ymat <- SMA2a$ymat

  SMA2b = SMA_PLSIMest(ymat, xmat, zmat, basis)
    EstBeta <- SMA2b$EstBeta

  # SAIE(s)
  EstAlpha = efitBetas[1,,]
  EstSAIE = matrix(c(EstC)*(EstAlpha*EstBeta), L0,1)

  return(EstSAIE)
}
