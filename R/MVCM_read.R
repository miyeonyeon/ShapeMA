#' MVCM_read is to read the fiber coordinate data, design matrix for functional linear model, and fiber bundle diffusion properties
#'   Can use first  3, 4 or 5 arguments.
#' Inputs:
#'     tractData: the text file containing (x, y, z) coordinates of all locations on a given fiber tract.
#'                The data set should start from one end to the other end.
#'                 tractData is a L0 x 3 matrix
#'                 L0 denotes the number of locations.
#'                 3 denotes the three coordinates.
#'     designData: the text file containing covariates of interest. Please always include the intercept in the first column.
#'                 is a n x p matrix.
#'                 n denotes the number of subjects.
#'                 p denotes the number of covariates.
#'     diffusionFiles:
#'                 diffusionFiles is a mx1 cell containing the names of all fiber diffusion properties files.
#'                 Each fiber bundle diffusion properties should contain a L_0 x n matrix.
#'                 Rows correspond to the rows in tractData, while columns correspond to the columns in designData.
#'     nofeatures: the number of diffusion properties that we want to joint model.
#'     scalediffusion: a m x 1 vector of scales for each properties
#' Output:
#'     NoSetup    - col vector of [n, L0, p, m]
#'     arclength  - col vector of the arclength from one end to the other end.
#'     Xdesign    - a n x p normalized design matrix.
#'     Ydesign    - a n x L0 x m matrix
#' ###################################################################################################################
#'
#'    Copyright (c) H. T. Zhu  2009
#'    Modifiedy by L. L. Kong 2010
#'
#' @references Hongtu Zhu. Runze Li. Linglong Kong. "Multivariate varying coefficient model for functional responses." Ann. Statist. 40 (5) 2634 - 2666, October 2012.
#' @export

MVCM_read <- function(tractdata, designdata, diffusionFiles, nofeatures) {

  ## Set parameters and defaults according to number of input arguments
  if (nargs() < 3 || nargs() > 5){ # less than 3 argument inputs,  quit
    quit()
  }

  if (nargs() == 3){ # 3 arguments input,  m=1
    nofeatures <- 1
  }

  ## Initializing matrices
  NoSetup <- matrix(0, 4, 1)

  ## Calculating arc length
  tractCON <- tractdata
  if (ncol(tractCON) > 3){
    tractCON <- t(tractCON)
  }

  NoSetup[2, ] <- nrow(tractCON)  # L0
  NoSetup[4, ] <- nofeatures  # m
  arclength <- matrix(0, nrow(tractCON), 1)
  temp1 <- tractCON - tractCON
  temp1[1, ] <- tractCON[1, ]
  temp1[2:nrow(tractCON), ] <- tractCON[1:(nrow(tractCON)-1), ]
  temp1 <- temp1 - tractCON
  for (i in 1:nrow(tractCON)){
    arclength[i, ] <- norm(data.matrix(temp1[i, ]))
    if (i > 1){
      arclength[i, ] <- arclength[i, ] + arclength[i-1, ]
    }
  }

  ## Normalize the design matrix
  Xdesign <- designdata
  if (nrow(Xdesign) < ncol(Xdesign)){
    Xdesign <- t(Xdesign)
  }
  NoSetup[1] <- nrow(Xdesign) # n
  NoSetup[3] <- ncol(Xdesign) # p
  for (i in 1:ncol(Xdesign)){
    if (sd(Xdesign[, i]) > 0){
      Xdesign[, i] <- (Xdesign[, i] - mean(Xdesign[, i])) / sd(Xdesign[, i])
    }
  }

  ## Combine all diffusion properties from the select fiber tract
  Ydesign0 <- array(0, c(nofeatures, NoSetup[2], NoSetup[1]))
  for (i in 1:nofeatures){
    tempt <- diffusionFiles[[i]]
    if (NoSetup[2] == nrow(tempt)){
      Ydesign0[i, , ] <- tempt[, 1:NoSetup[1]]
    } else {
      if (NoSetup[2] == ncol(tempt)){
        tempt <- t(tempt)
        Ydesign0[i, , ] <- tempt[, 1:NoSetup[1]]
      } else {
        exit(1)
      }
    }
  }

  ## Standardize all diffusion properties to a comparable scale
  meanYdata <- numeric(nofeatures)
  for (i in 1:nofeatures){
    meanYdata[i] <- mean(abs(Ydesign0[i, , ]))
  }
  maxMeanYdata <- base::max(meanYdata)
  for (i in 1:nofeatures){
    Ydesign0[i, , ] <- Ydesign0[i, , ] * floor(maxMeanYdata / meanYdata[i])
  }
  scalediffusion <- floor(meanYdata / maxMeanYdata)

  Ydesign <- aperm(Ydesign0, c(3,2,1))

  return(list(NoSetup=NoSetup, arclength=arclength, Xdesign=Xdesign, Ydesign=Ydesign, scalediffusion=scalediffusion))
}
