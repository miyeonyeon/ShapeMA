#' Estimate the direct and indirect effects using an unpenalized outcome model
#'
#' @param data a matrix including id, W (clinical confounders), X (genetic exposure), Psi (shape mediators), Y (neurocognitive outcome), and Phi (FPC scores).
#' @param mc_draws the number of Monte Carlo draws. Default is 1e2.
#' @param pt_est an option to get point estimators for the causal effects. Default is TRUE.
#' @param t the number of distinct mediators.
#' @return point estimators for the causal effects. ie.g0 and ie.g1 are the direct and indirect effects, respectively.
#' @references W. W. Loh, B. Moerkerke, T. Loeys, and S. Vansteelandt, “Nonlinear Mediation Analysis With High-Dimensional Mediators Whose Causal Structure Is Unknown,” Biometrics 78, no. 1 (2022): 46–59.
#' @import PLSiMCpp
#' @export

SMA_CEest <- function(data, mc_draws=1e2, pt_est=TRUE, t) {

  ## Fit outcome model
  Wnames <- colnames(data)[grep("W",colnames(data))]
  Phinames <- colnames(data)[grep("Phi",colnames(data))]
  fitY <- PLSiMCpp::plsim.est(data[,c("X",Wnames)], data[,Phinames], matrix(data$Y,nrow(data),1))

  npcval <- length(Phinames)
  Phiobs <- list(data[data$X==0, Phinames], data[data$X==1, Phinames])
  n_X0 <- sum(data$X==0)
  n_X1 <- sum(data$X==1)

  ## Duplicated data for each subject
  dat <- data.table(data[,c("id",Wnames),drop=FALSE])
  setkey(dat)
  xlevels <- dat[,as.data.table(SMA_CEdupdata(t)), by=id]
  setkey(xlevels)
  dat <- merge(dat,xlevels,all.x=TRUE)
  setkey(dat)
  rm(xlevels)

  SamplePhis <- function(mydt, av_mc) {
    # mydt is a data.table for each id containing the relevant columns.
    # av_mc is an option to get average Monte Carlo draws.

    mydt_mc <- mydt[rep(1:nrow(mydt), each=mc_draws)]
    mydt_mc <- cbind(mydt_mc, "mc"=rep(1:mc_draws, times=nrow(mydt)))
    setkey(mydt_mc)

    Phitilde <- data.frame(matrix(NA, nrow=nrow(mydt_mc), ncol=npcval))
    colnames(Phitilde) <- Phinames

    # sample *joint* mediator values
    # x1==0
    x1_0 <- mydt_mc$Marg==0 & mydt_mc$x1==0
    Phitilde[x1_0,] <- Phiobs[[1]][sample(n_X0, size=sum(x1_0), replace=TRUE),]
    rm(x1_0)
    # x1==1
    x1_1 <- mydt_mc$Marg==0 & mydt_mc$x1==1
    Phitilde[x1_1,] <- Phiobs[[2]][sample(n_X1, size=sum(x1_1), replace=TRUE),]
    rm(x1_1)

    # sample *marginal* mediator values (ignore)
    Phitilde[mydt_mc$Marg==1, ] <- 0
    mydt_mc <- cbind(mydt_mc, Phitilde); rm(Phitilde)
    setkey(mydt_mc)

    ## Predict potential outcomes using sampled mediator values
    mydt_mc[, "Y.x" := predict(fitY, as.matrix(mydt_mc)[,c("x0",Wnames),drop=FALSE], as.matrix(mydt_mc)[,Phinames])]
    setkey(mydt_mc)

    if (av_mc==TRUE) {
      # average over MC draws
      mydt_mc_means <- mydt_mc[,lapply(.SD, mean), by=c("Marg","x.i",paste0("x",0:t)), .SDcols="Y.x"]
      setkey(mydt_mc_means)
      return(mydt_mc_means[,Y.x])

    } else {
      return(mydt_mc)
    }
  }

  if (pt_est==TRUE) {
    # average potential outcomes for each subject only
    dat[, "Y.x" := SamplePhis(.SD, av_mc=TRUE), by=id]
    setkey(dat)

    # average over all subjects
    mu_hat <- dat[, lapply(.SD,mean), by=c("Marg","x.i",paste0("x",0:t)),.SDcols="Y.x"]
    setnames(mu_hat,old="Y.x",new="Y")
    setkey(mu_hat)
  }

  # delta transformation matrix for effects
  eff_mat <- matrix(0,nrow=t+4,ncol=nrow(mu_hat))
  eff_mat[1,2:3] <- c(-1,1)
  eff_mat[2,1:2] <- c(-1,1)
  for (s in 1:(t+1)) {
    eff_mat[s+2,c(4,s+4)] <- c(-1,1)
  }
  eff_mat[nrow(eff_mat),] <- eff_mat[2,] - eff_mat[nrow(eff_mat)-1,]
  est <- (eff_mat %*% mu_hat[,Y])[,1]
  names(est) <- c(paste0("g",0:1),paste0("t",1:t),"t_margsum","mu")

  if (pt_est==TRUE) {
    var_est <- NULL
  }

  return(list("pt"=est, "vcov"=var_est))
}
