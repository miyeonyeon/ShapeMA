#' Create duplicated data for one subject
#'
#' @param t the number of distinct mediators.
#' @return a (t+5) x (t+3) matrix of duplicated dataset.
#' @references  W. W. Loh, B. Moerkerke, T. Loeys, and S. Vansteelandt, “Nonlinear Mediation Analysis With High-Dimensional Mediators Whose Causal Structure Is Unknown,” Biometrics 78, no. 1 (2022): 46–59.
#' @export

SMA_CEdupdata <- function(t) {

  out <- diag(1, nrow=t, ncol=t)
  out <- rbind(0, out, 1)
  out <- cbind(0, out)
  out <- rbind(0, out[nrow(out),], 1, out)
  colnames(out) <- paste0("x", 0:t)
  out <- cbind("Marg"=c(rep(0,3), rep(1,nrow(out)-3)), "x.i"=1:nrow(out), out)

  return(out)
}
