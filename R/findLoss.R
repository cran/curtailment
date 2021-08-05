#' Find loss scores for a set of single-arm designs
#'
#' This function finds loss scores for single-arm designs, in particular, designs
#' found using the function singlearmdesigns(). Weights w0 and w1 are chosen,
#' and the loss for each design is found using the equation
#' \emph{loss = w0\*ESS(0) + w1\*ESS(1) + (1-w0-w1)\*N}
#' @param main.output Output from function singlearmDesigns().
#' @param w0 Choice of weight on ESS(0)
#' @param w1 Choice of weight on ESS(0)
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @return Output is a list identical to the output from function singlearmDesigns(),
#'  but with absolute and ranked loss scores added.
#' @examples output <- singlearmDesign(nmin = 30,
#'  nmax = 30,
#'  C = 5,
#'  p0 = 0.1,
#'  p1 = 0.4,
#'  power = 0.8,
#'  alpha = 0.05)
#'  output.w.loss <- findLoss(main.output=output,
#'  w0=0.2,
#'  w1=0.3)
#'@export
findLoss <- function(main.output, w0, w1){
  weight.vec <- data.frame(w0, w1, 1-w0-w1)
  loss <- apply(main.output$all.des, 1, findLossSingleDesign, weight.vec)
  loss.ranked <- rank(loss)
  loss.rounded <- round(loss, 1)
  main.output$all.des <- cbind(main.output$all.des, loss.rounded, loss.ranked)
  main.output
}
