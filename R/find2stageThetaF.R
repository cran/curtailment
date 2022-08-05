#' Find conditional power for futility
#'
#' Finds maximum conditional power at which trial will stop for futility
#'
#' @param r Number of responses to be exceeded for trial success
#' @param r1 Number of responses to be exceeded to continue to second stage
#' @param n2 Sample size of second stage
#' @param p Response rate under alternative hypothesis H1
#'
#' @return The maximum conditional power at which the trial will stop for futility. This is the CP when there are r1 responses out of n1 participants.
#'
#' @examples
#' find2stageThetaF(r=20, r1=8, n2=13, p=0.92)
#'
#' @export
find2stageThetaF <- function(r, r1, n2, p){
  response.vec <- (r-r1+1):n2
  thetaF <- sum(choose(n2, response.vec) * p^response.vec * (1-p)^(n2-response.vec))
  thetaF
}
