#' Find a single Simon design
#'
#' This function finds the operating characteristics of a single realisation of a Simon design. It returns
#' the inputted values along with the type-I error-rate (alpha), power,  expected sample size under p=p0 (EssH0)
#' and expected sample size under p=p1 (Ess).
#'
#' @param n1 Number of participants in stage 1
#' @param n2 Number of participants in stage 2
#' @param r1 Interim stopping boundary that must be exceeded to avoid no go stopping
#' @param r Final rejection boundary that must be exceeded to reject the null hypothesis.
#' @param p0 Anticipated response probability for inefficacious treatment
#' @param p1 Anticipated response probability for efficacious treatment
#' @return A vector containing all the inputted values and corresponding operating characteristics.
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @references
#' \doi{10.1016/0197-2456(89)90015-9}{Richard Simon,
#' Optimal two-stage designs for phase II clinical trials,
#' Controlled Clinical Trials,
#' Volume 10, Issue 1,
#' 1989,
#' Pages 1-10}
#' @examples findSingleSimonDesign(n1 = 15, n2 = 11, r1 = 1, r = 4, p0 = 0.1, p1 = 0.3)
#' @export
findSingleSimonDesign <- function(n1, n2, r1, r, p0, p1)
  {
  n <- n1+n2
  # Create Pascal's triangle for S1: these are the coefficients (before curtailment) A, where A * p^b * q*c
  pascal.list.s1 <- list(1)
  for (i in 2:(n1+1)) pascal.list.s1[[i]] <- c(0, pascal.list.s1[[i-1]]) + c(pascal.list.s1[[i-1]], 0)
  pascal.list.s1[[1]] <- NULL
  # For Simon's design, only need the final line:
  pascal.list.s1 <- pascal.list.s1[n1]

  # Curtail at n1 only:
  curtail.index <- 1:(r1+1)
  curtail.coefs.s1 <- pascal.list.s1[[1]][curtail.index] # from k=0 to k=r1

  # Use final column from S1:
  pascal.list.s2 <- pascal.list.s1
  pascal.list.s2[[1]][curtail.index] <- 0

  for (i in 2:(n2+1)) pascal.list.s2[[i]] <- c(0, pascal.list.s2[[i-1]]) + c(pascal.list.s2[[i-1]], 0)
  pascal.list.s2[[1]] <- NULL

  # Now obtain the rest of the probability -- the p^b * q^c :
  # S1
  q1 <- 1-p1
  coeffs <- p1^(0:n1)*q1^(n1:0)
  coeffs <- coeffs[curtail.index]

  q0 <- 1-p0
  coeffs.p0 <- p0^(0:n1)*q0^(n1:0)
  coeffs.p0 <- coeffs.p0[curtail.index]

  # Multiply the two vectors (A and p^b * q^c):
  prob.curt.s1 <- curtail.coefs.s1*coeffs

  # for finding type I error prob:
  prob.curt.s1.p0 <- curtail.coefs.s1*coeffs.p0

  # The (S1) curtailed paths:
  k.curt.s1 <- 0:r1
  n.curt.s1 <- rep(n1, length(k.curt.s1))
  curtail.s1 <- cbind(k.curt.s1, n.curt.s1, prob.curt.s1, prob.curt.s1.p0)


  ############## S2

  # Pick out the coefficients for the S2 paths (A, say):
  s2.index <- (r1+2):(n+1)
  curtail.coefs.s2 <- pascal.list.s2[[n2]][s2.index]

  # Now obtain the rest of the probability -- the p^b * q^c :
  coeffs.s2 <- p1^(0:n)*q1^(n:0)
  coeffs.s2 <- coeffs.s2[s2.index]

  coeffs.s2.p0 <- p0^(0:n)*q0^(n:0)
  coeffs.s2.p0 <- coeffs.s2.p0[s2.index]

  # Multiply the two vectors (A and p^b * q^c):
  prob.go <- curtail.coefs.s2*coeffs.s2

  # for finding type I error prob:
  prob.go.p0 <- curtail.coefs.s2*coeffs.s2.p0

  # Paths that reach the end:
  k.go <- (r1+1):n
  n.go <- rep(n, length(k.go))

  go <- cbind(k.go, n.go, prob.go, prob.go.p0)

  final <- rbind(curtail.s1, go)

  ############## WRAPPING UP THE RESULTS

  output <- data.frame(k=final[,1], n=final[,2], prob=final[,3], prob.p0=final[,4])
  output$success <- "Fail"
  output$success[output$k > r] <- "Success"
  power <- sum(output$prob[output$success=="Success"])
  sample.size.expd <- sum(output$n*output$prob)
  sample.size.expd.p0 <- sum(output$n*output$prob.p0)
  alpha <- sum(output$prob.p0[output$success=="Success"])
  to.return <- c(n1=n1, n2=n2, n=n, r1=r1, r=r, alpha=alpha, power=power, EssH0=sample.size.expd.p0, Ess=sample.size.expd)
  to.return
}
