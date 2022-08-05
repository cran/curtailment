# Simon's design: No curtailment -- only stopping is at end of S1:
simonEfficacy <- function(n1, n2, r1, r, e1, p0, p1)
{

  n <- n1+n2

  # Create Pascal's triangle for S1: these are the coefficients (before curtailment) A, where A * p^b * q*c
  pascal.list.s1 <- list(1)
  for (i in 2:(n1+1)) pascal.list.s1[[i]] <- c(0, pascal.list.s1[[i-1]]) + c(pascal.list.s1[[i-1]], 0)
  pascal.list.s1[[1]] <- NULL
  # For Simon's design, only need the final line:
  pascal.list.s1 <- pascal.list.s1[n1]

  # Curtail at n1 only:
  curtail.index <- c(1:(r1+1), (e1+2):(n1+1)) # We curtail at these indices -- Sm=[0, r1] and Sm=[e1+1, n1] (remembering row 1 is Sm=0)
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
  k.curt.s1 <- c(0:r1, (e1+1):n1)
  n.curt.s1 <- rep(n1, length(k.curt.s1))
  curtail.s1 <- cbind(k.curt.s1, n.curt.s1, prob.curt.s1, prob.curt.s1.p0)


  ############## S2 ###############

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


  ############## WRAPPING UP THE RESULTS ##############

  output <- data.frame(k=final[,1], n=final[,2], prob=final[,3], prob.p0=final[,4])

  output$success <- "Fail"
  output$success[output$k > r] <- "Success"
  output$success[output$n==n1 & output$k > e1] <- "Success"


  # Pr(early termination):
  #PET <- sum(output$prob[output$n < n])
  #PET.p0 <- sum(output$prob.p0[output$n < n])

  power <- sum(output$prob[output$success=="Success"])

  #output$obsd.p <- output$k/output$n

  #output$bias <- output$obsd.p - p

  #bias.mean <- wtd.mean(output$bias, weights=output$prob, normwt=TRUE)
  #bias.var <- wtd.var(output$bias, weights=output$prob, normwt=TRUE)

  #sample.size <- wtd.quantile(output$n, weights=output$prob, normwt=TRUE, probs=c(0.25, 0.5, 0.75))
  #sample.size.expd <- wtd.mean(output$n, weights=output$prob, normwt=TRUE)
  sample.size.expd <- sum(output$n*output$prob)

  #sample.size.p0 <- wtd.quantile(output$n, weights=output$prob.p0, normwt=TRUE, probs=c(0.25, 0.5, 0.75))
  #sample.size.expd.p0 <- wtd.mean(output$n, weights=output$prob.p0, normwt=TRUE)
  sample.size.expd.p0 <- sum(output$n*output$prob.p0)


  alpha <- sum(output$prob.p0[output$success=="Success"])

  #output <- list(output, mean.bias=bias.mean, var.bias=bias.var, sample.size=sample.size, expd.sample.size=sample.size.expd, PET=PET,
  #               sample.size.p0=sample.size.p0, expd.sample.size.p0=sample.size.expd.p0, PET.p0=PET.p0, alpha=alpha, power=power)
  to.return <- c(n1=n1, n2=n2, n=n, r1=r1, r=r, alpha=alpha, power=power, EssH0=sample.size.expd.p0, Ess=sample.size.expd, e1=e1)
  to.return
}
