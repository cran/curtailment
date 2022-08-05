findWald <- function(alpha, power, p0, p1){
  beta <- 1-power
  h0.top <- (1-alpha)*log(beta/(1-alpha))+alpha*log((1-beta)/alpha)
  h0.bottom <- p0*log(p1/p0)+(1-p0)*log((1-p1)/(1-p0))
  Ep0n <- h0.top/h0.bottom
  h1.top <- beta*log(beta/(1-alpha))+(1-beta)*log((1-beta)/alpha)
  h1.bottom <- p1*log(p1/p0)+(1-p1)*log((1-p1)/(1-p0))
  Ep1n <- h1.top/h1.bottom
  round(c(Ep0n, Ep1n),1)
}
