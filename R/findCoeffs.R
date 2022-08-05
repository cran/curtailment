findCoeffs <- function(n, p0, p1){
  ##### Matrix of coefficients: Probability of a SINGLE path leading to each point:
  coeffs <- matrix(NA, ncol=n, nrow=n+1)
  coeffs.p0 <- coeffs
  q0 <- 1-p0
  q1 <- 1-p1
  for(i in 1:(n+1)){
    coeffs[i, ] <- p1^(i-1) * q1^((2-i):(2-i+n-1))
    coeffs.p0[i, ] <- p0^(i-1) * q0^((2-i):(2-i+n-1))
  }
  coeffs.list <- list(coeffs.p0=coeffs.p0,
                      coeffs=coeffs)
}
