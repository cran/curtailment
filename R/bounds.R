bounds <- function(des, type=c("simon", "simon_e1", "nsc", "sc", "mstage"), p0, p){
  n1 <- des[["n1"]]
  n2 <- des[["n2"]]
  r1 <- des[["r1"]]
  r <- des[["r"]]
  N <- des[["n"]]
  alpha <- des[["alpha"]]
  beta <- 1-des[["power"]]

  if(type=="simon"){
    result <- list(J=2, n=c(n1, n2), N=N, a=c(r1, r), r=c(Inf, r+1), alpha=alpha, beta=beta, n1=n1) # Or N+1 instead of Inf -- either way, no stopping
  }
  if(type=="simon_e1"){
    result <- list(J=2, n=c(n1, n2), N=N, a=c(r1, r), r=c(des$e1, r+1), alpha=alpha, beta=beta, n1=n1)
  }
  if(type=="nsc"){
    eff.bound.Sm <- c(rep(Inf, r), rep(r+1, N-r)) # Stop for efficacy at r+1, at every point. Also, or N+1 instead of Inf -- either way, no stopping
    # Stage 1
    fut.bound.s1.m.partial <- (n1-r1):n1 # Values of m where S1 stopping for futility is possible.
    fut.bound.s1.m <- 1:n1
    fut.bound.s1.Sm <- c(rep(-Inf, fut.bound.s1.m.partial[1]-1), seq(from=0, length=length(fut.bound.s1.m.partial))) # or 0 instead of -Inf -- either way, no stopping
    # Stage 2:
    fut.bound.s2.m.partial <- max(N+r1-r+1, n1+1):N
    fut.bound.s2.m <- (n1+1):N
    fut.bound.s2.Sm <- c(rep(-Inf, length(fut.bound.s2.m)-length(fut.bound.s2.m.partial)), seq(from=r1+1, length=length(fut.bound.s2.m.partial))) # or 0 instead of -Inf -- either way, no stopping
    # Output
    result <- list(J=2, n=c(n1, n2), N=N, a=c(r1, r), r=c(Inf, r+1), a_curt=c(fut.bound.s1.Sm, fut.bound.s2.Sm), r_curt=eff.bound.Sm, alpha=alpha, beta=beta)
  }

  if(type=="sc"){
    result <- sc_bounds(n1=n1, n2=des$n2, r1=r1, r2=r, p0=p0, p=p, thetaF=des$thetaF, thetaE=des$thetaE, alpha=des$alpha, beta=1-des$power)
  }
  if(type=="mstage"){
    result <- mstage_bounds(n=des$n, r=des$r, thetaF=des$thetaF, thetaE=des$thetaE, p0=p0, p=p, power=des$power, alpha=des$alpha)
  }
  result
}
