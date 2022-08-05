# Constraining thetaF and thetaE: ####
countCPsSingleStageNSC <- function(r,N){
  total.CPs <- (r+1)*(N-r)-1
  total.CPs
}

countCPsTwoStageNSC <- function(r1,n1,r,N){
  total.CPs <- (r1+1)*(n1-r1) + (r-r1)*(N-r) - 1
  total.CPs
}

countOrderedPairsSimple <- function(total.CPs, SI.units=FALSE){
  total.ordered.pairs <- t(combn(c(1:total.CPs), 2))
  if(SI.units)  total.ordered.pairs <- format(total.ordered.pairs, scientific=TRUE)
  total.ordered.pairs
}

countOrderedPairsComplex <- function(r, n, p0, p1, thetaFmax=1, thetaEmin=0){
  CPmat <- findCPmatrix(r=r, n=n, Csize=1, p0=p0, p1=p1)
  CP.vec <- sort(unique(c(CPmat)))
  CP.vec <- CP.vec[CP.vec<=thetaFmax | CP.vec>=thetaEmin]
  ordered.pairs <- t(combn(CP.vec, 2))
  ordered.pairs <- ordered.pairs[ordered.pairs[,1]<=thetaFmax & ordered.pairs[,2]>=thetaEmin, ]
  output <- list(CPs=length(CP.vec),
                 OPs=nrow(ordered.pairs)
                 )
  return(output)
  }
