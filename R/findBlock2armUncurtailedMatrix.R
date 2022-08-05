################ Function for finding the uncurtailed CP matrix:
findBlock2armUncurtailedMatrix <- function(n, r, Bsize, pat.cols, prob.vec){
  cpmat <- matrix(3, ncol=2*n, nrow=min(n+r+Bsize+2, 2*n+1))
  rownames(cpmat) <- 0:(nrow(cpmat)-1)
  cpmat[(n+r+2):nrow(cpmat),] <- 1
  cpmat[1:(n+r+1),2*n] <- 0 # Fail at end
  for(i in (n+r+1):1){
    for(j in pat.cols){  # Only look every C patients (no need to look at final col)
      if(i-1<=j){ # Condition: Sm<=m
        cpmat[i,j] <- ifelse(test=j-(i-1) >= n-r+1, yes=0, no=sum(prob.vec*cpmat[i:(i+Bsize), j+Bsize]))
        # IF success is not possible (i.e. [total no. of pats-Xa+Ya-Xb] >= n-r+1), THEN set CP to zero. Otherwise, calculate it based on "future" CPs.
      }
    }
  }
  for(i in 3:nrow(cpmat)){
    cpmat[i, 1:(i-2)] <- NA
  }
  cpmat
}
