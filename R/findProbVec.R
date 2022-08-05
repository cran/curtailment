################ Function for finding the Prob(reponses on treatment + non-responses on control)=0, 1, 2,... Bsize:
findProbVec <- function(Bsize, pt, qt, pc, qc){
  prob.vec <- rep(NA, Bsize+1)
  for(i in 1:(Bsize+1)){
    positives <- i-1
    full.vec <- expand.grid(rep(list(0:1), Bsize))
    positive.mat <- full.vec[rowSums(full.vec) == positives,]
    negative.mat <- -1*(positive.mat-1)

    positive.vec <- rep(c(pt,qc), each=Bsize/2)
    negative.vec <- rep(c(qt,pc), each=Bsize/2)

    posneg.mat <- t(t(positive.mat)*positive.vec) + t(t(negative.mat)*negative.vec)
    prob.vec[i] <- sum(apply(posneg.mat, 1, prod))
  }
  if(sum(prob.vec)-1 > 1e-8) stop("Probabilities do not sum to 1.")
  prob.vec
}
