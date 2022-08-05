# Discard dominated designs:
rmDominatedDesigns <- function(data, ess0="EssH0", ess1="Ess", n="n"){
  discard <- rep(NA, nrow(data))
  for(i in 1:nrow(data)){
    discard[i] <- sum(data[i, ess0] > data[, ess0] & data[i, ess1] > data[, ess1] & data[i, n] > data[, n])
    }
  admissible.designs <- data[discard==0, ]
  admissible.designs
}
