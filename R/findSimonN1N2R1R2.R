#### Part 1: from "find_N1_N2_R1_R_df.R" ####
################################# THIS CODE CREATES ALL COMBNS OF N1, N2, N, R, R

###### OBJECTS
# ns: data frame containing all combinations of n1/n2/n
# n1.list2: Vector of n1's from this data frame
# n2.list2: Vector of n2's from this data frame
# n.list2:  Vector of n's from this data frame
#
# Note that we start from the greatest n and decrease. This is from a "deprecated" idea of mine re. searching,
#   but doesn't affect anything that follows.
#
# r1: List of vectors. Each vector contains all possibilities of r1 for the same row in ns. For example,
#     r1[[i]] is a vector containing all possibilities of r1 for the combination of n1/n2/n in the row ns[i, ].
#
# r: List of list of vectors. Each vector contains all possibilities of r for the same row in ns and single element in r1.
#    For example, r[[1]][[1]] contains all possibilities for r for ns[1, ] and r1[[1]][1], while
#    r[[1]] will contain all possibilities for r for n1[1] and n2[1] for all r1[[1]].
#
#
findSimonN1N2R1R2 <- function(nmin, nmax, e1=TRUE)
{

  nposs <- nmin:nmax

  n1.list <- list()
  n2.list <- list()

  for(i in 1:length(nposs))
  {
    n1.list[[i]] <- 1:(nposs[i]-1)
    n2.list[[i]] <- nposs[i]-n1.list[[i]]
  }

  # All possibilities together:

  n1 <- rev(unlist(n1.list))
  n2 <- rev(unlist(n2.list))

  n <- n1 + n2

  ns <- cbind(n1, n2, n)

  #colnames(ns) <- c("n1", "n2", "n")



  ################################ FIND COMBNS OF R1 AND R
  #
  # For each possible total N, there is a combination of possible n1 and n2's (above).
  #
  # For each possible combination of n1 and n2, there is a range of possible r1's.
  #
  # Note constraints for r1 and r:
  #
  # r1 < n1
  # r1 < r < n
  # r2-r1 <= n2
  #
  #
  # Note: Increasing r1 (or r) increases P(failure), and so decreases power and type I error;
  #       Decreasing r1 (or r) decreases P(failure), and so increases power and type I error.
  #
  # But how does power and type I error change WRT to each individually? Must find out to be able to search more intelligently.
  #



  r1 <- NULL

  for(i in 1:nrow(ns))
  {
    r1 <- c(r1, 0:(n1[i]-1)) # r1 values: 0 to n1-1
  }

  rownames(ns) <- 1:nrow(ns)

  ns <- ns[rep(row.names(ns), n1), ] # duplicate each row n1 times (0 to (n1-1))

  ns <- cbind(ns, r1)



  ######### Add possible r values


  r <- apply(ns, 1, function(x) {(x["r1"]+1):min(x["n"]-1, x["n"]-x["n1"]+x["r1"])})
  # r must satisfy r > r1 and r < n. Also, number of responses required in stage 2 (r2-r1) must be at most n2

  how.many.rs <- sapply(r, length)

  row.names(ns) <- 1:nrow(ns)

  ns <- ns[rep(row.names(ns), how.many.rs), ] # duplicate each row a certain number of times

  ns <- cbind(ns, unlist(r))

  colnames(ns)[5] <- "r"

  ### Finally, add e1 for stopping for benefit:

  if(e1==TRUE)
  {
    # If stopping for benefit, r1 values must run only from 0 to (n1-2), rather than 0 to (n1-1) if not stopping for benefit, otherwise there is no possible value for e1:
    ns <- ns[ns[, "n1"]-ns[, "r1"] >= 2, ]
    e1 <- apply(ns, 1, function(x) {(x["r1"]+1):(x["n1"]-1)})  # e1 must be constrained to the interval [r1+1, n1-1]
    how.many.e1s <- sapply(e1, length)

    row.names(ns) <- 1:nrow(ns)

    ns <- ns[rep(row.names(ns), how.many.e1s), ] # duplicate each row a certain number of times

    row.names(ns) <- 1:nrow(ns)

    ns <- data.frame(ns, e1=unlist(e1))
  } else {
    rownames(ns) <- 1:nrow(ns)
    ns <- data.frame(ns)
  }

  return(ns)
}
