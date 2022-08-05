findN1N2R1R2NSC <- function(nmin, nmax)
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

  n1.list2 <- rev(unlist(n1.list))
  n2.list2 <- rev(unlist(n2.list))
  length(n1.list2)
  # There are approx. 3000 combinations of n1 and n2.

  ns <- cbind(n1.list2, n2.list2)

  n.list2 <- apply(ns, 1, sum)

  ns <- cbind(ns, n.list2)

  colnames(ns) <- c("n1", "n2", "n")



  ########################################################################################
  ################################ FIND COMBNS OF R1 AND R ###############################
  ########################################################################################
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


  r1 <- list()

  for(i in 1:nrow(ns))
  {
    r1[[i]] <- 0:(ns[i,"n1"]-1) # r1 values: 0 to n1-1
  }

  # How many combinations so far?
  length(unlist(r1)) # Approx. 85,000 designs, before r is varied.

  # For each possible combination of n1, n2 and r1, there is a combination of possible r's:

  r <- vector("list", (nrow(ns)))
  # ^ r will be a list of lists: e.g., r[[1]] will contain all possibilities for r for n1[1] and n2[1],
  # while r[[1]][[1]] will contain all possibilities for r for n1[1] and n2[1] AND r1[[1]][1]

  for(i in 1:nrow(ns))
  {
    #print(i)
    for(j in 1:length(r1[[i]]))
    {
      r[[i]][[j]] <- (r1[[i]][j]+1):(ns[i,"n"]-1)  # r must satisfy r > r1 and r < n

      keep.these <- r[[i]][[j]]-r1[[i]][j] <= ns[i,"n2"] # <- Note constraint: The number of responses required in stage 2 (r2-r1) must be at most n2:
      r[[i]][[j]] <- r[[i]][[j]][keep.these]
    }
  }


  #length(unlist(r)) # Approx. 300,000 designs for n [2, 50]
  # Approx. 200,000 designs for n [61, 65] alone

  ######### REMOVE UNNEEDED OBJECTS:
  rm(nposs, keep.these)
  return(list(ns=ns, n1=n1.list2, n2=n2.list2, n=n.list2, r1=r1, r=r))
}
