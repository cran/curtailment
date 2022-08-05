findN1N2R1R2twoarm <- function(nmin, nmax, e1=FALSE){
    nposs <- nmin:nmax
    n1.list <- list()
    n2.list <- list()

    for(i in 1:length(nposs)){
      n1.list[[i]] <- 1:(nposs[i]-1)
      n2.list[[i]] <- nposs[i]-n1.list[[i]]
    }

    # All possibilities together:
    n1 <- rev(unlist(n1.list))
    n2 <- rev(unlist(n2.list))
    n <- n1 + n2
    ns <- cbind(n1, n2, n)

    ################################ FIND COMBNS OF R1 AND R ###############################

    r1.list <- vector("list")
    ns.list <- vector("list")
    for(i in 1:nrow(ns)){
      r1.list[[i]] <- -n1[i]:n1[i] # r1 values: -n1 to n1, for each possible n1
      #ns.list[[i]] <-
    }

    rownames(ns) <- 1:nrow(ns)
    ns <- ns[rep(row.names(ns), sapply(r1.list, length)), ] # duplicate each row so that there are sufficient rows for each r1 value
    ns <- cbind(ns, unlist(r1.list))
    colnames(ns)[4] <- "r1"

    ######### Add possible r values:
    r.list1 <- apply(ns, 1, function(x) {(x["r1"]-x["n2"]):x["n"]})  # r must satisfy r > r1 and r < n. Also, number of responses required in stage 2 (r2-r1) must be at most n2

    how.many.rs <- sapply(r.list1, length)

    row.names(ns) <- 1:nrow(ns)
    ns <- ns[rep(row.names(ns), how.many.rs), ] # duplicate each row a certain number of times
    ns <- cbind(ns, unlist(r.list1))
    colnames(ns)[5] <- "r2"

    ### Finally, add e1 for stopping for benefit:

    if(e1==TRUE)
    {

    } else {
      rownames(ns) <- 1:nrow(ns)
      ns <- data.frame(ns)
    }

    return(ns)
}
