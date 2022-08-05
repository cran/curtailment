findCarstenChenTypeITypeIIRmRows <- function(nr.list, pc, pt, runs, alpha, power, method){
  ########### Function for single row: Carsten #############
  carstenSim <- function(h0, n1, n, a1, r2, pc, pt, runs){

  if(h0==TRUE){
    pt <- pc
  }

  n2 <- n-n1
  nogo <- 0
  go <- 0
  ss <- rep(NA, runs)

  all.pairs <- rbinom(runs*n, 1, prob=pt) - rbinom(runs*n, 1, prob=pc)
  pairs.mat <- matrix(all.pairs, nrow=runs, ncol=n, byrow=TRUE)

  for(i in 1:runs){
    pair <- pairs.mat[i, ]
    successes <- 0
    fails <- 0
    y <- 0
    j <- 1

    ### Stage 1
    while(y<n1 & successes<a1 & fails<n1-a1+1){
      if(pair[j]==1){
        successes <- successes+1
      } else {
        fails <- fails+1
      }
      y <- y+1
      j <- j+1
    }

    ### Trial fails at stage 1:
    if(fails==n1-a1+1){
      nogo <- nogo+1
      ss[i] <- 2*y
    } else {
      ### Trial does not fail at stage 1 -- recruit the remaining participants until curtailment or end:
      while(y<n & successes<r2 & fails<n1+n2-r2+1){
        if(pair[j]==1){
          successes <- successes+1
        } else {
          fails <- fails+1
        }
        y <- y+1
        j <- j+1
      }
      ### Trial fails at stage 2:
      if(fails==n1+n2-r2+1){
        #print("fail at stage 2", q=F)
        nogo <- nogo+1
      } else {
        go <- go+1
      }
      ss[i] <- 2*y
    }
  }
  return(c(n1, n, a1, r2, go/runs, mean(ss)))
  }
########### End of function for single row #############

########### Function for single row: Chen #############
  chenSim <-  function(h0, n1, n, a1, r2, pc, pt, runs){

    n2 <- n-n1

    # Stopping rules for S1 and S2:
    s1.nogo <- n1-a1+1
    s2.go <- n+r2
    s2.nogo <- n-r2+1

    # h0: Set TRUE to estimate type I error and ESS|pt=pc, set to FALSE for power and ESS|pt=pt
    if(h0==TRUE){
      pt <- pc
    }

    # Simulate all successes together, on trt and on control. "Success" means reponse if on trt, non-response if on control:
    trt <- rbinom(n*runs, 1, prob=pt)
    con <- rbinom(n*runs, 1, prob=1-pc)

    ##### Build matrix of successes, both stage 1 and stage 2 #####
    # Allocate pats to trt or control. Note: Balance only required by end of trial.
    alloc <- vector("list", runs)
    n.times.i.minus.1 <- n*((1:runs)-1)
    success.s12 <- matrix(rep(0, 2*n*runs), nrow=runs)

    # TRUE for TREATMENT, FALSE for CONTROL:
    for(i in 1:runs){
      alloc[[i]] <- sample(rep(c(T, F), n), size=2*n, replace=F)
      s.index <- (n.times.i.minus.1[i]+1):(n.times.i.minus.1[i]+n)
      success.s12[i, alloc[[i]]] <-  trt[s.index]
      success.s12[i,!alloc[[i]]] <-  con[s.index]
    }

    success <- success.s12
    failure <- -1*(success-1)

    # Cumulative successes and failures over time:
    success.cum <- t(apply(success, 1, cumsum))
    failure.cum <- t(apply(failure, 1, cumsum))
    # Stage 1 only:
    success.s1.cum <- success.cum[,1:(2*n1)]
    failure.s1.cum <- failure.cum[,1:(2*n1)]

    # Split into "curtailed during S1" and "not curtailed during S1". Note: curtail for no go only.
    curtailed.s1.bin <- apply(failure.s1.cum, 1, function(x) any(x==s1.nogo))
    curtailed.s1.index <- which(curtailed.s1.bin) # Index of trials/rows that reach the S1 no go stopping boundary
    curtailed.s1.subset <- failure.s1.cum[curtailed.s1.index, , drop=FALSE]
    # Sample size of trials curtailed at S1:
    s1.curtailed.ss <- apply(curtailed.s1.subset, 1, function(x) which.max(x==s1.nogo))

    ########## All other trials progress to S2. Subset these:
    success.cum.nocurtail.at.s1 <- success.cum[-curtailed.s1.index, , drop=FALSE]
    failure.cum.nocurtail.at.s1 <- failure.cum[-curtailed.s1.index, , drop=FALSE]

    # Trials/rows that reach the S2 go stopping boundary (including trials that continue to the end):
    s2.go.bin <- apply(success.cum.nocurtail.at.s1, 1, function(x) any(x==s2.go))
    s2.go.index <- which(s2.go.bin)
    # Sample size of trials with a go decision:
    s2.go.ss <- apply(success.cum.nocurtail.at.s1[s2.go.index, , drop=FALSE], 1, function(x) which.max(x==s2.go))
    # Sample size of trials with a no go decision, conditional on not stopping in S1:
    s2.nogo.ss <- apply(failure.cum.nocurtail.at.s1[-s2.go.index, , drop=FALSE], 1, function(x) which.max(x==s2.nogo))

    ess <- sum(s1.curtailed.ss, s2.go.ss, s2.nogo.ss)/runs
    prob.reject.h0 <- length(s2.go.ss)/runs
    prob.accept.h0 <- (length(s2.nogo.ss)+length(s1.curtailed.ss))/runs

    return(c(n1, n, a1, r2, prob.reject.h0, ess))
  }


  chenSimX <- function(h0, n1, n, a1, r2, pc, pt, runs){

    n2 <- n-n1

    # Stopping rules for S1 and S2:
    s1.nogo <- n1-a1+1
    s2.go <- n+r2
    s2.nogo <- n-r2+1

    # h0: Set TRUE to estimate type I error and ESS|pt=pc, set to FALSE for power and ESS|pt=pt
    if(h0==TRUE){
      pt <- pc
    }

    # Simulate all successes together, on trt and on control. "Success" means reponse if on trt, non-response if on control:
    trt <- rbinom(n*runs, 1, prob=pt)
    con <- rbinom(n*runs, 1, prob=1-pc)

    ##### Build matrix of successes, both stage 1 and stage 2 #####
    # Allocate pats to trt or control. Note: Balance only required by end of trial.
    alloc <- vector("list", runs)
    n.times.i.minus.1 <- n*((1:runs)-1)
    success.s12 <- matrix(rep(0, 2*n*runs), nrow=runs)

    # TRUE for TREATMENT, FALSE for CONTROL:
    for(i in 1:runs){
      alloc[[i]] <- sample(rep(c(T, F), n), size=2*n, replace=F)
      s.index <- (n.times.i.minus.1[i]+1):(n.times.i.minus.1[i]+n)
      success.s12[i, alloc[[i]]] <-  trt[s.index]
      success.s12[i,!alloc[[i]]] <-  con[s.index]
    }

    success <- success.s12
    failure <- -1*(success-1)

    # Cumulative successes and failures over time:
    success.cum <- t(apply(success, 1, cumsum))
    failure.cum <- t(apply(failure, 1, cumsum))

    # Stage 1 only:
    success.s1.cum <- success.cum[,1:(2*n1)]
    failure.s1.cum <- failure.cum[,1:(2*n1)]

    # Which trials fail early due to hitting n1-a1+1 failures?
    fail.wrt.s1.boundary.index <- apply(failure.s1.cum, 1, function(x) any(x==s1.nogo))
    curtailed.s1.subset <- failure.s1.cum[fail.wrt.s1.boundary.index, , drop=FALSE]
    # Sample size of trials curtailed at S1:
    s1.curtailed.ss <- apply(curtailed.s1.subset, 1, function(x) which.max(x==s1.nogo))

    # Subset cumulative results to only the remaining trials:
    success.cum.non.s1.failure.trials.subset <- success.cum[-fail.wrt.s1.boundary.index, , drop=FALSE]
    failure.cum.non.s1.failure.trials.subset <- failure.cum[-fail.wrt.s1.boundary.index, , drop=FALSE]
    # Of these remaining trials, which trials reach the req'd no. of successes, and after how many participants?
    successful.trial.index <- apply(success.cum.non.s1.failure.trials.subset, 1, function(x) any(x==s2.go))
    success.cum.successful.trials.subset <- success.cum.non.s1.failure.trials.subset[successful.trial.index, , drop=FALSE]
    early.stop.success.ss <- apply(success.cum.successful.trials.subset, 1, function(x) which.max(x==s2.go))
    # The remaining trials must fail by reaching the required no. of failures. After how many participants do the trials end?
    failure.cum.s2.failure.trials.subset <- failure.cum.non.s1.failure.trials.subset[!successful.trial.index, , drop=FALSE]
    early.stop.failure.ss <- apply(failure.cum.s2.failure.trials.subset, 1, function(x) which.max(x==s2.nogo))

    # Find type I error, power, ESS:
    prob.reject.h0 <- length(early.stop.success.ss)/runs
    ess <- 0.5*sum(s1.curtailed.ss, early.stop.success.ss, early.stop.failure.ss)/runs



#
#     # Split into "curtailed during S1" and "not curtailed during S1". Note: curtail for no go only.
#     curtailed.s1.bin <- apply(failure.s1.cum, 1, function(x) any(x==s1.nogo))
#     curtailed.s1.index <- which(curtailed.s1.bin) # Index of trials/rows that reach the S1 no go stopping boundary
#     curtailed.s1.subset <- failure.s1.cum[curtailed.s1.index, , drop=FALSE]
#     # Sample size of trials curtailed at S1:
#     s1.curtailed.ss <- apply(curtailed.s1.subset, 1, function(x) which.max(x==s1.nogo))
#
#     ########## All other trials progress to S2. Subset these:
#     success.cum.nocurtail.at.s1 <- success.cum[-curtailed.s1.index, , drop=FALSE]
#     failure.cum.nocurtail.at.s1 <- failure.cum[-curtailed.s1.index, , drop=FALSE]
#
#     # Trials/rows that reach the S2 go stopping boundary (including trials that continue to the end):
#     s2.go.bin <- apply(success.cum.nocurtail.at.s1, 1, function(x) any(x==s2.go))
#     s2.go.index <- which(s2.go.bin)
#     # Sample size of trials with a go decision:
#     s2.go.ss <- apply(success.cum.nocurtail.at.s1[s2.go.index, , drop=FALSE], 1, function(x) which.max(x==s2.go))
#     # Sample size of trials with a no go decision, conditional on not stopping in S1:
#     s2.nogo.ss <- apply(failure.cum.nocurtail.at.s1[-s2.go.index, , drop=FALSE], 1, function(x) which.max(x==s2.nogo))
#
#     ess <- sum(s1.curtailed.ss, s2.go.ss, s2.nogo.ss)/runs
#     prob.reject.h0 <- length(s2.go.ss)/runs
#     prob.accept.h0 <- (length(s2.nogo.ss)+length(s1.curtailed.ss))/runs

    return(c(n1, n, a1, r2, prob.reject.h0, ess))
  }
########### End of function for single row ############


output <- vector("list", nrow(nr.list))

# n1, n, a1/r1 are the same for each row of the data frame:
n1 <- nr.list[,"n1"][1]
n <- nr.list[,"n"][1]
a1 <- nr.list[,"r1"][1]
r2.vec <- nr.list[,"r2"]

# Run simulations and keep only {n1,n,a1,r2} combns that are feasible in terms of power:
if(method=="carsten"){
  for(i in 1:nrow(nr.list)){
    output[[i]] <- carstenSim(h0=FALSE, n1=n1, n=n, a1=a1, r2=r2.vec[i], pc=pc, pt=pt, runs=runs)
    if(output[[i]][5]<power){ # stop as soon as power drops below fixed value (ie becomes unfeasible) and remove that row:
      output[[i]] <- NULL
      break
    }
  }
} else {
  for(i in 1:nrow(nr.list)){
    output[[i]] <- chenSim(h0=FALSE, n1=n1, n=n, a1=a1, r2=r2.vec[i], pc=pc, pt=pt, runs=runs)
    if(output[[i]][5]<power){ # stop as soon as power drops below fixed value (ie becomes unfeasible) and remove that row:
      output[[i]] <- NULL
      break
    }
  }
}

output <- do.call(rbind, output)


if(!is.null(output)){
  colnames(output) <- c("n1", "n", "r1", "r2", "pwr", "Ess")
  output <- subset(output, output[,"pwr"]>=power)
# Now type I error:
  typeIoutput <- vector("list", nrow(output))
  if(method=="carsten"){
    for(i in 1:nrow(output)){
      typeIerr <- 0
      # Reverse order of r2 values to start with greatest value and decrease, so that type I error increases the code proceeds:
      reversed.r2 <- rev(output[,"r2"])
      for(i in 1:nrow(output)){
        typeIoutput[[i]] <- carstenSim(h0=TRUE, n1=n1, n=n, a1=a1, r2=reversed.r2[i], pc=pc, pt=pt, runs=runs)
        if(typeIoutput[[i]][5]>alpha){ # stop as soon as type I error increases above fixed value (ie becomes unfeasible)and remove that row:
          typeIoutput[[i]] <- NULL
          break
        }
      }
    }
  }  else{
    for(i in 1:nrow(output)){
      typeIerr <- 0
      # Reverse order of r2 values to start with greatest value and decrease, so that type I error increases the code proceeds:
      reversed.r2 <- rev(output[,"r2"])
      for(i in 1:nrow(output)){
        typeIoutput[[i]] <- chenSim(h0=TRUE, n1=n1, n=n, a1=a1, r2=reversed.r2[i], pc=pc, pt=pt, runs=runs)
        if(typeIoutput[[i]][5]>alpha){ # stop as soon as type I error increases above fixed value (ie becomes unfeasible)and remove that row:
          typeIoutput[[i]] <- NULL
          break
        }
      }
    }
  }
} else{ # If there are no designs with pwr >= power, stop and return NULL:
  return(output)
  }

typeIoutput <- do.call(rbind, typeIoutput)

# If there are feasible designs, merge power and type I error results, o/w stop:
if(!is.null(typeIoutput)){
  colnames(typeIoutput) <- c("n1", "n", "r1", "r2", "typeIerr", "EssH0")
  all.results <- merge(output, typeIoutput, all=FALSE)
} else{
  return(typeIoutput)
}

# # Subset to feasible results:
# subset.results <- all.results[all.results[,"typeIerr"]<=alpha & all.results[,"pwr"]>=power, ]
#
# if(nrow(subset.results)>0){
#   # Discard all "inferior" designs:
#   discard <- rep(NA, nrow(subset.results))
#   for(i in 1:nrow(subset.results)){
#    discard[i] <- sum(subset.results[i, "EssH0"] > subset.results[, "EssH0"] & subset.results[i, "Ess"] > subset.results[, "Ess"] & subset.results[i, "n"] >= subset.results[, "n"])
#    #print(i)
#   }
#   subset.results <- subset.results[discard==0,,drop=FALSE]
# }
#   return(subset.results)
return(all.results)
}










 findSingle2arm2stageJungDesignFast <- function(n1, n2, n, a1, r2, p0, p1, alpha, power){

   #print(paste(n, n1, n2), q=F)

   k1 <- a1:n1

    y1.list <- list()
    for(i in 1:length(k1)){
      y1.list[[i]] <- max(0, -k1[i]):(n1-max(0, k1[i]))
    }

    k1 <- rep(k1, sapply(y1.list, length))
    y1 <- unlist(y1.list)

    combns <- cbind(k1, y1)
    colnames(combns) <- c("k1", "y1")
    rownames(combns) <- 1:nrow(combns)

    k2.list <- vector("list", length(k1))
    for(i in 1:length(k1)){
      k2.list[[i]] <- (r2-k1[i]):n2
    }

    combns2 <- combns[rep(row.names(combns), sapply(k2.list, length)), , drop=FALSE] # duplicate each row so that there are sufficient rows for each a1 value
    k2 <- unlist(k2.list)
    combns2 <- cbind(combns2, k2)
    rownames(combns2) <- 1:nrow(combns2)

    y2.list <- vector("list", length(k2))
    for(i in 1:length(k2)){
      current.k2 <- -k2[i]
      y2.list[[i]] <- max(0, current.k2):(n2-max(0, current.k2))
    }
    y2 <- unlist(y2.list)

    combns3 <- combns2[rep(row.names(combns2), sapply(y2.list, length)), , drop=FALSE] # duplicate each row so that there are sufficient rows for each a1 value
    all.combns <- cbind(combns3, y2)

    # Convert to vectors for speed:

    k1.vec <- all.combns[,"k1"]
    y1.vec <- all.combns[,"y1"]
    k2.vec <- all.combns[,"k2"]
    y2.vec <- all.combns[,"y2"]

    # Easier to understand, but slower:
     part1 <- choose(n1, y1.vec)*p0^y1.vec*(1-p0)^(n1-y1.vec) * choose(n2, y2.vec)*p0^y2.vec*(1-p0)^(n2-y2.vec)
     typeIerr <- sum(part1 * choose(n1, k1.vec+y1.vec)*p0^(k1.vec+y1.vec)*(1-p0)^(n1-(k1.vec+y1.vec)) * choose(n2, k2.vec+y2.vec)*p0^(k2.vec+y2.vec)*(1-p0)^(n2-(k2.vec+y2.vec)))
     pwr <- sum(part1 * choose(n1, k1.vec+y1.vec)*p1^(k1.vec+y1.vec)*(1-p1)^(n1-(k1.vec+y1.vec)) * choose(n2, k2.vec+y2.vec)*p1^(k2.vec+y2.vec)*(1-p1)^(n2-(k2.vec+y2.vec)))

    # Harder to understand, but faster:
    # q0 <- 1-p0
    # q1 <- 1-p1
    # n1.minus.y1 <- n1-y1.vec
    # n2.minus.y2 <- n1-y1.vec
    # k1.plus.y1 <- k1.vec+y1.vec
    # k2.plus.y2 <- k2.vec+y2.vec
    # n1.minus.k1.and.y1 <- n1-k1.plus.y1
    # n2.minus.k2.and.y2 <- n2-k2.plus.y2
    # choose.n1.k1.plus.y1 <- choose(n1, k1.plus.y1)
    # choose.n2.k2.plus.y2 <- choose(n2, k2.plus.y2)
    #
    # part1 <- choose(n1, y1.vec)*p0^y1.vec*q0^n1.minus.y1 * choose(n2, y2.vec)*p0^y2.vec*q0^n2.minus.y2 * choose.n1.k1.plus.y1 * choose.n2.k2.plus.y2
    # typeIerr <- sum(part1 * p0^k1.plus.y1*q0^n1.minus.k1.and.y1 * p0^k2.plus.y2*q0^n2.minus.k2.and.y2)
    # pwr <- sum(part1 * p1^k1.plus.y1*q1^n1.minus.k1.and.y1 * p1^k2.plus.y2*q1^n2.minus.k2.and.y2)



    # Find ESS under H0 and H1:
    if(typeIerr<=alpha & pwr>=power){
      k11 <- -n1:(a1-1)
      y11.list <- vector("list", length(k11))
      for(i in 1:length(k11)){
        y11.list[[i]] <- max(0, -k11[i]):(n1-max(0, k11[i]))
      }
      k11.vec <- rep(k11, sapply(y11.list, length))
      y11.vec <- unlist(y11.list)

      petH0 <- sum(choose(n1, y11.vec)*p0^y11.vec*(1-p0)^(n1-y11.vec) * choose(n1, k11.vec+y11.vec)*p0^(k11.vec+y11.vec)*(1-p0)^(n1-(k11.vec+y11.vec)))
      petH1 <- sum(choose(n1, y11.vec)*p1^y11.vec*(1-p1)^(n1-y11.vec) * choose(n1, k11.vec+y11.vec)*p1^(k11.vec+y11.vec)*(1-p1)^(n1-(k11.vec+y11.vec)))

      # choose.n1.y11 <- choose(n1, y11.vec)
      # n1.minus.y11 <- n1-y11.vec
      # choose.k11.y11 <- choose(n1, k11.vec+y11.vec)
      # k11.y11 <- k11.vec+y11.vec
      # n1.minus.k11.y11 <- n1-k11.y11
      #
      # pet.part1 <- choose.n1.y11 * choose.k11.y11
      # petH0 <- sum(pet.part1 * p0^y11.vec*q0^n1.minus.y11 * p0^k11.y11*q0^n1.minus.k11.y11)
      # petH1 <- sum(pet.part1 * p1^y11.vec*q1^n1.minus.y11 * p1^k11.y11*q1^n1.minus.k11.y11)

      essH0 <- n1*petH0 + n*(1-petH0)
      essH1 <- n1*petH1 + n*(1-petH1)

      return(c(n1, n2, n, a1, r2, typeIerr, pwr, essH0, essH1))
    } else {
        return(c(n1, n2, n, a1, r2, typeIerr, pwr, NA, NA))
      }
 }



######## Find stopping boundaries for one design #########
# The function findBounds can be used to find stopping boundaries for an SC design.







 # Plot rejection regions ####
 # Write a program that will find the sample size using our design and Carsten's design, for a given set of data #

 # for p0=0.1, p1=0.3, alpha=0.15, power=0.8, h0-optimal designs are:
 # Carsten:
 # n1=7; n=16; r1=1; r2=3
 # This design should have the following OCs:
 # Type I error: 0.148, Power: 0.809, EssH0: 21.27400, EssH1: 19.44200
 #
 # Our design (not strictly H0-optimal, but is within 0.5 of optimal wrt EssH0 and has N=31, vs N=71 for the actual H0-optimal):
 # n=31; r2=3; thetaF=0.1277766; thetaE=0.9300000

 findRejectionRegions <- function(n, r, thetaF=NULL, thetaE=NULL, pc, pt, method=NULL, Bsize=NULL){
   ######################## FUNCTION TO FIND BASIC CP MATRIX (CP=0/1/neither) FOR CARSTEN
   carstenCPonly <- function(n1, n, r1, r, pair=FALSE){
     cpmat <- matrix(0.5, ncol=2*n, nrow=r+1)
     rownames(cpmat) <- 0:(nrow(cpmat)-1)
     cpmat[r+1, (2*r):(2*n)] <- 1 # Success is a PAIR of results s.t. Xt-Xc=1
     cpmat[1:r,2*n] <- 0 # Fail at end
     pat.cols <- seq(from=(2*n)-2, to=2, by=-2)
     for(i in nrow(cpmat):1){
       for(j in pat.cols){  # Only look every C patients (no need to look at final col)
         if(i-1<=j/2){ # Condition: number of successes must be <= pairs of patients so far
           # IF success is not possible (i.e. [total no. of failures] > n1-r1+1 at stage 1 or [total no. of failures] > n-r+1), THEN set CP to zero:
           if((i<=r1 & j/2 - (i-1) >= n1-r1+1) | (j/2-(i-1) >= n-r+1)){
             cpmat[i,j] <- 0
           }
         } else{
           cpmat[i,j] <- NA # impossible to have more successes than pairs
         }
       }
     }
     if(pair==TRUE){
       cpmat <- cpmat[, seq(2, ncol(cpmat), by=2)]
     }
     return(cpmat)
   }

   #### Find CPs for block and Carsten designs:
   if(method=="block"){
     block.mat <- findBlockCP(n=n, r=r, Bsize=Bsize, pc=pc, pt=pt, thetaF=thetaF, thetaE=thetaE)
     #### Find lower and upper stopping boundaries for block and Carsten designs:
     lower <- rep(-Inf, ncol(block.mat)/2)
     upper <- rep(Inf, ncol(block.mat)/2)
     looks <- seq(2, ncol(block.mat), by=2)
     for(j in 1:length(looks)){
       if(any(block.mat[,looks[j]]==0, na.rm = T)){
         lower[j] <- max(which(block.mat[,looks[j]]==0))-1 # Minus 1 to account for row i == [number of responses-1]
       }
       if(any(block.mat[,looks[j]]==1, na.rm = T)){
         upper[j] <- which.max(block.mat[,looks[j]])-1 # Minus 1 to account for row i == [number of responses-1]
       }
     }

     m.list <- vector("list", max(looks))
     for(i in 1:length(looks)){
       m.list[[i]] <- data.frame(m=looks[i], successes=0:max(looks), outcome=NA)
       m.list[[i]]$outcome[m.list[[i]]$successes <= lower[i]] <- "No go"
       m.list[[i]]$outcome[m.list[[i]]$successes >= lower[i]+1 & m.list[[i]]$successes <= upper[i]-1] <- "Continue"
       m.list[[i]]$outcome[m.list[[i]]$successes >= upper[i] & m.list[[i]]$successes <= looks[i]] <- "Go"
       m.list[[i]]$outcome[m.list[[i]]$successes > looks[i]] <- NA
     }
   }

   if(method=="carsten"){
     carsten.mat <- carstenCPonly(n1=n[[1]], n=n[[2]], r1=r[[1]], r=r[[2]], pair=F)
     #### Find lower and upper stopping boundaries for block and Carsten designs:
     lower.carsten <- rep(-Inf, ncol(carsten.mat)/2)
     upper.carsten <- rep(Inf, ncol(carsten.mat)/2)
     looks.carsten <- seq(2, ncol(carsten.mat), by=2)
     for(j in 1:length(looks.carsten)){
       if(any(carsten.mat[,looks.carsten[j]]==0, na.rm = T)){
         lower.carsten[j] <- max(which(carsten.mat[,looks.carsten[j]]==0))-1 # Minus 1 to account for row i == [number of responses-1]
       }
       if(any(carsten.mat[,looks.carsten[j]]==1, na.rm = T)){
         upper.carsten[j] <- which.max(carsten.mat[,looks.carsten[j]])-1 # Minus 1 to account for row i == [number of responses-1]
       }
     }
     looks <- looks.carsten
     lower <- lower.carsten
     upper <- upper.carsten
     m.list <- vector("list", max(looks))
     for(i in 1:length(looks)){
       m.list[[i]] <- data.frame(m=looks[i], successes=0:(max(looks)/2), outcome=NA)
       m.list[[i]]$outcome[m.list[[i]]$successes <= lower[i]] <- "No go"
       m.list[[i]]$outcome[m.list[[i]]$successes >= lower[i]+1 & m.list[[i]]$successes <= upper[i]-1] <- "Continue"
       m.list[[i]]$outcome[m.list[[i]]$successes >= upper[i] & m.list[[i]]$successes <= looks[i]/2] <- "Go"
       m.list[[i]]$outcome[m.list[[i]]$successes > looks[i]/2] <- NA
     }
   }
   m.df <- do.call(rbind, m.list)
   m.df <- m.df[!is.na(m.df$outcome), ]
   return(m.df)
 }
