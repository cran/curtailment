findDesignsGivenCohortStage <- function(nmin,
                        nmax,
                        C=NA,
                        stages=NA,
                        p0,
                        p1,
                        alpha,
                        power,
                        minstop,
                        maxthetaF=p1,
                        minthetaE=p1,
                        bounds="wald",
                        return.only.admissible=TRUE,
                        max.combns=1e6,
                        maxthetas=NA,
                        fixed.r=NA,
                        exact.thetaF=NA,
                        exact.thetaE=NA,
                        progressBar=FALSE){
  grp <- V2 <- NULL
  q0 <- 1-p0
  q1 <- 1-p1

use.stages <- is.na(C)
if(!is.na(C) & !is.na(stages)) stop("Values given for both cohort/block size C and number of stages. Please choose one only.")

n.initial.range <- nmin:nmax
if(use.stages==TRUE){
  nposs.min.index <- min(which((nmin:nmax)%%stages==0)) # nmin must be a multiple of number of stages
  nposs.max.index <- max(which((nmin:nmax)%%stages==0)) # nmax must be a multiple of number of stages
  if(any(is.infinite(nposs.min.index), is.infinite(nposs.max.index))) stop("There must be at least one value within [nmin, nmax] that is a multiple of the number of stages")
  nposs.min <- n.initial.range[nposs.min.index]
  nposs.max <- n.initial.range[nposs.max.index]
  nposs <- seq(from=nposs.min, to=nposs.max, by=stages)
}else{
  nposs.min.index <- min(which((nmin:nmax)%%C==0)) # nmin must be a multiple of cohort size
  nposs.max.index <- max(which((nmin:nmax)%%C==0)) # nmax must be a multiple of cohort size
  if(any(is.infinite(nposs.min.index), is.infinite(nposs.max.index))) stop("There must be at least one value within [nmin, nmax] that is a multiple of the cohort size")
  nposs.min <- n.initial.range[nposs.min.index]
  nposs.max <- n.initial.range[nposs.max.index]
  nposs <- seq(from=nposs.min, to=nposs.max, by=C)
}



  r <- list()
  for(i in 1:length(nposs))
  {
    r[[i]] <- 0:(nposs[i]-1) # r values: 0 to nposs[i]-1
  }

  ns <- NULL

  for(i in 1:length(nposs)){
    ns <- c(ns, rep(nposs[i], length(r[[i]])))
  }

  sc.subset <- data.frame(n=ns, r=unlist(r))

  if(!is.na(bounds)){
    # Incorporate A'Hern's bounds:
    if(bounds=="ahern")  {
      sc.subset <- sc.subset[sc.subset$r >= p0*sc.subset$n & sc.subset$r <= p1*sc.subset$n, ]
    }

    if(bounds=="wald"){
    # Faster to incorporate Wald's bounds:
      denom <- log(p1/p0) - log((1-p1)/(1-p0))
      accept.null <- log((1-power)/(1-alpha)) / denom  + nposs * log((1-p0)/(1-p1))/denom
      accept.null <- floor(accept.null)

      reject.null <- log((power)/alpha) / denom  + nposs * log((1-p0)/(1-p1))/denom
      reject.null <- ceiling(reject.null)

      r.wald <- NULL
      ns.wald <- NULL
      for(i in 1:length(nposs)){
        r.wald <- c(r.wald, accept.null[i]:reject.null[i])
        ns.wald <- c(ns.wald, rep(nposs[i], length(accept.null[i]:reject.null[i])))
      }
      sc.subset <- data.frame(n=ns.wald, r=r.wald)
    }
  }

  # In case user wishes to specify values for r:
  if(!is.na(fixed.r)){
    sc.subset <- sc.subset[sc.subset$r %in% fixed.r,]
  }

  sc.subset <- sc.subset[sc.subset$r>=0, ]
  sc.subset <- sc.subset[sc.subset$r<sc.subset$n, ]

  if(use.stages==TRUE){
    sc.subset$C <- sc.subset$n/stages
  } else {
    sc.subset$C <- C
  }

  ###### Find thetas for each possible {r, N} combn:
  mat.list <- vector("list", nrow(sc.subset))

  sc.subset$eff.minstop <- NA
  for(i in 1:nrow(sc.subset)){
    C.vec <- seq(sc.subset$C[i], nposs.max, by=sc.subset$C[i]) # all possible TPs, ignoring minstop
    sc.subset$eff.minstop[i] <- C.vec[which.max(C.vec>=minstop)] # minstop must be multiple of C
    mat.list[[i]] <- findCPmatrix(n=sc.subset[i,"n"], r=sc.subset[i,"r"], Csize=sc.subset[i,"C"], p0=p0, p1=p1, minstop=sc.subset[i,"eff.minstop"])
  }

  store.all.thetas <- lapply(mat.list, function(x) {sort(unique(c(x))[unique(c(x)) <= 1])})

if(!is.na(max.combns) & is.na(maxthetas)){ # Only use max.combns if maxthetas is not specified.
  maxthetas <- sqrt(2*max.combns)
}

  if(is.na(exact.thetaF)){
    for(i in 1:nrow(sc.subset)){
      current.theta.vec <- store.all.thetas[[i]]
      store.all.thetas[[i]] <- current.theta.vec[current.theta.vec >= minthetaE | current.theta.vec <= maxthetaF] # Subset thetas BEFORE thinning
      while(length(store.all.thetas[[i]]) > maxthetas){
        every.other.element <- rep(c(FALSE, TRUE), 0.5*(length(store.all.thetas[[i]])-2))
        store.all.thetas[[i]] <-  store.all.thetas[[i]][c(TRUE, every.other.element, TRUE)]
      }
    }
    all.theta.combns <- lapply(store.all.thetas, function(x) {t(combn(x, 2)) })
    for(i in 1:nrow(sc.subset)){
      all.theta.combns[[i]] <- all.theta.combns[[i]][all.theta.combns[[i]][,2] >= minthetaE & all.theta.combns[[i]][,1] <= maxthetaF, ] # Remove ordered pairs where both are low or high
      if(length(all.theta.combns[[i]])==2) all.theta.combns[[i]] <- matrix(c(0, 0.999, 0, 1), nrow=2, byrow=T) # To avoid a crash. See (*) below
      all.theta.combns[[i]] <- all.theta.combns[[i]][order(all.theta.combns[[i]][,2], decreasing=T),]
    }
    # (*) If the only thetas for a design are 0 and 1, then some code below will crash. This could be resolved by adding an if or ifelse statement,
    # but probably at unnecessary computational expense. We deal with these exceptions here, by adding a "false" value for theta.
    }else{ # if exact thetas are given (to speed up a result check):
      for(i in 1:length(store.all.thetas)){
        keep <- abs(store.all.thetas[[i]]-exact.thetaF)<1e-3 | abs(store.all.thetas[[i]]-exact.thetaE)<1e-3
        store.all.thetas[[i]] <- store.all.thetas[[i]][keep]
      }
      all.theta.combns <- lapply(store.all.thetas, function(x) {t(combn(x, 2)) })
      }

  ##### Matrix of coefficients: Probability of a SINGLE path leading to each point:
  coeffs <- matrix(NA, ncol=nmax, nrow=nmax+1)
  coeffs.p0 <- coeffs

  for(i in 1:(nmax+1)){
    coeffs[i, ] <- p1^(i-1) * q1^((2-i):(2-i+nmax-1))
    coeffs.p0[i, ] <- p0^(i-1) * q0^((2-i):(2-i+nmax-1))
  }


  h.results.list <- vector("list", nrow(sc.subset)) #

  if(progressBar==TRUE)  pb <- txtProgressBar(min = 0, max = nrow(sc.subset), style = 3)


  # Now, find the designs, looping over each possible {r, N} combination, and within each {r, N} combination, loop over all combns of {thetaF, thetaE}:
  for(h in 1:nrow(sc.subset)){
    k <- 1
    h.results <- vector("list", nrow(all.theta.combns[[h]]))

    current.theta.combns <- all.theta.combns[[h]] # Take just the thetaF/1 combns for that design.
    # Don't need a matrix of all thetaF and thetaE combns -- potentially quicker to have thetaF as a vector (already defined), and the list can be a list of thetaE vectors:
    current.theta.combns <- data.table::data.table(current.theta.combns)
    current.theta.combns[, grp := .GRP, by = V2]
    data.table::setkey(current.theta.combns, grp)
    split.thetas <- current.theta.combns[, list(list(.SD)), by = grp]$V1
    thetaF.list <- lapply(split.thetas, function(x) x[, 1])
    all.thetas <- rev(store.all.thetas[[h]])[-length(store.all.thetas[[h]])] # All thetaE values, decreasing, not including the final value, thetaE=0.
    all.thetas <- all.thetas[all.thetas>=minthetaE]

    for(i in 1:length(all.thetas)){ # For each thetaE,
      thetaFs.current <- thetaF.list[[i]] # thetaF.list is a list of i vectors. The i'th vector in the list contains all possible values of thetaF for the i'th thetaE (stored as all.thetas[i])
      rows <- nrow(thetaFs.current)
      # Begin bisection method:
      a <- 1
      b <- nrow(thetaFs.current)
      d <- ceiling((b-a)/2)
      while((b-a)>1){
        output <- findDesignOCs(n=sc.subset$n[h], r=sc.subset$r[h], C=sc.subset$C[h], thetaF=as.numeric(thetaFs.current[d]), thetaE=all.thetas[i], mat=mat.list[[h]],
                                   power=power, alpha=alpha, minstop=sc.subset$eff.minstop[h], coeffs=coeffs, coeffs.p0=coeffs.p0, p0=p0, p1=p1)

        if(output[4] <= alpha) { # type I error decreases as index (and thetaF) increases. Jump backwards if type I error is smaller than alpha, o/w forwards.
          b <- d
        } else {
          a <- d
        }
        d <- a + floor((b-a)/2)
      }

      # Take care of "edge case" where feasible design exists in the first row:
      if(a==1) { # (and also by necessity, b==2)
        output <-  findDesignOCs(n=sc.subset$n[h], r=sc.subset$r[h], C=sc.subset$C[h], thetaF=as.numeric(thetaFs.current[1]), thetaE=all.thetas[i], mat=mat.list[[h]],
                                    power=power, alpha=alpha, minstop=sc.subset$eff.minstop[h], coeffs=coeffs, coeffs.p0=coeffs.p0, p0=p0, p1=p1)
        if(output[4] <= alpha) b <- 1 # Should we start at the first row or the second?
      }

      # We can now proceed moving sequentially from index==b (or cause a break if we wish).
      first.result <-  findDesignOCs(n=sc.subset$n[h], r=sc.subset$r[h], C=sc.subset$C[h], thetaF=as.numeric(thetaFs.current[b]), thetaE=all.thetas[i],
                                        mat=mat.list[[h]], coeffs.p0=coeffs.p0, coeffs=coeffs, power=power, alpha=alpha, minstop=sc.subset$eff.minstop[h], p0=p0, p1=p1)
          if((first.result[4]!=0 | first.result[4]!=2) ) {
        pwr <- first.result[5]
        while(pwr >= power & b <= rows) # Keep going until power drops below 1-beta, i.e. no more feasible designs, or we reach the end of the data frame.
        {
          h.results[[k]] <-  findDesignOCs(n=sc.subset$n[h], r=sc.subset$r[h], C=sc.subset$C[h], thetaF=as.numeric(thetaFs.current[b]), thetaE=all.thetas[i],
                                              mat=mat.list[[h]], coeffs.p0=coeffs.p0, coeffs=coeffs, power=power, alpha=alpha, minstop=sc.subset$eff.minstop[h], p0=p0, p1=p1)
          pwr <- h.results[[k]][5] ####
          k <- k+1
          b <- b+1
        }
      } else { # if first.result[4]==0 or 2, i.e., there are no feasible designs for this thetaE (and hence no remaining feasible designs for this n/r combn), break and move on to next n/r combn:
        break
      }

      # If there is only one possible combination of thetaF and thetaE, this is the result:
      if(rows==1){
        h.results[[k]] <- output
        k <- k+1
      }
    } # end of "i" loop

    if(progressBar==TRUE) setTxtProgressBar(pb, h)

    h.results.df <- do.call(rbind, h.results)

    if(!is.null(h.results.df)){
      # Remove all "skipped" results:
      colnames(h.results.df) <- c("n", "r", "C", "alpha", "power", "EssH0", "Ess", "thetaF", "thetaE", "eff.n", "eff.minstop")
      h.results.df <- h.results.df[!is.na(h.results.df[, "Ess"]),]
      if(nrow(h.results.df)>0){
        # Remove dominated and duplicated designs:
        discard <- rep(NA, nrow(h.results.df))
        for(i in 1:nrow(h.results.df)){
          discard[i] <- sum(h.results.df[i, "EssH0"] > h.results.df[, "EssH0"] & h.results.df[i, "Ess"] > h.results.df[, "Ess"] & h.results.df[i, "n"] >= h.results.df[, "n"])
        }
        h.results.df <- h.results.df[discard==0,, drop=FALSE]
        h.results.list[[h]] <- h.results.df
      }
    }
  } # End of "h" loop


  full.results <- do.call(rbind, h.results.list)
  if(length(full.results)>0) {
    if(return.only.admissible==TRUE){
      # Discard all "inferior" designs:
      discard <- rep(NA, nrow(full.results))
      for(i in 1:nrow(full.results)){
        discard[i] <- sum(full.results[i, "EssH0"] > full.results[, "EssH0"] & full.results[i, "Ess"] > full.results[, "Ess"] & full.results[i, "n"] >= full.results[, "n"])
        }
      full.results <- full.results[discard==0,,drop=FALSE]
    }
    # Remove duplicates:
    duplicates <- duplicated(full.results[, c("n", "Ess", "EssH0"), drop=FALSE])
    final.results <- full.results[!duplicates,,drop=FALSE]
    stage <- final.results[,"n"]/final.results[,"C"]
    final.results <- cbind(final.results, stage)
  } else {
    final.results <- NULL
    if(nmin!=nmax){
      print("There are no feasible designs for this combination of design parameters")
      return(final.results)
    }
  }
  return(final.results)
}

