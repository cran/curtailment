#### Part 5: from "SC_bisection_oct2020.R" ####

############## THIS CODE FINDS ALPHA AND POWER FOR ALL COMBNS OF n1/n2/n/r1/r, given p/p0, UNDER NSC

#' findSCDesigns
#'
#' This function finds admissible design realisations for single-arm binary outcome trials, using stochastic curtailment.
#' This function differs from singlearmDesign in that it includes a Simon-style interim analysis after some n1 participants.
#' The output is a data frame of admissible design realisations.
#' @param nmin Minimum permitted sample size.
#' @param nmax Maximum permitted sample size.
#' @param p0 Probability for which to control the type-I error-rate
#' @param p1 Probability for which to control the power
#' @param alpha Significance level
#' @param power Required power (1-beta).
#' @param maxthetaF Maximum value of lower CP threshold theta_F_max. Defaults to p.
#' @param minthetaE Minimum value of upper threshold theta_E_min. Defaults to p.
#' @param bounds choose what final rejection boundaries should be searched over: Those of A'Hern ("ahern"), Wald ("wald") or no constraints (NA). Defaults to "wald".
#' @param max.combns Provide a maximum number of ordered pairs (theta_F, theta_E). Defaults to 1e6.
#' @param maxthetas Provide a maximum number of CP values used to create ordered pairs (theta_F, theta_E). Can be used instead of max.combns. Defaults to NA.
#' @param fixed.r1 Choose what interim rejection boundaries should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param fixed.r Choose what final rejection boundaries should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param fixed.n1 Choose what interim sample size values n1 should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param exact.thetaF Provide an exact value for lower threshold theta_F. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param exact.thetaE Provide an exact value for upper threshold theta_E. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param progressBar Logical. If TRUE, shows progress bar. Defaults to FALSE.
#' @return Output is a list of two dataframes. The first, $input, is a one-row
#' data frame that contains important arguments used in the call. The second,
#' $all.des,contains the operating characteristics of all admissible designs found.
#' @export
#' @examples \donttest{findSCdesigns(nmin = 20, nmax = 20, p0 = 0.1, p1 = 0.4, power = 0.8, alpha = 0.1, max.combns=1e2)}
findSCdesigns <- function(nmin,
                 nmax,
                 p0,
                 p1,
                 alpha,
                 power,
                 minthetaE=p1,
                 maxthetaF=p1,
                 bounds="wald",
                 fixed.r1=NA,
                 fixed.r=NA,
                 fixed.n1=NA,
                 max.combns=1e6,
                 maxthetas=NA,
                 exact.thetaF=NA,
                 exact.thetaE=NA,
                 progressBar=FALSE){
  grp <- V2 <- NULL
  #sc.all <- findN1N2R1R2_df(nmin=nmin, nmax=nmax, e1=FALSE)
  sc.all <- findSimonN1N2R1R2(nmin=nmin, nmax=nmax, e1=FALSE)

  if(is.numeric(bounds)){
    sc.subset <- sc.all[sc.all$r %in% bounds,]
  }else{
    if(bounds=="ahern"){
      ### AHERN'S BOUNDS: ###
      sc.subset <- sc.all[sc.all$r>=p0*sc.all$n & sc.all$r<=p1*sc.all$n, ]
    }else{
      if(bounds=="wald"){
        # Even better to incorporate Wald's bounds:
        nposs <- nmin:nmax
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
        #sc.subset <- data.frame(n=ns.wald, r=r.wald)
        keep.index <- NULL
        for(i in 1:nrow(sc.all)){
          keep.index[i] <- sum(sc.all$n[i]==ns.wald & sc.all$r[i]==r.wald) > 0
        }
        sc.subset <- sc.all[keep.index,]
      }
    }
  }


  if(!is.na(fixed.r1)){
    sc.subset <- sc.subset[sc.subset$r1 %in% fixed.r1,]
  }
  if(!is.na(fixed.r)){
    sc.subset <- sc.subset[sc.subset$r %in% fixed.r,]
  }
  if(!is.na(fixed.n1)){
    sc.subset <- sc.subset[sc.subset$n1 %in% fixed.n1,]
  }

  # For each design, find the potential values of theta:
  n.design.rows <- nrow(sc.subset)
  store.all.thetas <- vector("list", n.design.rows)

  for(i in 1:n.design.rows)
  {
    store.all.thetas[[i]] <- findThetas(n1=sc.subset[i,1], n2=sc.subset[i,2], n=sc.subset[i,3], r1=sc.subset[i,4], r2=sc.subset[i,5], p0=p0, p=p1)
  }
  # Takes seconds




  # Much of the required computing is independent of thetaF/1, and therefore only needs to be computed once per design.
  # The function
  #
  #       outsideTheta
  #
  # undertakes this computation, and stores all the necessary objects for obtaining power, type I error, ESS, etc. for later
  # use in the function
  #
  #       insideTheta
  #


  outside.theta.data <- vector("list", nrow(sc.subset))
  all.theta.combns <- vector("list", nrow(sc.subset))

if(!is.na(max.combns) & is.na(maxthetas)){ # Only use max.combns if maxthetas is not specified.
    maxthetas <- sqrt(2*max.combns)
}

  if(is.na(exact.thetaF)){ # if exact theta values have NOT been supplied:
  ##### To cut down on computation, try cutting down the number of thetas used.
  ##### If a design has > maxthetas theta values (CP values), cut this down:
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
    outside.theta.data[[i]] <- outsideTheta(n1=sc.subset[i,1], n2=sc.subset[i,2], r1=sc.subset[i,4], r2=sc.subset[i,5], p0=p0, p=p1)
    all.theta.combns[[i]] <- t(combn(store.all.thetas[[i]], 2))
    all.theta.combns[[i]] <- all.theta.combns[[i]][all.theta.combns[[i]][,2] >= minthetaE & all.theta.combns[[i]][,1] <= maxthetaF, ] # Reduce number of combinations
    if(length(all.theta.combns[[i]])==2) all.theta.combns[[i]] <- matrix(c(0, 0.999, 0, 1), nrow=2, byrow=T) # To avoid a crash. See (*) below
    all.theta.combns[[i]] <- all.theta.combns[[i]][order(all.theta.combns[[i]][,2], decreasing=T),]  ### XXX
  }
  # (*) If the only thetas for a design are 0 and 1, then some code below will crash. This
  # could be elegantly resolved by adding an if or ifelse statement, but probably at some
  # computational expense. Will deal with these exceptions here, by adding a "false" value
  # for theta.
  }else{ # if exact theta values have been supplied:
    for(i in 1:length(store.all.thetas)){
      keep <- abs(store.all.thetas[[i]]-exact.thetaF)<1e-2 | abs(store.all.thetas[[i]]-exact.thetaE)<1e-2
      store.all.thetas[[i]] <- store.all.thetas[[i]][keep]
    }
    for(i in 1:nrow(sc.subset)){
      outside.theta.data[[i]] <- outsideTheta(n1=sc.subset[i,1], n2=sc.subset[i,2], r1=sc.subset[i,4], r2=sc.subset[i,5], p0=p0, p=p1)
      all.theta.combns[[i]] <- t(combn(store.all.thetas[[i]], 2))
      all.theta.combns[[i]] <- all.theta.combns[[i]][order(all.theta.combns[[i]][,2], decreasing=T),]  ### XXX
    }
  }







  # NOTE: thetaE is fixed while thetaF increases. As thetaF increases, power decreases.

  if(progressBar==TRUE) pb <- txtProgressBar(min = 0, max = length(outside.theta.data), style = 3)
  h.results.list <- vector("list", length(outside.theta.data)) #
  y <- NULL

  for(h in 1:length(outside.theta.data)) # For every possible combn of n1/n2/n/r1/r:
  {
    k <- 1
    h.results <- vector("list", nrow(all.theta.combns[[h]]))
    current.theta.combns <- all.theta.combns[[h]] # Take just the thetaF/1 combns for that design.
    # current.theta.combns <- current.theta.combns[current.theta.combns[,2]>0.7, ] # remove combns where thetaE < 0.7 (as not justifiable)

    # Don't need a matrix of all thetaF and thetaE combns -- potentially quicker to have thetaF as a vector (already defined), and the list can be a list of thetaE vectors:
    current.theta.combns <- data.table::data.table(current.theta.combns)
    current.theta.combns[, grp := .GRP, by = V2]
    data.table::setkey(current.theta.combns, grp)
    split.thetas <- current.theta.combns[, list(list(.SD)), by = grp]$V1
    thetaF.list <- lapply(split.thetas, function(x) x[, 1])
    all.thetas <- rev(store.all.thetas[[h]])[-length(store.all.thetas[[h]])] # All thetaE values, decreasing, not including the final value, thetaE=0.
    all.thetas <- all.thetas[all.thetas>=minthetaE]

    # Of course we need the other characteristics and so on of the design, which are stored in the list "outside.theta.data":
    y <- outside.theta.data[[h]]
    eff.n <- y$n

    for(i in 1:length(all.thetas)) # For each thetaE,
    {
      thetaFs.current <- thetaF.list[[i]] # thetaF.list is a list of i vectors. The i'th vector in the list contains all possible values of thetaF for the i'th thetaE (stored as all.thetas[i])
      rows <- nrow(thetaFs.current)

      # Begin bisection method: Use the bisection method to find where to "start":
      a <- 1
      b <- nrow(thetaFs.current)
      d <- ceiling((b-a)/2)

      while((b-a)>1)
      {
        output <- insideTheta(n1=y$n1, n2=y$n2, n=y$n, r1=y$r1, r2=y$r2, thetaF=as.numeric(thetaFs.current[d]), thetaE=all.thetas[i], cp.sm=y$cp.sm, cp.m=y$cp.m, p0=p0, p=p1,
                              mat=y$mat, coeffs.p0=y$coeffs.p0, coeffs=y$coeffs, power=power, alpha=alpha)
        if(output[6] <= alpha) { # type I error decreases as index (and thetaF) increases. Jump backwards if type I error is smaller than alpha, o/w forwards.
          b <- d } else {
            a <- d
          }
        #  print(paste("a is ", a, "... b is ", b, "... d is ", d, sep=""))
        d <- a + floor((b-a)/2)
      }

      # Take care of "edge case" where feasible design exists in the first row:
      if(a==1) # (and also by necessity, b==2)
      {
        output <- insideTheta(n1=y$n1, n2=y$n2, n=y$n, r1=y$r1, r2=y$r2, thetaF=as.numeric(thetaFs.current[1]), thetaE=all.thetas[i], cp.sm=y$cp.sm, cp.m=y$cp.m, p0=p0, p=p1,
                              mat=y$mat, coeffs.p0=y$coeffs.p0, coeffs=y$coeffs, power=power, alpha=alpha)
        if(output[6] <= alpha) b <- 1 # Should we start at the first row or the second?
      }

      # We can now proceed, moving sequentially from index==b (or cause a break if we wish):

      first.result <-  insideTheta(n1=y$n1, n2=y$n2, n=y$n, r1=y$r1, r2=y$r2, thetaF=as.numeric(thetaFs.current[b]), thetaE=all.thetas[i], cp.sm=y$cp.sm, cp.m=y$cp.m, p0=p0, p=p1,
                                   mat=y$mat, coeffs.p0=y$coeffs.p0, coeffs=y$coeffs, power=power, alpha=alpha)
      # If the type I error equals zero, then the trial is guaranteed to fail at m=1. Further, this will also be the case for all subsequent thetaF values
      # (as thetaF is increasing). Further still, there are no feasible designs for subsequent thetaE values, and so we can skip to the next n/r combination
      # If type I error doesn't equal zero, then feasible designs exist and they should be recorded. Note: type I error equals 2 means success guaranteed
      # at m=1, and the conclusion is the same.
      if(first.result[6]!=0 | first.result[6]!=2)
      {
        pwr <- first.result[7]
        while(pwr >= power & b <= rows) # Keep going until power drops below 1-beta, i.e. no more feasible designs, or we reach the end of the data frame.
        {
          h.results[[k]] <-  insideTheta(n1=y$n1, n2=y$n2, n=y$n, r1=y$r1, r2=y$r2, thetaF=as.numeric(thetaFs.current[b]), thetaE=all.thetas[i], cp.sm=y$cp.sm, cp.m=y$cp.m, p0=p0, p=p1,
                                         mat=y$mat, coeffs.p0=y$coeffs.p0, coeffs=y$coeffs, power=power, alpha=alpha)
          pwr <- h.results[[k]][7] #### Added 3rd July -- should have been there all the time ¬_¬
          k <- k+1
          b <- b+1
        }
      } else { # if first.result[3]==0 or 2, i.e., there are no feasible designs for this thetaE (and hence no remaining feasible designs for this n/r combn), break:
        break
      }

    } # end of "i" loop

    if(progressBar==TRUE) setTxtProgressBar(pb, h)

    h.results.df <- do.call(rbind.data.frame, h.results)

    # Before moving on to the next set of values, avoid creating a huge object by cutting out all dominated and duplicated designs for this combn of n1/n2/n/r1/r:
    if(nrow(h.results.df)>0)
    {
      names(h.results.df) <- c("n1", "n2", "n", "r1", "r2", "alpha", "power", "EssH0", "Ess", "thetaF", "thetaE", "eff.n")

      # Remove all "skipped" results:
      h.results.df <- h.results.df[!is.na(h.results.df$Ess),]
      # Remove dominated and duplicated designs:
      discard <- rep(NA, nrow(h.results.df))
      for(i in 1:nrow(h.results.df))
      {
        discard[i] <- sum(h.results.df$EssH0[i] > h.results.df$EssH0 & h.results.df$Ess[i] > h.results.df$Ess & h.results.df$n[i] >= h.results.df$n)
      }
      h.results.df <- h.results.df[discard==0,]

      duplicates <- duplicated(h.results.df[, c("n", "Ess", "EssH0")])
      h.results.df <- h.results.df[!duplicates,]

      h.results.list[[h]] <- h.results.df
    }

  }

  # All designs are combined in the dataframe "results".
  results <- do.call(rbind.data.frame, h.results.list)

  # Remove all "skipped" results:
  results <- results[!is.na(results$Ess),]

  # Remove dominated and duplicated designs:
  discard <- rep(NA, nrow(results))
  for(i in 1:nrow(results))
  {
    discard[i] <- sum(results$EssH0[i] > results$EssH0 & results$Ess[i] > results$Ess & results$n[i] >= results$n)
  }
  subset.results <- results[discard==0,]

  # Remove duplicates:
  duplicates <- duplicated(subset.results[, c("n", "Ess", "EssH0")])
  subset.results <- subset.results[!duplicates,]
  sc.input <- data.frame(nmin=nmin, nmax=nmax, p0=p0, p1=p1, alpha=alpha, power=power)
  sc.output <- list(input=sc.input,
                       all.des=subset.results)
  return(sc.output)
}
