#' Find two-arm trial designs that use stochastic curtailment
#'
#' This function finds admissible design realisations for two-arm binary outcome trials, using stochastic curtailment.
#' The output can be used as the sole argument in the function 'drawDiagram', which will return the stopping boundaries for the
#' admissible design of your choice. Monitoring frequency can set in terms of block size.
#' @param nmin.arm Smallest maximum sample size *per arm*. Should be a multiple of block size.
#' @param nmax.arm Largest maximum sample size *per arm*. Should be a multiple of block size.
#' @param block.size Block size.
#' @param minstop Minimum permitted sample size at the first interim analysis.
#' @param pc Anticipated response rate on the control arm.
#' @param pt Anticipated response rate on the treatment arm.
#' @param alpha Significance level
#' @param power Required power (1-beta).
#' @param maxthetaF Maximum value of lower CP threshold theta_F_max.
#' @param minthetaE Minimum value of upper threshold theta_E_min.
#' @param bounds choose what final rejection boundaries should be searched over: Those of A'Hern ("ahern"), Wald ("wald") or no constraints (NA). Defaults to "wald".
#' @param max.combns Provide a maximum number of ordered pairs (theta_F, theta_E). Defaults to 1e6.
#' @param fixed.r Choose what final rejection boundaries should be searched over. Useful for reproducing a particular design realisation. Defaults to NULL.
#' @param rm.dominated.designs Logical. If TRUE, dominated designs will be
#' removed from final output. Defaults to TRUE.
#' @param exact.thetaF Provide an exact value for lower threshold theta_F. Useful for reproducing a particular design realisation. Defaults to NULL.
#' @param exact.thetaE Provide an exact value for upper threshold theta_E. Useful for reproducing a particular design realisation. Defaults to NULL.
#' @param fast.method Logical. If FALSE, design search is conducted over all combinations of
#' (theta_F, theta_E). If TRUE, a much faster, though less thorough, design search is undertaken.
#' Defaults to FALSE.
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @return Output is a list of two dataframes. The first, $input, is a one-row data frame that contains all the arguments used in the call.
#' The second, $all.des, contains the operating characteristics of all admissible designs found.
#' @examples
#' \donttest{
#' des <- twoarmDesign(nmin.arm=20,
#' nmax.arm=24,
#' block.size=8,
#' pc=0.1,
#' pt=0.4,
#' alpha=0.1,
#' power=0.8,
#' maxthetaF=0.4,
#' minthetaE=0.7,
#' max.combns=1e4)
#' }
#' @export
twoarmDesign <- function(nmin.arm,
                            nmax.arm,
                            block.size,
                            minstop=NULL,
                            pc,
                            pt,
                            alpha,
                            power,
                            maxthetaF=NULL,
                            minthetaE=0.7,
                            bounds="ahern",
                            fixed.r=NULL,
                            max.combns=1e6,
                            rm.dominated.designs=TRUE,
                            exact.thetaF=NULL,
                            exact.thetaE=NULL,
                            fast.method=FALSE)
{
  Bsize <- block.size

  if(Bsize%%2!=0) stop("Block size must be an even number")

  if((2*nmin.arm)%%Bsize!=0) stop("2*nmin.arm must be a multiple of block size")
  if((2*nmax.arm)%%Bsize!=0) stop("2*nmax.arm must be a multiple of block size")

  if(is.null(minstop)) minstop <- block.size

  if(minstop>=2*nmin.arm) stop("earliest stopping point (minstop) must be smaller than all potential maximum sample sizes, i.e. smaller than 2*nmin.arm.")

  nposs <- seq(from=nmin.arm, to=nmax.arm, by=Bsize/2)

  qc <- 1-pc
  qt <- 1-pt

  prob.vec <- findProbVec(Bsize=Bsize, pt=pt, qt=qt, pc=pc, qc=qc)
  prob.vec.p0 <- findProbVec(Bsize=Bsize, pt=pc, qt=qc, pc=pc, qc=qc)

  prob.vec.minstop <- findProbVec(Bsize=minstop, pt=pt, qt=qt, pc=pc, qc=qc)
  prob.vec.p0.minstop <- findProbVec(Bsize=minstop, pt=pc, qt=qc, pc=pc, qc=qc)

  possible.tps <- seq(Bsize, 2*nmin.arm, by=Bsize) # all possible TPs, ignoring minstop
  eff.minstop <- possible.tps[which.max(possible.tps>=minstop)] # minstop must be multiple of block size

  pat.cols.list <- lapply(nposs, function(x) seq(from=2*x, to=eff.minstop, by=-Bsize)[-1])

  names(pat.cols.list) <- nposs

  if(is.null(maxthetaF)){
    maxthetaF <- pt
  }

  r.list <- list()
  for(i in 1:length(nposs))
  {
    r.list[[i]] <- 0:(nposs[i]-2) # r values: 0 to nposs[i]-2
  }

  ns <- NULL
  for(i in 1:length(nposs)){
    ns <- c(ns, rep(nposs[i], length(r.list[[i]])))
  }

  sc.subset <- data.frame(n=ns, r=unlist(r.list))

  if(!is.null(bounds)){
    # Incorporate A'Hern's bounds:
    if(bounds=="ahern")  {
      #sc.subset <- sc.subset[sc.subset$r >= pc*sc.subset$n & sc.subset$r <= pt*sc.subset$n, ] # One-arm case
      sc.subset <- sc.subset[sc.subset$r >= 1 & sc.subset$r <= pt*sc.subset$n, ] # Try this for two-arm case -- interval [1, pt*Narm]
    }

    if(bounds=="wald"){
      # Even better to incorporate Wald's bounds:
      denom <- log(pt/pc) - log((1-pt)/(1-pc))
      accept.null <- log((1-power)/(1-alpha)) / denom  + nposs * log((1-pc)/(1-pt))/denom
      accept.null <- floor(accept.null)

      reject.null <- log((power)/alpha) / denom  + nposs * log((1-pc)/(1-pt))/denom
      reject.null <- ceiling(reject.null)

      r.wald <- NULL
      ns.wald <- NULL
      for(i in 1:length(nposs)){
        r.wald <- c(r.wald, accept.null[i]:reject.null[i])
        ns.wald <- c(ns.wald, rep(nposs[i], length(accept.null[i]:reject.null[i])))
      }
      sc.subset <- data.frame(n=ns.wald, r=r.wald)
      sc.subset <- sc.subset[sc.subset$n - sc.subset$r >=2, ]
    }
  }

  # In case you want to specify values for r:
  if(!is.null(fixed.r))  {
    sc.subset <- sc.subset[sc.subset$r %in% fixed.r,]
  }

  ###### Find thetas for each possible {r, N} combn:
  mat.list <- vector("list", nrow(sc.subset))
  for(i in 1:nrow(sc.subset)){
    mat.list[[i]] <- findBlock2armUncurtailedMatrix(n=sc.subset[i,"n"], r=sc.subset[i,"r"], Bsize=Bsize, pat.cols=pat.cols.list[[paste(sc.subset$n[i])]], prob.vec=prob.vec)
  }

  store.all.thetas <- lapply(mat.list, function(x) {sort(unique(c(x))[unique(c(x)) <= 1])})


  ##### To cut down on computation, try cutting down the number of thetas used:
  ##### max.combns:=max. number of (thetaF, thetaE) combinations.
  ##### n.thetas*(n.thetas-1)/2 = n.combns, so if n.thetas > sqrt(2*max.combns), take out every other value, excluding 0 and 1.
  ##### Note: further below, more combns are removed if constraints on maxthetaF and minthetaE are specified.
  # check ####
  if(max.combns!=Inf){
    maxthetas <- sqrt(2*max.combns)
    for(i in 1:nrow(sc.subset))
    {
      while(length(store.all.thetas[[i]]) > maxthetas)
      {
        every.other.element <- rep(c(FALSE, TRUE), 0.5*(length(store.all.thetas[[i]])-2))
        store.all.thetas[[i]] <-  store.all.thetas[[i]][c(TRUE, every.other.element, TRUE)]
      }
    }
  }

  if(!is.null(exact.thetaF) & !is.null(exact.thetaE)){ # if exact thetas are given (to speed up a result check):
    for(i in 1:length(store.all.thetas)){
      keep <- abs(store.all.thetas[[i]]-exact.thetaF)<1e-3 | abs(store.all.thetas[[i]]-exact.thetaE)<1e-3
      store.all.thetas[[i]] <- store.all.thetas[[i]][keep]
    }
  }

  h.results.list <- vector("list", nrow(sc.subset)) #

  pb <- txtProgressBar(min = 0, max = nrow(sc.subset), style = 3)

  # Now, find the designs, looping over each possible {r, N} combination, and within each {r, N} combination, loop over all combns of {thetaF, thetaE}:
  for(h in 1:nrow(sc.subset)){
    #k <- 1
    blank.mat <- matrix(NA, nrow=nrow(mat.list[[h]]), ncol=ncol(mat.list[[h]]))
    rownames(blank.mat) <- 0:(nrow(blank.mat)-1)
    zero.mat <- matrix(0, nrow=nrow(mat.list[[h]]), ncol=ncol(mat.list[[h]]))
    rownames(zero.mat) <- rownames(blank.mat)
    pat.cols.single <- pat.cols.list[[paste(sc.subset$n[h])]]

    if(fast.method==TRUE){
      h.results <- fastSearch(thetas=store.all.thetas[[h]],
                              maxthetaF=maxthetaF,
                              minthetaE=minthetaE,
                              sc.h=sc.subset[h,],
                              Bsize=Bsize,
                              mat.h=mat.list[[h]],
                              blank.mat=blank.mat,
                              zero.mat=zero.mat,
                              power=power,
                              alpha=alpha,
                              pat.cols.single=pat.cols.single,
                              prob.vec=prob.vec,
                              prob.vec.p0=prob.vec.p0,
                              prob.vec.minstop=prob.vec.minstop,
                              prob.vec.p0.minstop=prob.vec.p0.minstop,
                              minstop=eff.minstop)
    }else{
      h.results <- slowSearch(thetas=store.all.thetas[[h]],
                              maxthetaF=maxthetaF,
                              minthetaE=minthetaE,
                              sc.h=sc.subset[h,],
                              Bsize=Bsize,
                              mat.h=mat.list[[h]],
                              blank.mat=blank.mat,
                              zero.mat=zero.mat,
                              power=power,
                              alpha=alpha,
                              pat.cols.single=pat.cols.single,
                              prob.vec=prob.vec,
                              prob.vec.p0=prob.vec.p0,
                              prob.vec.minstop=prob.vec.minstop,
                              prob.vec.p0.minstop=prob.vec.p0.minstop,
                              minstop=eff.minstop)
    }


    setTxtProgressBar(pb, h)
    h.results.df <- do.call(rbind, h.results)

    if(!is.null(h.results.df)){
      # Remove all "skipped" results:
      colnames(h.results.df) <- c("n", "r", "block", "alpha", "power", "EssH0", "Ess", "thetaF", "thetaE", "eff.n")
      h.results.df <- h.results.df[!is.na(h.results.df[, "Ess"]),]
      if(nrow(h.results.df)>0){
        # Remove dominated designs:
        if(rm.dominated.designs==TRUE){
          discard <- rep(NA, nrow(h.results.df))
          for(i in 1:nrow(h.results.df)){
            discard[i] <- sum(h.results.df[i, "EssH0"] > h.results.df[, "EssH0"] & h.results.df[i, "Ess"] > h.results.df[, "Ess"] & h.results.df[i, "n"] >= h.results.df[, "n"])
            }
          h.results.df <- h.results.df[discard==0,, drop=FALSE]
        }
        # Remove duplicated designs:
        if(is.matrix(h.results.df)){ # i.e. if there is more than one design (if not, h.results.df is a vector)
        duplicates <- duplicated(h.results.df[, c("n", "Ess", "EssH0"), drop=FALSE])
        h.results.df <- h.results.df[!duplicates,, drop=FALSE]
        }
        h.results.list[[h]] <- h.results.df
      }
    }
  } # End of "h" loop


  full.results <- do.call(rbind, h.results.list)
  #if(length(full.results)==0) stop("There are no feasible designs for this combination of design parameters" , call. = FALSE)
  if(length(full.results)>0){
    # Discard all "inferior" designs:
    discard <- rep(NA, nrow(full.results))
    for(i in 1:nrow(full.results)){
      discard[i] <- sum(full.results[i, "EssH0"] > full.results[, "EssH0"] & full.results[i, "Ess"] > full.results[, "Ess"] & full.results[i, "n"] >= full.results[, "n"])
      #print(i)
    }
    subset.results <- full.results[discard==0,,drop=FALSE]


    # Remove duplicates:
    duplicates <- duplicated(subset.results[, c("n", "EssH0", "Ess"), drop=FALSE])
    all.des <- subset.results[!duplicates,,drop=FALSE]
    all.des$stage <- all.des[,"eff.n"]/all.des[,"block"]
    all.des$eff.minstop <- eff.minstop
    names(all.des)[names(all.des)=="n"] <- "n.arm"
    names(all.des)[names(all.des)=="eff.n"] <- "n"
    input <- data.frame(nmin.arm=nmin.arm,
                        nmax.arm=nmax.arm,
                        block=block.size,
                        minstop=minstop,
                        pc=pc,
                        pt=pt,
                        alpha=alpha,
                        power=power,
                        maxthetaF=maxthetaF,
                        minthetaE=minthetaE,
                        bounds=bounds,
                        max.combns=max.combns,
                        fast.method=fast.method)
    final.output <- list(input=input,
                         all.des=all.des)
    class(final.output) <- append(class(final.output), "curtailment_twoarm")
    return(final.output)
  }
}


# rmDominatedDesigns <- function(df, essh0="EssH0", essh1="Ess", n="n"){
#     discard <- rep(NA, nrow(df))
#     if("tbl_df" %in% class(df)){
#         essh0.vec <- df[[essh0]]
#         essh1.vec <- df[[essh1]]
#         n.vec <- df[[n]]
#       for(i in 1:nrow(df)){
#         discard[i] <-  any(essh0.vec[i] > essh0.vec & essh1.vec[i] > essh1.vec & n.vec[i] >= n.vec)
#       }
#     } else {
#         essh0.vec <- df[, essh0]
#         essh1.vec <- df[, essh1]
#         n.vec <- df[, n]
#         for(i in 1:nrow(df)){
#           discard[i] <- any(essh0.vec[i] > essh0.vec & essh1.vec[i] > essh1.vec & n.vec[i] >= n.vec)
#         }
#       }
#     newdf <- df[discard==FALSE,,drop=FALSE]
#     newdf
# }
