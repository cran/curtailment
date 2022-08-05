
#### Plotting ####
#### Plotting omni-admissible designs (Fig 6 in manuscript):

design.plot <- function(results.df, loss.df, scenario, design)
{
  qs <- q0 <- q1 <- NULL
  index <- which(results.df$design==design)
  all.loss.subset <- t(loss.df[index,])
  designs <- results.df[index,]
  designs.vec <- as.character(as.vector(apply(all.loss.subset, 1, which.min)))
  #code <- c("1"=design.details[1], "2"=design.details[2], "3"=design.details[3])
  #simon.designs.scen1.vec <- code[simon.designs.scen1.vec]
  designs.used <- sort(as.numeric(unique(designs.vec)))

  if(design=="simon" | design=="nsc")  design.details <- apply(results.df[index[designs.used],], 1, function(x) paste("{", x["r1"], "/", x["n1"], ", ", x["r"], "/", x["n"], "}", sep=""))
  if(design=="simon_e1") design.details <- apply(results.df[index[designs.used],], 1, function(x) paste("{(", x["r1"], " ", x["e1"], ")/", x["n1"], ", ", x["r"], "/", x["n"], "}", sep=""))
  if(design=="sc") design.details <- apply(results.df[index[designs.used],], 1, function(x) paste("{", x["r1"], "/", x["n1"], ", ", x["r"], "/", x["n"], ", ", format(round(as.numeric(x["thetaF"]),2), nsmall=2), "/", format(round(as.numeric(x["thetaE"]),2), nsmall=2), "}", sep=""))
  if(design=="mstage")  design.details <- apply(results.df[index[designs.used],], 1, function(x) paste("{", x["r"], "/", x["n"], ", ", format(round(as.numeric(x["thetaF"]),2), nsmall=2), "/", format(round(as.numeric(x["thetaE"]),2), nsmall=2), "}", sep=""))

  designs.df <- data.frame(design=designs.vec, q0=qs[,"w0"], q1=qs[, "w1"])
  designs.df$design <- factor(designs.df$design, levels = as.character(sort(as.numeric(levels(designs.df$design)))))
  design.type <- switch(design, "simon"="Simon", "simon_e1"="Mander Thompson", "nsc"="NSC", "sc"="SC", "mstage"="m-stage")

  output <-  ggplot2::ggplot(designs.df,  ggplot2::aes(x=q1, y=q0)) +
     ggplot2::geom_raster(ggplot2::aes(fill = design)) +
     ggplot2::scale_fill_discrete(name="Design",labels=as.vector(design.details)) +
     ggplot2::ggtitle(paste(design.type, ", scenario ", scenario, sep="")) +
     ggplot2::xlab(expression(w[1])) +
     ggplot2::ylab(expression(w[0])) +
     ggplot2::theme_bw()+
     ggplot2::theme(#axis.text.x=element_text(),
      legend.text= ggplot2::element_text(size=7),
      plot.title =  ggplot2::element_text(size =  ggplot2::rel(1.4)),
      axis.ticks= ggplot2::element_blank(),
      axis.line= ggplot2::element_blank(),
      axis.title= ggplot2::element_text(size=rel(1.5)),
      axis.title.y= ggplot2::element_text(angle = 0, vjust=0.5),
      axis.text= ggplot2::element_text(size=rel(1.5)),
      legend.key.size =  ggplot2::unit(0.4, "cm"),
      legend.position = c(0.8,0.75),
      panel.border= ggplot2::element_blank(),
      panel.grid.major= ggplot2::element_line(color='#eeeeee'))
  return(output)
}

plot.all <- function(results, loss, scen){
  design.types <- c("simon", "simon_e1", "nsc", "sc", "mstage")
  plots.list <- vector("list", 5)
  for(i in 1:5){
    plots.list[[i]] <- design.plot(results.df=results, loss.df=loss, scenario=scen, design=design.types[i])
  }
  plots.list
}




# Expected loss:

plotExpdLoss <- function(loss.df, design.type, scenario, design.loss.range){
  title <- paste("Expected loss: ", design.type, ", Scenario ", scenario, sep="")
  ggplot2::ggplot(loss.df, ggplot2::aes_string("w1", "w0")) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::xlab(expression(w[1])) +
    ggplot2::ylab(expression(w[0])) +
    ggplot2::geom_raster(ggplot2::aes_string(fill="loss")) +
    ggplot2::scale_fill_gradient(low="red", high="darkblue", limits=design.loss.range) +
    ggplot2::theme(#axis.text.x=element_text(),
      axis.text=ggplot2::element_text(size=15),
      axis.title.x = ggplot2::element_text(size=15),
      axis.title.y = ggplot2::element_text(size=15, angle = 0, vjust = 0.5),
      legend.text=ggplot2::element_text(size=12),
      plot.title = ggplot2::element_text(size = rel(1.5)),
      axis.ticks=ggplot2::element_blank(),
      axis.line=ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(0.75, "cm"),
      legend.position = c(0.9,0.65),
      legend.title = ggplot2::element_text(size=15),
      panel.border=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_line(color='#eeeeee'))
}

plotAllExpdLoss <- function(loss.list, scen, loss.range){
  design.types <- c("Simon", "MT", "NSC", "SC", "m-stage")
  plots.list <- vector("list", 5)
  for(i in 1:5){
    plots.list[[i]] <- plotExpdLoss(loss.df=loss.list[[i]],
                                    design.type = design.types[i],
                                    scenario=scen,
                                    design.loss.range=loss.range)
  }
  plots.list
}




##### Plotting quantiles of ESS #####
cohortFindMatrix <- function(n, r, Csize, p0, p1){
  q1 <- 1-p1
  mat <- matrix(3, ncol=n, nrow=min(r+Csize, n)+1) #  nrow=r+Csize+1 unless number of stages equals 1.
  rownames(mat) <- 0:(nrow(mat)-1)
  mat[(r+2):nrow(mat),] <- 1
  mat[1:(r+1),n] <- 0

  pat.cols <- seq(n, 1, by=-Csize)[-1]

  for(i in (r+1):1){
    for(j in pat.cols){  # Only look every C patients (no need to look at final col)
      if((r+1-(i-1) > n-j)) { # IF success is not possible [Total responses needed (r+1) - responses so far (i-1)] > [no. of patients remaining (n-j)], THEN
        mat[i,j] <- 0
      }else{
        if(i-1<=j){ # Condition: Sm<=m
          newcp <- sum(choose(Csize,0:Csize)*p1^(0:Csize)*q1^(Csize:0)*mat[i:(i+Csize),j+Csize])
          # NOTE: No rounding of CP up to 1 or down to 0: This is the basic matrix to which many pairs of (thetaF, thetaE) will be applied in the function cohortFindDesign(). Lines below can be deleted.
          # if(newcp > thetaE) mat[i,j] <- 1
          # if(newcp < thetaF) mat[i,j] <- 0
          # if(newcp <= thetaE & newcp >= thetaF) mat[i,j] <- newcp
          mat[i,j] <- newcp
        }
      }
    }
  }


  for(i in 3:nrow(mat)) {
    mat[i, 1:(i-2)] <- NA
  }
  return(mat)
}


cohortFindQuantileSS <- function(n, r, Csize, thetaF, thetaE, p0, p1, coarseness=0.01, lower=0.1, upper=0.9, return.only.quantiles=TRUE){

  mat <- cohortFindMatrix(n = n, r=r, Csize=Csize, p0=p0, p1=p1)

  q1 <- 1-p1
  q0 <- 1-p0

  # Note: Only need the first r+Csize+1 rows -- no other rows can be reached.
  coeffs <- matrix(NA, ncol=n, nrow=r+Csize+1)
  coeffs.p0 <- coeffs

  for(i in 1:(r+Csize+1)){
    coeffs[i, ] <- p1^(i-1) * q1^((2-i):(2-i+n-1))
    coeffs.p0[i, ] <- p0^(i-1) * q0^((2-i):(2-i+n-1))
  }


  ##### LIST of matrix of coefficients: Probability of a path leading to each point:
  p.vec <- seq(0, 1, by=coarseness)
  coeffs.list <- vector("list", length=length(p.vec))

  for(k in 1:length(coeffs.list)){
    coeffs.list[[k]] <- matrix(NA, ncol=n, nrow=r+Csize+1)
    for(i in 1:(r+Csize+1)){
      coeffs.list[[k]][i, ] <- p.vec[k]^(i-1) * (1-p.vec[k])^((2-i):(2-i+n-1))
    }
  }


  pat.cols <- seq(n, 1, by=-Csize)[-1]

  # Amend CP matrix, rounding up to 1 when CP>theta_1 and rounding down to 0 when CP<thetaF:

  for(i in (r+1):1){
    for(j in pat.cols){  # Only look every Csize patients (no need to look at final col)
      if(i-1<=j){ # Condition: Sm<=m
        newcp <- sum(choose(Csize,0:Csize)*p1^(0:Csize)*q1^(Csize:0)*mat[i:(i+Csize),j+Csize])
        if(newcp > thetaE) mat[i,j] <- 1
        if(newcp < thetaF) mat[i,j] <- 0
        if(newcp <= thetaE & newcp >= thetaF) mat[i,j] <- newcp
      }
    }
  }



  ############# Number of paths to each point:

  pascal.list <- list(c(1,1))
  if(Csize>1){
    for(i in 2:Csize){
      pascal.list[[i]] <- c(0, pascal.list[[i-1]]) + c(pascal.list[[i-1]], 0)
    }
  }

  for(i in (Csize+1):n){
    if(i %% Csize == 1 | Csize==1){
      column <- as.numeric(mat[!is.na(mat[,i-1]), i-1])
      newnew <- pascal.list[[i-1]]
      CPzero.or.one <- which(column==0 | column==1)
      newnew[CPzero.or.one] <- 0
      pascal.list[[i]] <- c(0, newnew) + c(newnew, 0)
    } else {
      pascal.list[[i]] <- c(0, pascal.list[[i-1]]) + c(pascal.list[[i-1]], 0)
    }
  }

  for(i in 1:length(pascal.list)){
    pascal.list[[i]] <- c(pascal.list[[i]], rep(0, length(pascal.list)+1-length(pascal.list[[i]])))
  }


  pascal.t <- t(do.call(rbind, args = pascal.list))
  pascal.t <- pascal.t[1:(r+Csize+1), ]

  # Multiply the two matrices (A and p^b * q^c). This gives the probability of reaching each point:
  # Note: Only need the first r+Csize+1 rows -- no other rows can be reached.
  pascal.t <- pascal.t[1:(r+Csize+1), ]

  final.probs.mat <- pascal.t*coeffs
  final.probs.mat.p0 <- pascal.t*coeffs.p0

  final.probs.mat.list <- vector("list", length = length(p.vec))
  for(i in 1:length(p.vec)){
    final.probs.mat.list[[i]] <- pascal.t*coeffs.list[[i]]
  }


  ##### Now find the terminal points: m must be a multiple of Csize, and CP must be 1 or 0:
  ##### SUCCESS:
  ## Only loop over rows that contain CP=1:
  rows.with.cp1 <- which(apply(mat, 1, function(x) {any(x==1, na.rm = T)}))

  m.success <- NULL
  Sm.success <- NULL
  prob.success <- NULL
  prob.success.p0 <- NULL
  prob.success.list <- vector("list", length(p.vec))

  interims <- seq(Csize, n, by=Csize)

  for(i in rows.with.cp1){
    for(j in interims[interims >= i-1]){ # condition ensures Sm >= m
      if(mat[i,j]==1){
        m.success <- c(m.success, j)
        Sm.success <- c(Sm.success, i-1)
        prob.success <- c(prob.success, final.probs.mat[i, j])
        prob.success.p0 <- c(prob.success.p0, final.probs.mat.p0[i, j])
        for(k in 1:length(p.vec)){
          prob.success.list[[k]] <- c(prob.success.list[[k]], final.probs.mat.list[[k]][i,j])
        }
      }
    }
  }

  prob.success.mat <- do.call(cbind, prob.success.list)
  success <- data.frame(Sm=Sm.success, m=m.success, prob=prob.success, prob.p0=prob.success.p0, success=rep("Success", length(m.success)), prob.success.mat)
  names(success)[6:ncol(success)]


  pwr <- sum(success[,"prob"])
  typeI <- sum(success[,"prob.p0"])


  ##### FAILURE:
  first.zero.in.row <- apply(mat[1:(r+1),], 1, function(x) {which.max(x[interims]==0)})
  m.fail <- Csize * as.numeric(first.zero.in.row)
  Sm.fail <- as.numeric(names(first.zero.in.row))

  fail.deets <- data.frame(Sm=Sm.fail, m=m.fail)
  fail.deets$prob <- apply(fail.deets, 1, function(x) {final.probs.mat[x["Sm"]+1, x["m"]]})
  fail.deets$prob.p0 <- apply(fail.deets, 1, function(x) {final.probs.mat.p0[x["Sm"]+1, x["m"]]})

  prob.fail.list <- vector("list", length(p.vec))
  for(i in 1:length(p.vec)){
    prob.fail.list[[i]] <- apply(fail.deets, 1, function(x) {final.probs.mat.list[[i]][x["Sm"]+1, x["m"]]})
  }
  prob.fail.mat <- do.call(cbind, prob.fail.list)

  fail.deets$success <- rep("Fail", nrow(fail.deets))

  fail.deets <- cbind(fail.deets, prob.fail.mat)

  names(fail.deets) <- names(success)
  output <- rbind(fail.deets, success)

  ###### ADDED

  ordered.outp <- output[order(output$m),]
  cum.probs.mat <- apply(ordered.outp[,c(6:ncol(ordered.outp))], 2, cumsum)

  index.median <- apply(cum.probs.mat, 2, function(x) which.max(x>=0.5))
  actual.median <- ordered.outp$m[index.median]

  index.lower.percentile <- apply(cum.probs.mat, 2, function(x) which.max(x>=lower))
  actual.lower.percentile <- ordered.outp$m[index.lower.percentile]

  index.upper.percentile <- apply(cum.probs.mat, 2, function(x) which.max(x>=upper))
  actual.upper.percentile <- ordered.outp$m[index.upper.percentile]

  ###### END OF ADDED


  sample.size.expd <- sum(output$m*output$prob)
  sample.size.expd.p0 <- sum(output$m*output$prob.p0)


  ############ EFFECTIVE N. THE "EFFECTIVE N" OF A STUDY IS THE "REAL" MAXIMUM SAMPLE SIZE
  ######## The point is where every Sm for a given m equals zero or one is necessarily where a trial stops

  cp.colsums <- apply(mat, 2, function(x) { sum(x==0, na.rm=TRUE)+sum(x==1, na.rm=TRUE)} ) # Sum the CP values that equal zero or one in each column
  possible.cps <- apply(mat, 2, function(x) {sum(!is.na(x))})

  effective.n <- min(which(cp.colsums==possible.cps))
  if(return.only.quantiles==FALSE){
    to.return <- list(row=c(n, r, Csize, typeI, pwr, sample.size.expd.p0, sample.size.expd, thetaF, thetaE, effective.n),
                      quantiles=data.frame(p.vec=p.vec, median=actual.median, lower=actual.lower.percentile, upper=actual.upper.percentile, design="m-stage", cohort.size=Csize, stages=n/Csize))
  } else {
    to.return <- data.frame(p.vec=p.vec, median=actual.median, lower=actual.lower.percentile, upper=actual.upper.percentile, design="m-stage")
    to.return$cohort.size <- Csize
    to.return$stages <- n/to.return$cohort.size
  }
  return(to.return)
}


simonQuantiles <- function(r1, n1, n2, coarseness=0.1, lower=0.1, upper=0.9){
  p.vec <- seq(from=0, to=1, by=coarseness)
  pet <- pbinom(r1, n1, p.vec)
  actual.median <- rep(NA, length(p.vec))
  actual.median[pet > 0.5] <- n1
  actual.median[pet == 0.5] <- n1 + 0.5*n2
  actual.median[pet < 0.5] <- n1 + n2
  actual.10.percentile <- rep(NA, length(p.vec))
  actual.10.percentile[pet > lower] <- n1
  actual.10.percentile[pet == lower] <- n1 + 0.5*n2
  actual.10.percentile[pet < lower] <- n1 + n2
  actual.90.percentile <- rep(NA, length(p.vec))
  actual.90.percentile[pet > upper] <- n1
  actual.90.percentile[pet == upper] <- n1 + 0.5*n2
  actual.90.percentile[pet < upper] <- n1 + n2
  output <- data.frame(p.vec=p.vec, median=actual.median, lower=actual.10.percentile, upper=actual.90.percentile,
                       design="simon", cohort.size=0.5*(n1+n2), stages=2)
  return(output)
}


plotMedianSSDesign <- function(dataf, p0, p1, factor, title){
  required.columns <- c("p.vec", "median", "lower", "upper")
  if(!all(required.columns %in% names(dataf))){
    stop(cat("Data frame must have the following columns: ", required.columns))
  }
  dataf[, factor] <- factor(dataf[, factor])

  ggplot2::ggplot(dataf, ggplot2::aes_string(x="p.vec", y="median", ymin="lower", ymax="upper"), ggplot2::aes_string(group="factor")) +
    ggplot2::geom_line(ggplot2::aes_string(color="factor"), size=1) +
    ggplot2::geom_ribbon(ggplot2::aes_string(fill="factor"), alpha = 0.3) +
    ggplot2::labs(title=title,
         x ="p", y = "Sample size", fill="Design\n(block size B)", color="Design\n(block size B)")+
    ggplot2::geom_vline(xintercept = p0, size=0.5, color="grey60", linetype="longdash") +
    ggplot2::geom_vline(xintercept = p1, size=0.5, color="grey60", linetype="longdash") +
    ggplot2::theme_bw()+
    ggplot2::theme(legend.text=ggplot2::element_text(size=rel(1.5)),
          legend.title=ggplot2::element_text(size=rel(1.5)),
          plot.title = ggplot2::element_text(size = rel(1.5)),
          axis.ticks=ggplot2::element_blank(),
          axis.line=ggplot2::element_blank(),
          axis.text.x=ggplot2::element_text(angle=0),
          axis.title=ggplot2::element_text(size=rel(1.5)),
          axis.text=ggplot2::element_text(size=rel(1.5))
    )
}
