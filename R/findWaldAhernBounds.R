# Compare boundaries to Wald/A'Hern ####
# Find the boundaries:
findWaldAhernBounds <- function(N.range, beta, alpha, p0, p1, round=F){
  findWaldBounds <- function(N.vec, bet, alph, p0., p1., rounding){
    denom <- log(p1./p0.) - log((1-p1.)/(1-p0.))
    accept.null <- log((1-bet)/alph) / denom  +  N.vec * log((1-p0.)/(1-p1.))/denom
    reject.null <- log(bet/(1-alph)) / denom  +  N.vec * log((1-p0.)/(1-p1.))/denom
    if(rounding==TRUE){
      accept.null <- ceiling(accept.null)
      reject.null <- floor(reject.null)
    }
    output <- data.frame(N=N.vec, lower=reject.null, upper=accept.null, design="Wald")
    output
  }
  findAhernBounds <- function(N.vec, p0., p1., rounding){
    ahern.lower <- N.vec*p0.
    ahern.upper <- N.vec*p1.
    if(rounding==TRUE){
      ahern.lower <- floor(ahern.lower)
      ahern.upper <- ceiling(ahern.upper)
    }
    output <- data.frame(N=N.vec, lower=ahern.lower, upper=ahern.upper, design="A'Hern")
    output
  }

  wald.output <- findWaldBounds(N.vec=N.range, bet=beta, alph=alpha, p0.=p0, p1.=p1, rounding=round)
  ahern.output <- findAhernBounds(N.vec=N.range, p0.=p0, p1.=p1, rounding=round)
  output <- rbind(ahern.output, wald.output)
  output
}

# Plot the boundaries along with the m-stage boundaries (two different data frames):
plotBounds <- function(wald.ahern.bounds, mstage.bounds, title="XXX"){
  print(range(wald.ahern.bounds$N))
  plot.output <- ggplot2::ggplot()+
    ggplot2::geom_ribbon(data=wald.ahern.bounds,
                         ggplot2::aes(x=wald.ahern.bounds$N,
                                      y=NULL,
                                      ymin=wald.ahern.bounds$lower,
                                      ymax=wald.ahern.bounds$upper,
                                      fill=wald.ahern.bounds$design),
                         alpha = 0.3)+
    ggplot2::geom_point(data=mstage.bounds,
                        ggplot2::aes(x=mstage.bounds$n,
                                     y=mstage.bounds$r))+
    ggplot2::labs(title=title,
                  x ="Maximum sample size",
                  y = "Final rejection boundary",
                  fill="Design")+
    ggplot2::theme_bw()+
    ggplot2::scale_x_continuous(limits=range(wald.ahern.bounds$N), expand=c(0,1))
  plot.output
}
