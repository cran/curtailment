#' drawDiagramGeneric
#'
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @param n Maximum sample size
#' @param go Two-column matrix detailing stopping boundaries for a go decision
#' or efficacy stopping, where the first column contains the number of responses and
#' the second column contains the corresponding number of participants so far.
#' @param nogo Two-column matrix detailing stopping boundaries for a no go decision
#' or futility stopping, where the first column contains the number of responses and
#' the second column contains the corresponding number of participants so far.
#' @param xmax Optional. Maximum value for x axis.
#' @param ymax Optional. Maximum value for y axis.
#' @param ylab Optional. Label for y axis.
#' @return An object of class ``ggplot``, showing a visualisation of the
#' maximum sample size and inputted stopping boundaries.
#' @export
#' @examples
#'  go <- cbind(6:8, rep(8,3))
#'  nogo <- cbind(0:5, rep(8,6))
#'  drawDiagramGeneric(n=8, go=go, nogo=nogo, ylab="Number of successes")
drawDiagramGeneric <- function(n,
                               go,
                               nogo,
                               xmax=NULL,
                               ymax=NULL,
                               ylab=NULL){
  m <- Sm <- decision <- NULL
  if(length(go)==2) go <- matrix(go, ncol=2)
  if(length(nogo)==2) nogo <- matrix(nogo, ncol=2)

  tp.success <- as.data.frame(x=go)
  names(tp.success) <- c("Sm", "m")

  tp.fail <- as.data.frame(x=nogo)
  names(tp.fail) <- c("Sm", "m")

  analysis <- sort(unique(c(tp.success$m, tp.fail$m)))

  coords <- expand.grid(0:n, 1:n)
  diag.df <- data.frame(Sm=as.numeric(coords[,1]),
                        m=as.numeric(coords[,2]),
                        decision=rep("Continue", nrow(coords)))
  diag.df$decision <- as.character(diag.df$decision)
  diag.df$decision[coords[,1]>coords[,2]] <- NA

  success.index <- apply(diag.df, 1, function(y) any(as.numeric(y[1])==tp.success$Sm & as.numeric(y[2])==tp.success$m))
  diag.df$decision[success.index] <- "Go decision"
  fail.index <- apply(diag.df, 1, function(y) any(as.numeric(y[1])==tp.fail$Sm & as.numeric(y[2])==tp.fail$m))
  diag.df$decision[fail.index] <- "No go decision"

  for(i in 1:nrow(tp.fail)){
    not.poss.fail.index <- diag.df$Sm==tp.fail$Sm[i] & diag.df$m>tp.fail$m[i]
    diag.df$decision[not.poss.fail.index] <- NA
  }

  for(i in 1:nrow(tp.success)){
    not.poss.pass.index <- diag.df$m-diag.df$Sm==tp.success$m[i]-tp.success$Sm[i] & diag.df$m>tp.success$m[i]
    diag.df$decision[not.poss.pass.index] <- NA
  }


  # Add shading:
  diag.df.subset <- diag.df[!is.na(diag.df$decision),]
  diag.df.subset$analysis <- "No"
  stop.index <- diag.df.subset$m %in% analysis
  diag.df.subset$analysis[stop.index] <- "Yes"

  plot.title <- "Stopping boundaries"
  plot.subtitle2 <- paste("Max no. of analyses: ", length(analysis), ". Max(N): ", n, ".", sep="")

  diagram <- pkgcond::suppress_warnings(ggplot2::ggplot(data=diag.df.subset, mapping = ggplot2::aes(x=m, y=Sm, fill=decision, alpha=analysis))+
                                          ggplot2::scale_alpha_discrete(range=c(0.5, 1)),
                                        "Using alpha for a discrete variable is not advised")
  diagram <- diagram +
    ggplot2::geom_tile(color="white")+
    ggplot2::labs(fill="Decision",
                  alpha="Analysis",
                  x="Number of participants",
                  y="Number of responses",
                  title=plot.title,
                  subtitle = plot.subtitle2)+
    ggplot2::coord_cartesian(expand = 0)+
    ggplot2::theme_minimal()

  xbreaks <- unique(c(analysis, n))

  if(!is.null(xmax)){
    diagram <- diagram +
      ggplot2::expand_limits(x=xmax)
    xbreaks <- c(xbreaks, xmax)
  }

  if(!is.null(ymax)){
    diagram <- diagram +
      expand_limits(y=ymax)
  }

  if(!is.null(ylab)){
    diagram <- diagram +
      ggplot2::labs(y=ylab)
  }

  diagram <- diagram +
    ggplot2::scale_x_continuous(breaks=xbreaks)+
    ggplot2::scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
  return(diagram)
}
