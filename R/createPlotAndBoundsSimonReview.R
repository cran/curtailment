createPlotAndBoundsSimonReview <- function(des, xmax=NULL, ymax=NULL){
  m <- Sm <- decision <- analysis <- NULL
  des <- as.data.frame(t(des))
  coords <- expand.grid(0:des$n, 1:des$n)
  diag.df <- data.frame(Sm=as.numeric(coords[,1]),
                        m=as.numeric(coords[,2]),
                        decision=rep("Continue", nrow(coords)))
  diag.df$decision <- as.character(diag.df$decision)
  diag.df$decision[coords[,1]>coords[,2]] <- NA

  fails.sm <- c(0:des$r1, (des$r1+1):des$r)
  fails.m <- c(rep(des$n1, length(0:des$r1)),
               rep(des$n, length((des$r1+1):des$r))
               )
  tp.fail <- data.frame(Sm=fails.sm,
                        m=fails.m)
  tp.success <- data.frame(Sm=(des$r+1):des$n,
                           m=rep(des$n, length((des$r+1):des$n))
                           )
  # If stopping for benefit:
  if("e1" %in% names(des)){
    tp.success.s1 <- data.frame(Sm=(des$e1+1):des$n1,
                                m=rep(des$n1, length((des$e1+1):des$n1))
                                )
    tp.success <- rbind(tp.success, tp.success.s1)
  }

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

  diag.df.subset <- diag.df[!is.na(diag.df$decision),]
  diag.df.subset$analysis <- "No"
  diag.df.subset$analysis[diag.df.subset$m %in% tp.fail$m] <- "Yes"

  # plot.title <- "Stopping boundaries"
  # sub.text1 <- paste("Max no. of analyses: 2. Max(N): ", des$n, ". ESS(p", sep="")
  # sub.text2 <- paste("): ", round(des$EssH0, 1), ". ESS(p", sep="")
  # sub.text3 <- paste("):", round(des$Ess, 1), sep="")
  # plot.subtitle2 <- bquote(.(sub.text1)[0]*.(sub.text2)[1]*.(sub.text3))
  diagram <- pkgcond::suppress_warnings(ggplot2::ggplot(data=diag.df.subset, mapping = aes(x=m, y=Sm, fill=decision, alpha=analysis))+
                                          scale_alpha_discrete(range=c(0.5, 1)),
                                        "Using alpha for a discrete variable is not advised")
  diagram <- diagram +
    geom_tile(color="white")+
     labs(fill="Decision",
          alpha="Analysis",
          x="Number of participants",
          y="Number of responses"#,
    #      title=plot.title,
    #      subtitle = plot.subtitle2
    )+
    coord_cartesian(expand = 0)+
    theme_minimal()

  xbreaks <- c(des$n1, des$n)

  if(!is.null(xmax)){
    diagram <- diagram +
      expand_limits(x=xmax)
    xbreaks <- c(xbreaks, xmax)
  }

  if(!is.null(ymax)){
    diagram <- diagram +
      expand_limits(y=ymax)
  }

  diagram <- diagram +
    scale_x_continuous(breaks=xbreaks)+
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))

  print(diagram)

  stop.bounds <- data.frame(m=c(des$n1, des$n),
                           success=c(Inf, des$r+1),
                           fail=c(des$r1, des$r))
  return(list(diagram=diagram,
              bounds.mat=stop.bounds))
}
