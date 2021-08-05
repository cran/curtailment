createPlotAndBounds <- function(des, des.input, rownum, xmax, ymax){
  m <- Sm <- decision <- analysis <- NULL
  rownum <- as.numeric(rownum)
  des <- as.data.frame(t(des[rownum, ]))
  coef.list <- findCoeffs(n=des$n, p0 = des.input$p0, p1 = des.input$p1)
  matt <- findCPmatrix(n=des$n, r=des$r, Csize=des$C, p0=des.input$p0, p1=des.input$p1) # uncurtailed matrix
  findo <- findDesignOCs(n=des$n, r = des$r, C=des$C, thetaF = des$thetaF, thetaE = des$thetaE, p1 = des.input$p1, p0 = des.input$p0, alpha = des.input$alpha, power = des.input$power,
                         return.tps = TRUE, coeffs = coef.list$coeffs, coeffs.p0 = coef.list$coeffs.p0, mat=matt)
  tp.success <- findo$tp[findo$tp$success=="Success", ]
  # Order by Sm (though this should already be the case):
  tp.success <- tp.success[order(tp.success$Sm),]
  tp.success.unneeded <- tp.success[duplicated(tp.success$m), ]
  tp.fail <- findo$tp[findo$tp$success=="Fail", ]

  coords <- expand.grid(0:des$n, 1:des$n)
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

  # for(i in 1:nrow(tp.success.unneeded)){
  #   unneeded.success.index <- diag.df$Sm==tp.success.unneeded$Sm[i] & diag.df$m==tp.success.unneeded$m[i]
  #   diag.df$decision[unneeded.success.index] <- NA
  # }

  # Add shading:
  diag.df.subset <- diag.df[!is.na(diag.df$decision),]
  diag.df.subset$analysis <- "No"
  stop.index <- diag.df.subset$m %in% unique(findo$tp$m)
  diag.df.subset$analysis[stop.index] <- "Yes"

  plot.title <- "Stopping boundaries"
  sub.text1 <- paste("Max no. of analyses: ", des$stage, ". Max(N): ", des$n, ". ESS(p", sep="")
  sub.text2 <- paste("): ", round(des$EssH0, 1), ". ESS(p", sep="")
  sub.text3 <- paste("):", round(des$Ess, 1), sep="")
  plot.subtitle2 <- bquote(.(sub.text1)[0]*.(sub.text2)[1]*.(sub.text3))
  diagram <- pkgcond::suppress_warnings(ggplot2::ggplot(data=diag.df.subset, mapping = aes(x=m, y=Sm, fill=decision, alpha=analysis))+
                                          scale_alpha_discrete(range=c(0.5, 1)),
                                        "Using alpha for a discrete variable is not advised")
  diagram <- diagram +
    geom_tile(color="white")+
    labs(fill="Decision",
         alpha="Analysis",
         x="Number of participants",
         y="Number of responses",
         title=plot.title,
         subtitle = plot.subtitle2)+
    coord_cartesian(expand = 0)+
    theme_minimal()

  xbreaks <- seq(from=des$C, to=des$n, by=des$C)

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

  tp.success.unique.m <- tp.success[!duplicated(tp.success$m), ]
  stop.bounds <- data.frame(m=seq(from=des$C, to=des$n, by=des$C),
                            success=Inf,
                            fail=-Inf)
  stop.bounds$success[match(tp.success.unique.m$m, stop.bounds$m)] <- tp.success.unique.m$Sm
  stop.bounds$fail[match(tp.fail$m, stop.bounds$m)] <- tp.fail$Sm
  return(list(diagram=diagram,
              bounds.mat=stop.bounds))
}
