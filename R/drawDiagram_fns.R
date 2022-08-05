
#' drawDiagram
#'
#' This function produces both a data frame and a diagram of stopping boundaries.
#' The function takes a single argument: the output from the function singlearmDesign.
#' If the supplied argument contains more than one admissible designs, the user is offered a choice of which design to use.
#' @param  findDesign.output Output from either the function singlearmDesign or find2stageDesigns
#' @param  print.row Choose a row number to directly obtain a plot and stopping boundaries for a particular design realisation. Default is NULL.
#' @param  xmax,ymax Choose values for the upper limits of the x- and y-axes respectively. Helpful for comparing two design realisations. Default is NULL.
#' @export
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @return The output is a list of two elements. The first, $diagram, is a ggplot2 object showing how the trial should proceed: when to to undertake an interim analysis, that is, when to check if a stopping boundary has been reached (solid colours) and what decision to make at each possible point (continue / go decision / no go decision). The second list element, $bounds.mat, is a data frame containing three columns: the number of participants at which to undertake an interim analysis (m), and the number of responses at which the trial should stop for a go decision (success) or a no go decision (fail).
#' @examples output <- singlearmDesign(nmin = 30,
#'  nmax = 30,
#'  C = 5,
#'  p0 = 0.1,
#'  p1 = 0.4,
#'  power = 0.8,
#'  alpha = 0.05)
#'  dig <- drawDiagram(output, print.row = 2)
drawDiagram <- function(findDesign.output, print.row=NULL, xmax=NULL, ymax=NULL){
  UseMethod("drawDiagram")
}

#' @export
drawDiagram.curtailment_single <- function(findDesign.output, print.row=NULL, xmax=NULL, ymax=NULL){
    des <- findDesign.output$all.des
  row.names(des) <- 1:nrow(des)
  if(!is.null(print.row)){
    des <- des[print.row, , drop=FALSE]
  }
  ##print(des)
  des.input <- findDesign.output$input
  if(nrow(des)>1){
    rownum <- 1
    while(is.numeric(rownum)){
      rownum <- readline("Input a row number to choose a design and see the trial design diagram. Press 'q' to quit: ")
      if(rownum=="q"){
        if(exists("plot.and.bounds")){
          return(plot.and.bounds)
        }else{
          print("No designs selected, nothing to return", q=F)
          return()
        }
      }else{
        rownum <- as.numeric(rownum)
        plot.and.bounds <- createPlotAndBounds(des=des, des.input=des.input, rownum=rownum, xmax=xmax, ymax=ymax)
      }
    } # end of while
  }else{
    #print("Returning diagram and bounds for single design.", quote = F)
    plot.and.bounds <- createPlotAndBounds(des=des, des.input=des.input, rownum=1, xmax=xmax, ymax=ymax)
    return(plot.and.bounds)
  }
} # end of drawDiagram()

#' @export
drawDiagram.curtailment_simon <- function(findDesign.output, print.row=NULL, xmax=NULL, ymax=NULL){
  des <- findDesign.output$all.des
  row.names(des) <- 1:nrow(des)
  if(!is.null(print.row)){
    des <- des[print.row, , drop=FALSE]
  }
  #print(des)
  des.input <- findDesign.output$input
  if(nrow(des)>1){
    rownum <- 1
    while(is.numeric(rownum)){
      rownum <- readline("Input a row number to choose a design and see the trial design diagram. Press 'q' to quit: ")
      if(rownum=="q"){
        if(exists("plot.and.bounds")){
          return(plot.and.bounds)
        }else{
          print("No designs selected, nothing to return", q=F)
          return()
        }
      }else{
        rownum <- as.numeric(rownum)
        plot.and.bounds <- createPlotAndBoundsSimon(des=des, des.input=des.input, rownum=rownum, xmax=xmax, ymax=ymax)
      }
    } # end of while
  }else{
    #print("Returning diagram and bounds for single design.", quote = F)
    plot.and.bounds <- createPlotAndBoundsSimon(des=des, des.input=des.input, rownum=1, xmax=xmax, ymax=ymax)
    return(plot.and.bounds)
  }
}
