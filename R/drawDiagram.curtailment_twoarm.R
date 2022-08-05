#' @export
drawDiagram.curtailment_twoarm <- function(findDesign.output, print.row=NULL, xmax=NULL, ymax=NULL){
  des <- findDesign.output$all.des
  row.names(des) <- 1:nrow(des)
  if(!is.null(print.row)){
    des <- des[print.row, , drop=FALSE]
  }
  print(des)
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
        plot.and.bounds <- createPlotAndBounds2arm(des=des, des.input=des.input, rownum=rownum, xmax=xmax, ymax=ymax)
      }
    } # end of while
  }else{
    print("Returning diagram and bounds for single design.", quote = F)
    plot.and.bounds <- createPlotAndBounds2arm(des=des, des.input=des.input, rownum=1, xmax=xmax, ymax=ymax)
    return(plot.and.bounds)
  }
} # end of drawDiagram()
