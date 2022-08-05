#function generator
defunct <- function(msg = "This function is depreciated") function(...) return(stop(msg))
# @export
#' SCfn = defunct("SCfn changed name to findSCdesigns")
