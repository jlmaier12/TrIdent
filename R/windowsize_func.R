#' Change the read coverage window size
#'
#' Re-averages window sizes of read coverage averages. Start with 100bp windows always. Cannot make window size less than 100bp.
#'
#' @param read_dataset A read coverage dataset that has been cleaned and reformatted by the readcovdf_formatter function
#' @param windowsize The number of base pairs to average coverage values over. Options are 100, 500, 1000, or 2000 only!
#' @keywords internal
windowsize_func <- function(read_dataset, windowsize){
  coverage <- vector()
  X <- 0
  Y <- windowsize/100
  repeat{
    coverage <- c(coverage, mean(read_dataset[c(X:Y),2]))
    X <- X+(windowsize/100)
    Y <- Y+(windowsize/100)
    if (Y > nrow(read_dataset)) break
  }
  position <- seq(windowsize, length(coverage)*windowsize, windowsize)
  ref_name <- rep(read_dataset[1,1], length(position))
  newdataset <- cbind.data.frame(ref_name, coverage, position) %>% as.data.frame()
  newdataset[is.nan.data.frame(newdataset)] <- 0
  return(newdataset)
}
