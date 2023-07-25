#' NA remover
#'
#' Removes NAs from dataframe. Function found on Stack Exchange.
#'
#' @seealso
#' \url{https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame/18143097#18143097}
#'
#' @param x dataset with potential NAs
#' @keywords internal
is.nan.data.frame <- function(x) {
  do.call(cbind, lapply(x, is.nan))}

