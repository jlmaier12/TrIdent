#' Counts zero values to the left and right of prophage-like borders
#'
#' Checks to see at which point the number of consecutive zero values to the
#' left and right of the prophage-like pattern match borders equals the
#' noReadCov parameter
#'
#' @param startOrEnd
#' searching the start (left side) or end (right side) of the prophage-like
#' pattern match
#' @param viralSubsetZoom
#' viralSubset dataframe subsetted to 50,000 bp outside the pattern match
#' borders
#' @param startOrEndPosRow
#' The row index of the start or end position of the prophage-like pattern match
#' @param noReadCov
#' How many bp of no read coverage are encountered before specialized
#' transduction searching stops? Default is 500.
#' @return List
#' @keywords internal
zeroCountSearch <-
    function(startOrEnd,
            viralSubsetZoom,
            startOrEndPosRow,
            noReadCov) {
        X <- ifelse(startOrEnd == "start", 1, -1)
        zeroCount <- 0
        repeat {
            if ((startOrEndPosRow - X) <= 1) {
                break
            }
            if (startOrEndPosRow - X >= nrow(viralSubsetZoom)) {
                break
            }
            zeroCount <-
                ifelse((viralSubsetZoom[(startOrEndPosRow - X), 4] == 0),
                    (zeroCount + 1), 0
                )
            if (zeroCount == noReadCov / 100) {
                break
            }
            X <- ifelse(startOrEnd == "start", X + 1, X - 1)
        }
        return(list(zeroCount, abs(X)))
    }
