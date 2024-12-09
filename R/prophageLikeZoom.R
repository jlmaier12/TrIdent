#' Prophage-like pattern zoom
#'
#' 'Zoom-in' on (aka subset) desired region surrounding block pattern.
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param classificationPatterns The pattern match information associated with
#'   each contig classified as Prophage-like, sloping, or HighCovNoPattern
#' @param i The index for the contig currently being assessed
#' @param zoom The number of rows outside the start and stop positions of the
#'   block pattern to zoom-in on
#' @param windowSize The window size used to re-average read coverage pileups
#' @return Dataframe
#' @keywords internal
prophageLikeZoom <-
    function(viralSubset,
            classificationPatterns,
            i,
            zoom,
            windowSize) {
        startRowIdx <- classificationPatterns[[i]][[5]] * windowSize / 100
        endRowIdx <- classificationPatterns[[i]][[6]] * windowSize / 100
        zoomStartRowIdx <-
            ifelse((startRowIdx - zoom) < 1, 1, (startRowIdx - zoom))
        zoomEndRowIdx <- ifelse((endRowIdx + zoom > nrow(viralSubset)),
            nrow(viralSubset),
            (endRowIdx + zoom)
        )
        zoomedViralSubset <- viralSubset[c(zoomStartRowIdx:zoomEndRowIdx), ]
        return(zoomedViralSubset)
    }
