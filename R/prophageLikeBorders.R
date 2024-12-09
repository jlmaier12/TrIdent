#' Prophage-like border finder
#'
#' Find borders of Prophage-like patterns with more specificity than
#' pattern-matching using 100 bp window pileups and sliding standard deviation
#' technique.
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param classificationPatterns The pattern match information associated with
#'   each contig classified as Prophage-like, Sloping, or HighCovNoPattern
#' @param i The index for the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage pileups
#' @return List
#' @importFrom roll roll_sd
#' @keywords internal
prophageLikeBorders <-
    function(viralSubset,
            classificationPatterns,
            i,
            windowSize) {
        startRowIdx <- classificationPatterns[[i]][[5]] * windowSize / 100
        endRowIdx <- classificationPatterns[[i]][[6]] * windowSize / 100
        if (classificationPatterns[[i]][[5]] == 1) {
            leftBorderRowIdx <- 1
        } else {
            searchStartRowIdx <-
                ifelse((startRowIdx - 25) < 1, 1, (startRowIdx - 25))
            leftSubset <-
                viralSubset[c(searchStartRowIdx:(startRowIdx + 25)), 2]
            rollingSd <- roll_sd(leftSubset, width = 2)
            rollingSdNoNA <- rollingSd[-is.na(rollingSd)]
            leftBorderRowIdx <-
                which(rollingSd ==
                        sort(rollingSdNoNA,
                             TRUE)[1])[length(which(rollingSd ==
                                                      sort(rollingSdNoNA,
                                                           TRUE)[1]))]
            leftBorderRowIdx <- searchStartRowIdx + leftBorderRowIdx
        }
        if (endRowIdx + 1 == (nrow(viralSubset))) {
            rightBorderRowIdx <- nrow(viralSubset)
        } else {
            searchEndRowIdx <- ifelse((endRowIdx + 25 > nrow(viralSubset)),
                nrow(viralSubset), (endRowIdx + 25)
            )
            rightSubset <-
                viralSubset[c((endRowIdx - 25):searchEndRowIdx), 2]
            rollingSd <- roll_sd(rightSubset, width = 2)
            rollingSdNoNA <- rollingSd[-is.na(rollingSd)]
            rightBorderRowIdx <-
                which(rollingSd == sort(rollingSdNoNA, TRUE)[1])[1]
            rightBorderRowIdx <- ifelse((endRowIdx - 25) +
                rightBorderRowIdx > nrow(viralSubset),
            nrow(viralSubset),
            (endRowIdx - 25) + rightBorderRowIdx
            )
        }
        leftBorderBpPos <- viralSubset[leftBorderRowIdx, 3]
        rightBorderBpPos <- viralSubset[rightBorderRowIdx, 3]
        margins <- list(leftBorderBpPos, rightBorderBpPos)
        return(margins)
    }
