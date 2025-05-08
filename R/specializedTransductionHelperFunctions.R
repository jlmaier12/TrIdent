#' Specialized transduction search and plot
#'
#' Search contigs classified as prophage-like for potential specialized
#' transduction and return the plot visualizing the search results.
#'
#' @param VLPpileup A table containing contig names, coverages averaged over 100
#'   bp windows, and contig positions associated with mapping VLP-fraction reads
#'   to whole-community contigs
#' @param classifPatternMatches The pattern match information associated with
#'   each contig classified as prophage-like,  sloping, or HighCovNoPattern
#' @param windowSize The window size used to re-average read coverage pileups
#' @param i The index for the contig currently being assessed
#' @param noReadCov How many bp of no read coverage are encountered before
#'   searching stops? Default is 500.
#' @param specTransLength How many bp of read coverage to look for outside of
#'   prophage borders? Default is 2000.
#' @param classifSumm The summary information associated with each contig
#'   classified as Prophage-like, Sloping, or HighCovNoPattern
#' @param contigName The reference name of the contig currently being assessed
#'   (i.e "NODE_1")
#' @param logScale If TRUE, coverage is plotted in log10. If FALSE, raw coverage
#'   values are plotted. Default is FALSE.
#' @return List containing two objects
#' @keywords internal
specTransductionSearch <- function(contigName,
                                   VLPpileup,
                                   classifPatternMatches,
                                   classifSumm,
                                   windowSize,
                                   i,
                                   noReadCov,
                                   specTransLength,
                                   logScale) {
  specTransSumm <- c(contigName, rep(NA, 5))
  viralSubsetZoom <- prophageLikeZoom(
    VLPpileup[which(VLPpileup[, 1] ==
      contigName), ],
    classifPatternMatches, i, 500, windowSize
  )
  borders <-
    prophageLikeBorders(
      VLPpileup[which(VLPpileup[, 1] == contigName), ],
      classifPatternMatches, i, windowSize
    )
  startPosBp <- borders[[1]]
  startPosRow <- which(viralSubsetZoom[, 3] == startPosBp)
  endPosBp <- borders[[2]]
  endPosRow <- which(viralSubsetZoom[, 3] == endPosBp)
  viralSubsetZoom$logcoverage <- abs(log10(viralSubsetZoom[, 2]))
  viralSubsetZoom[viralSubsetZoom == Inf] <- 0
  startSearch <-
    zeroCountSearch("start", viralSubsetZoom, startPosRow, noReadCov)
  endSearch <-
    zeroCountSearch("end", viralSubsetZoom, endPosRow, noReadCov)
  if (startSearch[[2]] - startSearch[[1]] >= specTransLength / 100) {
    SpecTransLeft <- viralSubsetZoom[startPosRow - (startSearch[[2]] -
      startSearch[[1]]), 3]
    specTransSumm[c(3, 5)] <- c("yes", (startSearch[[2]] -
      startSearch[[1]]) * 100)
  } else {
    SpecTransLeft <- startPosBp
    specTransSumm[c(3, 5)] <- c("no", NA)
  }
  if (endSearch[[2]] - endSearch[[1]] >= specTransLength / 100) {
    specTransRight <- viralSubsetZoom[endPosRow + (endSearch[[2]] -
      endSearch[[1]]), 3]
    specTransSumm[c(4, 6)] <-
      c("yes", (endSearch[[2]] - endSearch[[1]]) * 100)
  } else {
    specTransRight <- endPosBp
    specTransSumm[c(4, 6)] <- c("no", NA)
  }
  specTransSumm[2] <-
    ifelse((
      startSearch[[2]] - startSearch[[1]] >= specTransLength / 100 |
        endSearch[[2]] - endSearch[[1]] >= specTransLength / 100
    ),
    "yes",
    "no"
    )
  plot <- specTransductionPlot(
    viralSubsetZoom,
    startPosBp,
    endPosBp,
    SpecTransLeft,
    specTransRight,
    contigName,
    classifPatternMatches,
    i,
    specTransSumm,
    logScale,
    classifSumm
  )
  return(list(specTransSumm, plot))
}


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
#' @param windowSize The window size used to re-average read coverage pileups
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
          sort(
            rollingSdNoNA,
            TRUE
          )[1])[length(which(rollingSd ==
          sort(
            rollingSdNoNA,
            TRUE
          )[1]))]
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



#' Counts zero values to the left and right of prophage-like borders
#'
#' Checks to see at which point the number of consecutive zero values to the
#' left and right of the prophage-like pattern match borders equals the
#' noReadCov parameter
#'
#' @param startOrEnd searching the start (left side) or end (right side) of the
#'   prophage-like pattern-match
#' @param viralSubsetZoom viralSubset dataframe subsetted to 50,000 bp outside
#'   the pattern match borders
#' @param startOrEndPosRow The row index of the start or end position of the
#'   prophage-like pattern match
#' @param noReadCov How many bp of no read coverage are encountered before
#'   specialized transduction searching stops? Default is 500.
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


#' Specialized transduction plot
#'
#' Plot search results of `specializedTransductionID()`
#'
#' @param viralSubsetZoom contig subset surrounding Prophage-like pattern-match
#' @param startPosBp Left border position
#' @param endPosBp Right border position
#' @param SpecTransLeft End position of spec transduction on left border
#' @param specTransRight End position of spec transduction on right border
#' @param contigName The reference name of the contig currently being assessed
#'   (i.e "NODE_1")
#' @param classifPatternMatches The pattern match information associated with
#'   each contig classified as prophage-like,  sloping, or HighCovNoPattern
#' @param i The index for the contig currently being assessed
#' @param specTransSumm Results for spec transduction search
#' @param logScale If TRUE, coverage is plotted in log10. If FALSE, raw coverage
#'   values are plotted. Default is FALSE.
#' @param classifSumm The summary information associated with each contig
#'   classified as Prophage-like, Sloping, or HighCovNoPattern
#' @return ggplot object
#' @keywords internal
specTransductionPlot <-
  function(viralSubsetZoom,
           startPosBp,
           endPosBp,
           SpecTransLeft,
           specTransRight,
           contigName,
           classifPatternMatches,
           i,
           specTransSumm,
           logScale,
           classifSumm) {
    position <- logcoverage <- NULL
    if (classifSumm[
      which(classifSumm[, 1] == contigName),
      8
    ] == "Elevated") {
      prophageLikeInfo <- "Highly active/abundant prophage-like element"
    } else if (classifSumm[
      which(classifSumm[, 1] == contigName),
      8
    ] == "Depressed") {
      prophageLikeInfo <-
        "Not homogenously present/integrated prophage-like element"
    } else {
      prophageLikeInfo <- NULL
    }
    fill <-
      ifelse(specTransSumm[2] == "yes", "seagreen", "deepskyblue3")
    alphaL <- ifelse(specTransSumm[3] == "yes", 1, 0)
    alphaR <- ifelse(specTransSumm[4] == "yes", 1, 0)
    coverageType <-
      if (logScale) {
        viralSubsetZoom$logcoverage
      } else {
        viralSubsetZoom$coverage
      }
    plot <- (
      ggplot(data = as.data.frame(viralSubsetZoom), aes(
        x = position,
        y = coverageType
      )) +
        geom_area(fill = fill) +
        geom_vline(
          xintercept = c(startPosBp, endPosBp),
          linewidth = 1
        ) +
        geom_vline(
          xintercept = SpecTransLeft,
          color = "red",
          alpha = alphaL,
          linewidth = 1
        ) +
        geom_vline(
          xintercept = specTransRight,
          color = "red",
          alpha = alphaR,
          linewidth = 1
        ) +
        scale_x_continuous(expand = c(0, 0)) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 12),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 15),
          plot.title = element_text(size = 14)
        ) +
        labs(
          title = paste(
            contigName, prophageLikeInfo
          ),
          subtitle = paste(
            "Specialized transduction on left:",
            specTransSumm[3],
            ",",
            "on right:",
            specTransSumm[4]
          ),
          x = "Contig Position (bp)",
          y = paste(
            "VLP-Fraction Read Coverage",
            ifelse(logScale,
              "(Log10)", ""
            )
          )
        )
    )
    return(plot)
  }
