#' Summarizes pattern-match information
#'
#' Summarizes the classifications made in the patternMatcher() function into a
#' dataframe.
#'
#' @param bestMatchList Classifications made with patternMatcher function.
#' @return dataframe
#' @keywords internal
contigClassSumm <- function(bestMatchList) {
  if (length(bestMatchList) == 0) {
    stop("NO TRANSDUCTION EVENTS FOUND")
  }
  contigName <- rep(NA, length(bestMatchList))
  classifications <- rep(NA, length(bestMatchList))
  normMatchScore <- rep(NA, length(bestMatchList))
  lapply(seq_along(bestMatchList), function(i) {
    normMatchScore[i] <<- bestMatchList[[i]][[9]]
    contigName[i] <<- bestMatchList[[i]][[8]]
    classifications[i] <<- bestMatchList[[i]][[7]]
  })
  classifSumm <-
    cbind.data.frame(contigName, classifications, normMatchScore)
  return(classifSumm)
}


#' Summarize slopes for sloping classifications
#'
#' Add slope information for sloping classifications to summary table
#'
#' @param classifSumm Classification summary table
#' @param slopingClassifList A list containing pattern match information
#'   associated with all contigs classified as sloping.
#' @param windowSize The window size used to re-average read coverage pileups
#' @return dataframe
#' @keywords internal
slopeSumm <- function(classifSumm, slopingClassifList, windowSize) {
  classifSumm$slope <- rep(NA, nrow(classifSumm))
  if (length(slopingClassifList) == 0) {
    return(classifSumm)
  }
  lapply(seq_along(slopingClassifList), function(i) {
    classifSumm[which(classifSumm[, 1] ==
      slopingClassifList[[i]][[8]]), 10] <<-
      round(slopingClassifList[[i]][[4]] / windowSize, digits = 4)
  })
  return(classifSumm)
}

#' Pattern-match size calculator
#'
#' Calculate the size (bp) of the matching region for Prophage-like and Sloping
#' patterns
#'
#' @param classifSumm Classification summary table
#' @param classifList A list containing pattern match information associated
#'   with all contig classifications
#' @param windowSize The window size used to re-average read coverage pileups
#' @return dataframe
#' @keywords internal
patternMatchSize <- function(classifSumm, classifList, windowSize) {
  message("Determining sizes (bp) of pattern matches")
  classifSumm <- as.data.frame(classifSumm)
  classifSumm$matchSize <- rep(NA, nrow(classifSumm))
  classifSumm$startPosBp <- rep(NA, nrow(classifSumm))
  classifSumm$endPosBp <- rep(NA, nrow(classifSumm))
  lapply(seq_along(classifList), function(i) {
    contigName <- classifList[[i]][[8]]
    startPos <- classifList[[i]][[5]]
    endPos <- classifList[[i]][[6]]
    classifSumm[which(classifSumm[, 1] == contigName), 5] <<-
      (length(c(startPos:endPos)) - 1) * windowSize
    classifSumm[which(classifSumm[, 1] == contigName), 6] <<-
      startPos * windowSize
    classifSumm[which(classifSumm[, 1] == contigName), 7] <<-
      endPos * windowSize
  })
  return(classifSumm)
}

#' VLP-fraction:whole-community read coverage ratio calculator
#'
#' Calculate the VLP-fraction:whole-community read coverage ratio for every
#' contig using the median read coverage values. If the ratio is greater than 2
#' (i.e VLP-fraction read coverage is, on average, at least double the
#' whole-community read coverage), then the contig is classified as
#' HighCovNoPattern
#'
#' @param classifSumm Classification summary table
#' @param WCpileup A table containing contig names, coverages averaged over 100
#'   bp windows, and contig positions associated with mapping whole-community
#'   reads to whole-community contigs
#' @param VLPpileup A table containing contig names, coverages averaged over 100
#'   bp windows, and contig positions associated with mapping VLP-fraction reads
#'   to whole-community contigs
#' @return dataframe
#' @keywords internal
VLPtoWCRatioCalc <- function(classifSumm, WCpileup, VLPpileup) {
  classifSumm$VLPWCRatio <- rep(NA, nrow(classifSumm))
  noneClassifIdxs <- which(classifSumm[, 2] == "NoPattern")
  if (length(noneClassifIdxs) == 0) {
    return(classifSumm)
  }
  lapply(seq_along(noneClassifIdxs), function(p) {
    i <- noneClassifIdxs[[p]]
    contigName <- classifSumm[i, 1]
    viralSubset <- VLPpileup[which(VLPpileup[, 1] == contigName), ]
    viralSubset[NARemover(viralSubset)] <- 0
    microbialSubset <-
      WCpileup[which(WCpileup[, 1] == contigName), ]
    microbialSubset[NARemover(microbialSubset)] <- 0
    VLPtoWCratio <-
      round(median(viralSubset[, 2] / median(microbialSubset[, 2])),
        digits = 4
      )
    classifSumm[i, 2] <<-
      ifelse(VLPtoWCratio > 2, "HighCovNoPattern", "NoPattern")
    classifSumm[i, 4] <<- VLPtoWCratio
  })
  return(classifSumm)
}


#' Determine Prophage-like read coverage elevation in whole-community
#'
#' Determines whether a detected Prophage-like genetic element has read coverage
#' in the whole-community that is either elevated or depressed compared to the
#' average read coverage of the non-prophage region.
#'
#' @param classifSummTable Classification summary table
#' @param prophageLikeClassifList A list containing pattern match information
#'   associated with all contigs classified as Prophage-like.
#' @param VLPpileup A table containing contig names, coverages averaged over 100
#'   bp windows, and contig positions associated with mapping VLP-fraction reads
#'   to whole-community contigs
#' @param WCpileup A table containing contig names, coverages averaged over 100
#'   bp windows, and contig positions associated with mapping whole-community
#'   reads to whole-community contigs
#' @param windowSize The window size used to re-average read coverage pileups
#' @return dataframe
#' @keywords internal
prophageLikeElevation <-
  function(classifSummTable,
           prophageLikeClassifList,
           VLPpileup,
           WCpileup,
           windowSize) {
    message(
      "Identifying highly active/abundant or heterogenously integrated
        Prophage-like elements"
    )
    classifSummTable$proLikeWCReadCov <- rep(NA, nrow(classifSummTable))
    classifSummTable$proLikeWCReadCovRatio <-
      rep(NA, nrow(classifSummTable))
    if (length(prophageLikeClassifList) == 0) {
      return(classifSummTable)
    }
    lapply(seq_along(prophageLikeClassifList), function(i) {
      viralSubset <- changeWindowSize(
        VLPpileup[which(VLPpileup[, 1] ==
          prophageLikeClassifList[[i]][[8]]), ],
        windowSize
      )
      startPos <- prophageLikeClassifList[[i]][[5]]
      endPos <- prophageLikeClassifList[[i]][[6]]
      contigName <- prophageLikeClassifList[[i]][[8]]
      blockLengthBp <- abs(endPos - startPos) * windowSize
      nonBlockLengthBp <-
        (nrow(viralSubset) * windowSize) - blockLengthBp
      microbialSubset <-
        changeWindowSize(
          WCpileup[which(WCpileup[, 1] ==
            contigName), ],
          windowSize
        )
      prophageLikeRegion <- microbialSubset[c(startPos:endPos), 2]
      nonProphageLikeRegion <-
        microbialSubset[which(!microbialSubset[, 2] %in%
          prophageLikeRegion), 2]
      ratio <-
        round(mean(prophageLikeRegion) / mean(nonProphageLikeRegion),
          digits = 4
        )
      classifSummTable[which(classifSummTable[, 1] ==
        contigName), 9] <<- ratio
      if (ratio > 1.15) {
        classifSummTable[which(classifSummTable[, 1] ==
          contigName), 8] <<- "Elevated"
      } else if (ratio < 0.75) {
        classifSummTable[which(classifSummTable[, 1] ==
          contigName), 8] <<- "Depressed"
      } else {
        classifSummTable[which(classifSummTable[, 1] ==
          contigName), 8] <<- "None"
      }
    })
    return(classifSummTable)
  }
