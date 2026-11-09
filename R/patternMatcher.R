#' Main pattern-matching function
#'
#' Creates the viralSubset, representative of one contig, that is used as input
#' for each individual pattern-matching function. After the information
#' associated with the best match for each pattern is obtained, the pattern with
#' the smallest match score is used to classify the contig being assessed. Prior
#' to the pattern-matching, contigs smaller than the minContigLength and contigs
#' without 5,000 bp of 10x read coverage are removed.
#'
#' @param VLPpileup A table containing contig names, coverages averaged over 100
#'  bp windows, and contig positions associated with mapping VLP-fraction reads
#'  to whole-community contigs
#' @param WCpileup A table containing contig names, coverages averaged over 100
#'  bp windows, and contig positions associated with mapping whole-community
#'  reads to whole-community contigs
#' @param windowSize The window size used to re-average read coverage datasets
#' @param minBlockSize The minimum size of the prophage-like block pattern.
#'  Default is 10,000 bp.
#' @param maxBlockSize The maximum size of the prophage-like block pattern.
#'  Default is NA
#' @param minContigLength The minimum contig size (in bp) to perform
#'  pattern-matching on. Must be at least 20,000 bp. Default is 30,000 bp.
#' @param minSlope The minimum slope value to test for sloping patterns
#' @param minSlopeSize The minimum width of sloping patterns. 
#' @param verbose TRUE or FALSE. Print progress messages to console. Default is TRUE.
#' @param searchMethod Search method to use. Either "grid" for the original grid
#'   search or "direct" for DIRECT global optimization.
#' @param DirectMaxEval Maximum number of DIRECT evaluations to make.
#' @param globalLocal Use global or local DIRECT search. Default is local. 
#' @return List containing three objects.
#' @keywords internal
patternMatcher <-
  function(VLPpileup,
           WCpileup,
           windowSize,
           minBlockSize,
           maxBlockSize,
           minContigLength,
           minSlope,
           minSlopeSize,
           verbose,
           searchMethod,
           DirectMaxEval,
           globalLocal) {
    contigNames <- unique(VLPpileup[, 1])
    bestMatchList <- list()
    filteredOutContigs <- rep(NA, length(contigNames))
    reason <- rep(NA, length(contigNames))
    refs <- rep(NA, length(contigNames))
    A <- 1
    B <- 1
    C <- 1
    for (p in seq_along(contigNames)) {
      i <- contigNames[[p]]
      viralSubset <- VLPpileup[which(VLPpileup[, 1] == i), ]
      if(verbose){
      if (B == floor(length(contigNames) / 4)) {
        message("A quarter of the way done with pattern-matching")
      }
      if (B == floor(length(contigNames) / 2)) {
        message("Half of the way done with pattern-matching")
      }
      if (B == floor((length(contigNames) * 3) / 4)) {
        message("Almost done with pattern-matching!")
      }
      }
      B <- B + 1
      if (viralSubset[nrow(viralSubset), 3] < minContigLength) {
        filteredOutContigs[C] <- i
        reason[C] <- "Contig length too small"
        C <- C + 1
        next
      } else if (viralSubset[(order(viralSubset[, 2],
        decreasing = TRUE))[minBlockSize / 100], 2] <= 5) {
        filteredOutContigs[C] <- i
        reason[C] <- "Low VLP-fraction read cov"
        C <- C + 1
        next
      }
      viralSubset <- changeWindowSize(viralSubset, windowSize)
      if (length(unique(viralSubset[, 2])) != 1) {
        blocksList <- blockBuilder(
          viralSubset,
          windowSize,
          minBlockSize,
          maxBlockSize,
          searchMethod,
          DirectMaxEval,
          globalLocal
        )
      }
      if (length(unique(viralSubset[, 2])) == 1) {
        bestMatchSumm <- list(noPattern(viralSubset, searchMethod, DirectMaxEval, globalLocal))
        bestMatchScoreSumm <-
          c(bestMatchSumm[[1]][[1]]) %>% as.numeric()
      } else {
        slopeList <- slopeWithStart(viralSubset, windowSize, minSlope, 
                                    minSlopeSize, searchMethod, DirectMaxEval, globalLocal)
        slopeListNoStart <-
          fullSlope(viralSubset, windowSize, minSlope, 
                    minSlopeSize, searchMethod, DirectMaxEval, globalLocal)
        bestMatchSumm <- list(
          noPattern(viralSubset, searchMethod, DirectMaxEval, globalLocal),
          blocksList[[1]],
          blocksList[[2]],
          blocksList[[3]],
          slopeListNoStart[[1]],
          slopeListNoStart[[2]],
          slopeList[[1]],
          slopeList[[2]]
        )
        bestMatchScoreSumm <-
          c(
            bestMatchSumm[[1]][[1]],
            bestMatchSumm[[2]][[1]],
            bestMatchSumm[[3]][[1]],
            bestMatchSumm[[4]][[1]],
            bestMatchSumm[[5]][[1]],
            bestMatchSumm[[6]][[1]],
            bestMatchSumm[[7]][[1]],
            bestMatchSumm[[8]][[1]]
          ) %>% as.numeric()
      }
      #
      # diffVec <- abs(bestMatchScoreSumm[1] - bestMatchScoreSumm[c(2:8)])
      # diffVecNorm <- (diffVec-min(diffVec))/(max(diffVec)-min(diffVec))*100
      bestMatch <- bestMatchSumm[[which(bestMatchScoreSumm ==
                   min(bestMatchScoreSumm))[1]]]
      # if(bestMatch[[7]]=="NoPattern" & any(diffVecNorm<=0.5) ) {
      #     newMatch <- min(which(diffVecNorm<=0.5))[[1]]
      #     bestMatch <- bestMatchSumm[[newMatch]]
      # }
      #
      bestMatchList[[A]] <- c(bestMatch, i)
      names(bestMatchList)[A] <- i
      names(bestMatchList[[A]])[1:8] <- c("match_score", "min_cov", "max_cov", "slope", "start_pos", "end_pos", "classif.", "contig")
      refs[A] <- i
      A <- A + 1
    }
    filteredOutContigs <- filteredOutContigs[!is.na(filteredOutContigs)]
    reason <- reason[!is.na(reason)]
    refs <- refs[!is.na(refs)]
    filteredOutSummaryTable <-
      cbind.data.frame(filteredOutContigs, reason)
    patternMatchSummaryList <-
      list(
        bestMatchList,
        filteredOutSummaryTable
      )
    return(patternMatchSummaryList)
  }
