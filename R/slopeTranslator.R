#' Sloping pattern translator
#'
#' Translates a sloping pattern containing the initial jump-up in read coverage across a contig.
#' Translate the pattern 1000 bp at a time. Stop translating when the pattern left on the contig reaches 20,000 bp.
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to the contig currently being assessed
#' @param bestMatchInfo The pattern-match information associated with the current best pattern match.
#' @param windowSize The window size used to re-average read coverage pileups
#' @param pattern A vector containing the values associated with the sloping pattern
#' @param leftOrRight The direction of the sloping pattern. Either "Left" for left to right (neg) slopes or "Right" for right to left (pos) slopes.
#' @keywords internal
slopeTranslator <- function(viralSubset, bestMatchInfo, windowSize, slopeChange, leftOrRight){
  pattern <- slopeChange[[1]]
  minPattern <- min(pattern)
  repeat {
    if (leftOrRight == "Left") pattern <- c(rep(minPattern, 1000 / windowSize), pattern[-c((length(pattern) - ((1000 / windowSize) - 1)):length(pattern))])
    if (leftOrRight == "Right") pattern <- c(pattern[-c(1:(1000 / windowSize))], rep(minPattern, 1000 / windowSize))
    slopeBottomIdx <- min(pattern[pattern != min(pattern)])
    startRowIdx <- ifelse(leftOrRight == "Left", which(pattern == max(pattern)), which(pattern == slopeBottomIdx))
    endRowIdx <- ifelse(leftOrRight == "Left", which(pattern == slopeBottomIdx), which(pattern == max(pattern)))
    if((length(pattern[!(pattern %in% minPattern)]) * windowSize) < 20000) break
    diff <- mean(abs(viralSubset[,2] - pattern))
    if (diff < bestMatchInfo[[1]]) {
      covSteps <- ((max(pattern) - slopeBottomIdx) / abs(endRowIdx - startRowIdx))
      covSteps <- ifelse(leftOrRight == "Left", -covSteps, covSteps)
      bestMatchInfo <- list(diff, slopeBottomIdx, max(pattern), covSteps, startRowIdx, endRowIdx, slopeChange[[2]], "Sloping")
    }
  }
  return(bestMatchInfo)
}

