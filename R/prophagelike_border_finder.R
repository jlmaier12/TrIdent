#' Prophage-like border finder
#'
#' Find borders of prophage-like patterns with more specificity than pattern-matching
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param transductionclassificationpatterns The pattern match information associated with each contig classified as prophage-like, gen/lat/GTA, or none with high VLP-fraction:whole-community read coverage ratios
#' @param i The index for the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @keywords internal
prophagelike_border_finder <- function(viral_subset,transductionclassificationpatterns,i, windowsize) {
  max_value <- viral_subset[order(viral_subset[,2], decreasing=TRUE),2][2]
  left_margin_rowpos <- transductionclassificationpatterns[[i]][[5]] *windowsize/100
  right_margin_rowpos <- transductionclassificationpatterns[[i]][[6]] *windowsize/100
  left_margin_bppos <- transductionclassificationpatterns[[i]][[5]] *windowsize
  right_margin_bppos <- transductionclassificationpatterns[[i]][[6]] *windowsize
  X <- 1
  repeat {
    if (viral_subset[X,2] >= 0.2*max_value) {
      left_margin_bppos <- viral_subset[X,3]
      left_margin_rowpos <- X
      break
    }
    X <- X+1
    if(X >= nrow(viral_subset)) break
  }
  Y <- nrow(viral_subset)
  repeat {
    if (viral_subset[Y,2] >= 0.2*max_value) {
      right_margin_bppos <- viral_subset[Y,3]
      right_margin_rowpos <- Y
      break
    }
    Y <- Y-1
    if(Y <= 1) break
  }
  margins <- list(left_margin_bppos, right_margin_bppos, left_margin_rowpos, right_margin_rowpos)
  return(margins)
}
