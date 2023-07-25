#' Prophage-like zoom-in
#'
#' Subset contigs with prophage-like patterns to 'zoom-in' on desired region surrounding pattern.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param transductionclassificationpatterns The pattern match information associated with each contig classified as prophage-like, gen/lat/GTA, or none with high VLP-fraction:whole-community read coverage ratios
#' @param i The index for the contig currently being assessed
#' @param zoom The number of rows outside the start and stop positions of the prophage-like pattern to zoom-in on
#' @param windowsize The window size used to re-average read coverage datasets
#' @keywords internal
viral_subset_zoom <- function (viral_subset, transductionclassificationpatterns, i, zoom, windowsize) {
  start_pos_row <- transductionclassificationpatterns[[i]][[5]] *windowsize/100
  end_pos_row <- transductionclassificationpatterns[[i]][[6]] *windowsize/100
  start_pos <- ifelse((start_pos_row - zoom) < 1, 1, (start_pos_row -zoom))
  end_pos <- ifelse((end_pos_row +zoom > nrow(viral_subset)), nrow(viral_subset), (end_pos_row+ zoom))
  return(list(start_pos, end_pos))
}
