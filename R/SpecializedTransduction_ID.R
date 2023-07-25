#' Identify potential specialized transduction events on contigs classified as Prophage-like
#'
#' Search contigs classified as prophage-like for potential specialized transduction. Returns a list with the first object containing a summary table and the second object containing a list of plots of prophage-like patterns with associated specialzied transduction search results. If the plot is green, it has been identified as having potential specialized transduction.
#'
#' @param phageread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#' @param transductionclassification Output from TrIdent_Classifier
#' @param windowsize The window size used to re-average read coverage datasets. Default is 1000
#' @param noreadcov How many bp of no read coverage are encountered before specialized transduction searching stops? Default is 500.
#' @param spectranslength How many bp of read coverage are needed for specialized transduction to be considered? Default is 2000.
#' @param specificcontig Provide the name of a specific contig if you would like to search only that contig.i.e. "NODE_1"
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Spec_Transduction <- SpecializedTransduction_ID(VLP_fracreadcov, TrIdent_results, windowsize=1000)
#' Summary_Table <- spec_transduction[[1]]
#' NODE_8_plot <- spec_transduction[[2]]$NODE_8
#'
#' Spec_Transduction_NODE1 <- SpecializedTransduction_ID(VLP_fracreadcov, TrIdent_results, windowsize=1000, noreadcov=1000, spectranslength=4000, specificcontig="NODE_1")
#' }
 SpecializedTransduction_ID <- function(phageread_dataset, transductionclassification, windowsize=1000, noreadcov=500, spectranslength=2000, specificcontig){
  transductionclassificationpatterns <- transductionclassification[[3]]
  transductionclassificationsummary <- transductionclassification[[1]]
  phageread_dataset <- readcovdf_formatter(phageread_dataset)
  if(missing(specificcontig)) {
    potential_spec_transduction <- 0
    specialized_transduction_summary <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(specialized_transduction_summary) <- c("ref_name", "Specialized_transduction", "Left", "Right")
    plots <- list()
    J <- 1
    for (i in  seq(1, length(transductionclassificationpatterns), 1)) {
      classification <- transductionclassificationpatterns[[i]][[7]]
      if (classification != "Prophage-like") next
      ref_name <- transductionclassificationpatterns[[i]][[8]]
      spec_trans_pred <- spec_transduction_search_and_plot(ref_name, phageread_dataset, transductionclassificationpatterns, transductionclassificationsummary, windowsize, i, noreadcov, spectranslength)
      if(spec_trans_pred[[1]][[2]] == "yes") {
        potential_spec_transduction <- potential_spec_transduction +1}
      specialized_transduction_summary[J,c(1:4)] <- spec_trans_pred[[1]]
      plots[[J]] <- spec_trans_pred[[2]]
      J <- J+1
    }
    print(paste(potential_spec_transduction, "contigs had potential specialized transduction"))
    names(plots) <- specialized_transduction_summary[,1]
    return(list(specialized_transduction_summary,plots))
  } else {
    for (i in  seq(1, length(transductionclassificationpatterns), 1)) {
      ref_name <- transductionclassificationpatterns[[i]][[8]]
      if (ref_name == specificcontig) {
        spec_trans_pred <- spec_transduction_search_and_plot(ref_name, phageread_dataset, transductionclassificationpatterns, transductionclassificationsummary, windowsize, i, noreadcov, spectranslength)
        return(spec_trans_pred[[2]])
      }
    }
  }
}
