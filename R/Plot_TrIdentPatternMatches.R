#' Plot read coverage graphs of contigs classified as Prophage-like, Gen/Lat/GTA, or HighVLPWCRC
#'
#' Plot the read coverages of a contig and its associated pattern-match for each contig classified as prophage-like, gen/lat/GTA or none with high VLP-fraction:whole-community read coverage ratios. Returns a list of plot objects.
#'
#' @param phageread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#' @param microbialread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param transductionclassifications Output from TrIdent_Classifier
#' @param windowsize The window size used to re-average read coverage datasets in the TrIdent_Classifier. Default is 1000.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pattern_matches <- Plot_TrIdentPatternMatches(VLP_fracreadcov, whole_commreadcov, TrIdent_results, windowsize=1000)
#' pattern_matches$NODE_1
#' }
Plot_TrIdentPatternMatches <- function(phageread_dataset, microbialread_dataset, transductionclassifications, windowsize=1000) {
  position <- coverage <- NULL
  transductionclassifications_wlowratios <- transductionclassifications[[3]]
  final_summary_table <- transductionclassifications[[1]]
  phageread_dataset <- readcovdf_formatter(phageread_dataset)
  microbialread_dataset <- readcovdf_formatter(microbialread_dataset)
  plots <- list()
  ref_names <- c()
  for (i in seq(1,length(transductionclassifications_wlowratios),1)) {
    ref_name <- transductionclassifications_wlowratios[[i]][[8]]
    ref_names <- c(ref_names, ref_name)
    viral_subset <- phageread_dataset[which(phageread_dataset[,1] == ref_name),]
    viral_subset <- windowsize_func(viral_subset,windowsize)
    viral_subset[is.nan.data.frame(viral_subset)] <- 0
    microbial_subset <- microbialread_dataset[which(microbialread_dataset[,1] == ref_name),]
    microbial_subset <- windowsize_func(microbial_subset,windowsize)
    microbial_subset[is.nan.data.frame(microbial_subset)] <- 0
    match_info <- final_summary_table[which(final_summary_table[,1]==ref_name),]
    classification <- match_info[,2]
    pattern <- pattern_builder(viral_subset, transductionclassifications_wlowratios, classification, i)
    pattern_match <- cbind(viral_subset, pattern)
    match_length <- match_info[,3]
    classification <- match_info[,2]
    if (is.na(match_info[4]) == TRUE){
      prophage_activity <- NULL
    } else if (match_info[4]=="YES"){
      prophage_activity <- "Active/highly abundant prophage"
    } else {
      prophage_activity <- NULL
    }
    wholecomm_plot <- ggplot(data=microbial_subset, aes(x=position, y=coverage))+
      geom_area(fill="deepskyblue3") +
      labs(title=paste(ref_name, "Classification", classification ),subtitle=paste("Matching-region size (bp):", match_length, prophage_activity), x=" ", y="Whole-community \n read coverage") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15))

    Overlay_plot <- ggplot(data=pattern_match, aes(x=position, y=coverage))+
      geom_area(fill="deepskyblue3") +
      geom_line(aes(y=pattern), color="black", size=1)+
      labs(x="Contig position (bp)", y="VLP-fraction \n read coverage")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15))

    combined_plot <- (wholecomm_plot/Overlay_plot)
    plots[[i]] <- combined_plot
  }
  names(plots) <- ref_names
  return(plots)
}
