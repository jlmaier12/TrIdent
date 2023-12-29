#' Specialized transduction search and plot
#'
#' Search contigs classified as prophage-like for potential specialized transduction and return the plot visualizing the search results. Plots in green have potential specialized transduction events.
#'
#' @param phageread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#' @param transductionclassifications The pattern match information associated with each contig classified as prophage-like, gen/lat/GTA, or none with high VLP-fraction:whole-community read coverage ratios
#' @param windowsize The window size used to re-average read coverage datasets
#' @param i The index for the contig currently being assessed
#' @param noreadcov How many bp of no read coverage are encountered before specialized transduction searching stops? Default is 500.
#' @param spectranslength  How many bp of read coverage are needed for specialized transduction to be considered? Default is 2000.
#' @param transductionclassificationsummary The summary information associated with each contig classified as prophage-like, gen/lat/GTA, or none with high VLP-fraction:whole-community read coverage ratios
#' @param ref_name The reference name of the contig currently being assessed (i.e "NODE_1")
#' @keywords internal
spec_transduction_search_and_plot <- function(ref_name, phageread_dataset, transductionclassifications, transductionclassificationsummary, windowsize, i, noreadcov, spectranslength){
  position <- logcoverage <- NULL
  specialized_transduction_summary <- c(ref_name, rep(NA, 5))
  active_prophage <- ifelse(transductionclassificationsummary[which(transductionclassificationsummary[,1]== ref_name),6] == "YES", "Highly active/abundant", NA)
  viral_subset <- phageread_dataset[which(phageread_dataset[,1] == ref_name),]
  classification <- transductionclassifications[[i]][[8]]
  viral_zoom_pos <- viral_subset_zoom(viral_subset, transductionclassifications, i, 500, windowsize)
  start_pos <- viral_zoom_pos[[1]]
  end_pos <- viral_zoom_pos[[2]]
  viral_zoom <- viral_subset[c(start_pos:end_pos),]
  viral_zoom_pos <- viral_subset_zoom(viral_subset, transductionclassifications, i, 100, windowsize)
  start_pos <- viral_zoom_pos[[1]]
  end_pos <- viral_zoom_pos[[2]]
  viral_zoom_margin <- viral_subset[c(start_pos:end_pos),]
  margins <- prophagelike_border_finder(viral_zoom_margin, transductionclassifications, i, windowsize)
  start_pos_bp <-margins[[1]]
  start_pos_row <- which(viral_zoom[,3]==start_pos_bp)
  end_pos_bp <- margins[[2]]
  end_pos_row <- which(viral_zoom[,3]==end_pos_bp)
  viral_zoom$logcoverage <- abs(log10(viral_zoom[,2]))
  viral_zoom[viral_zoom == Inf] <- 0
  X <- 1
  zero_countX <- 0
  repeat{
    if((start_pos_row - X) <= 1) break
    zero_countX <- ifelse((viral_zoom[start_pos_row-X,4] == 0), (zero_countX+1), 0)
    if(zero_countX == noreadcov/100) break
    X <- X+1
  }
  Y <- 1
  zero_countY <- 0
  repeat{
    if(end_pos_row +Y >= nrow(viral_zoom)) break
    zero_countY <- ifelse((viral_zoom[(end_pos_row+Y),4] == 0), (zero_countY+1), 0)
    if(zero_countY == noreadcov/100) break
    Y <- Y+1
  }
  fill <- "deepskyblue3"
  if (X-zero_countX >= spectranslength/100) {
    transduction_start_left <- viral_zoom[start_pos_row-(X-zero_countX),3]
    transduction_left <- "Yes"
    alpha_l <- 1
    fill <- "seagreen"
    specialized_transduction_summary[3] <- "yes"
    specialized_transduction_summary[5] <- (X-zero_countX)*100
  } else {
    transduction_left <- "No"
    transduction_start_left <- start_pos_bp
    alpha_l <- 0
    specialized_transduction_summary[3] <- "no"
    specialized_transduction_summary[5] <- NA
  }
  if (Y-zero_countY >= spectranslength/100) {
    transduction_start_right <- viral_zoom[end_pos_row+(Y-zero_countY),3]
    transduction_right <- "Yes"
    alpha_r <- 1
    fill <- "seagreen"
    specialized_transduction_summary[4] <- "yes"
    specialized_transduction_summary[6] <- (Y-zero_countY)*100
  } else {
    transduction_right <- "No"
    transduction_start_right <- end_pos_bp
    alpha_r <- 0
    specialized_transduction_summary[4] <- "no"
    specialized_transduction_summary[6] <- NA
  }
  specialized_transduction_summary[2] <- ifelse((X-zero_countX >=spectranslength/100 | Y-zero_countY >=spectranslength/100), "yes", "no")
  plot <- (ggplot(data=viral_zoom, aes(x=position, y=logcoverage))+
             geom_area(fill=fill) +
             geom_vline(xintercept=c(start_pos_bp, end_pos_bp), linewidth=1)+
             geom_vline(xintercept=transduction_start_left, color="red", alpha=alpha_l, linewidth=1)+
             geom_vline(xintercept=transduction_start_right, color="red", alpha=alpha_r, linewidth=1)+
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15))+
             labs(title=paste(ref_name,classification, active_prophage), subtitle=paste0("Specialized transduction on left: ", transduction_left,", ", "on right: ", transduction_right), x="Contig Position (bp)", y="VLP-fraction \n log10 read coverage"))
  return(list(specialized_transduction_summary, plot))
}
