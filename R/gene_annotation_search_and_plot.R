#' Gene annotation search and plot
#'
#' Search contigs for gene annotations that match the provided keywords and plot the contigs with matches and their respective locations on the contig
#'
#' @param subset_gene_annots The gene annotations associated with the contig being searched
#' @param keywords The key-word(s) you would like to search for. Case independent. Searches will return the full gene annotations that contain the matching key-word. Key-word(s) must be in quotes, comma-separated, and surrounded by c() i.e( c("antibiotic", "resistance", "drug") )
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param classification The classification of the contig being searched
#' @param ref_name The reference name for the contig being searched
#' @param colindex The column index of either 'gene' or 'product'
#' @param start_bprange The basepair at which the search is started if a 'specific' search is used
#' @param end_bprange The basepair at which the search is ended if a 'specific' search is used
#' @keywords internal
gene_annotation_search_and_plot <- function(subset_gene_annots, keywords, viral_subset, classification, ref_name, colindex, start_bprange, end_bprange){
  position <- logcoverage <- V4 <- NULL
  if (missing(start_bprange)==TRUE & missing(end_bprange) == TRUE){
    start_bprange <- NULL
    end_bprange <- NULL
  }
  match_indexes <- str_which(subset_gene_annots[,colindex], regex(paste(keywords, collapse="|"), ignore_case=T))
  gene_annot_subset <- subset_gene_annots[match_indexes,]
  gene_start_pos <- gene_annot_subset[,2]
  gene_annots_labels <- paste0("#", c(1:nrow(gene_annot_subset)), ": ",gene_annot_subset[,colindex], sep=" ", collapse=" \n ")
  viral_subset$logcoverage <- abs(log10(viral_subset[,2]))
  viral_subset[viral_subset == Inf] <- 0
  plot <- ggplot(data=viral_subset, aes(x=position, y=logcoverage))+
    geom_area(fill="deepskyblue3") +
    geom_vline(xintercept=gene_start_pos, size=1)+
    geom_vline(xintercept=c(start_bprange, end_bprange), color="red", size=0.5)+
    geom_label(data=gene_annot_subset, aes(x=V4,y=(max(viral_subset$logcoverage)/2),label=paste0("#", c(1:nrow(gene_annot_subset)))), size=2.5)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 15))+
    labs(title=paste(ref_name,classification), x="Contig position (bp)", caption=gene_annots_labels, y="VLP-fraction \n read coverage")
  return(plot)
}
