#' Search/explore gene annotations of contigs classified as having transduction events
#'
#' Search contigs classified as prophage-like, gen/lat/GTA or none with high VLP-fraction:whole-community read coverage ratios for gene-annotations that match a provided key-word(s). Outputs read coverage plots for contigs with matching annotations.
#'
#' @param transductionclassification Output from TrIdent_Classifier.
#' @param phageread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#' @param windowsize The window size used to re-average read coverage datasets in TrIdent_Classifier. Default is 1000.
#' @param genelocation Would you like to search for gene-annotations in a specific region surrounding the location of the predicted transduction event? If yes, "specific", if no "nonspecific". Default is "nonspecific" (i.e the entire contig, including any regions not associated with the matching region of the predicted transduction event, will be searched for the gene annotation key-words)
#' @param bprange If you are searching for gene annotations in specific ("specific") locations, you may specify the region (in basepairs) that should be searched to the left and right of the predicted transduction event location on the contig. Default is 0 (if you choose specific and don't set the bprange, it will only search gene annotations within the pattern-match region)
#' @param gene_annots GFF3 file cleaned with the gff3_cleanup function
#' @param geneorproduct Are you searching for a gene name or gene product? "gene" or "product"
#' @param keywords The key-word(s) you would like to search for. Case independent. Searches will return the full gene annotations that contain the matching key-word. Key-word(s) must be in quotes, comma-separated, and surrounded by c() i.e( c("antibiotic", "resistance", "drug") )
#' @param specificcontig Provide the name of a specific contig if you would like to search only that contig.i.e. "NODE_1"
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' Virulence_matches <- GeneAnnotationSearch(TrIdent_results, whole_commreadcov, VLP_fracreadcov, cleaned_gff, "product", windowsize=1000, keywords=c("resistance", "antibiotic", "toxin"), genelocation="specific", bprange=10000)
#' Virulence_matches$NODE_450
#'
#' CmR_matches <- GeneAnnotationSearch(TrIdent_results, whole_commreadcov, VLP_fracreadcov, cleaned_gff, "gene", windowsize=1000, keywords=c("cmR"), genelocation="specific", bprange=0,)
#' CmR_matches$NODE_10
#' Contig_1_GeneAnnots <- GeneAnnotationSearch(TrIdent_results, whole_commreadcov, VLP_fracreadcov, cleaned_gff, "gene", windowsize=1000, keywords=c("cmR"), genelocation="nonspecific",specificcontig="NODE_1")
#' }
GeneAnnotationSearch <- function(transductionclassification, phageread_dataset, gene_annots, geneorproduct, windowsize, keywords, genelocation="nonspecific", bprange = 0, specificcontig) {
  transductionclassificationpatterns <- transductionclassification[[3]]
  phageread_dataset <- readcovdf_formatter(phageread_dataset)
  if(missing(specificcontig)==TRUE) {
    plots <- list()
    ref_names <- c()
    X <- 1
    for (i in  seq(1, length(transductionclassificationpatterns), 1)) {
      ref_name <- transductionclassificationpatterns[[i]][[8]]
      viral_subset <- phageread_dataset[which(phageread_dataset[,1] == ref_name),]
      subset_gene_annots <- gene_annots[which(gene_annots[,1]== ref_name),]
      classification <- transductionclassificationpatterns[[i]][[7]]
      if (genelocation== "specific" & transductionclassificationpatterns[[i]][[5]] =="NA") next
      colindex <- ifelse(geneorproduct=="gene", 7, 8)
      if (TRUE%in% (str_detect(subset_gene_annots[,colindex], regex(paste(keywords, collapse="|"),ignore_case=T)))==TRUE) {
        if (genelocation=="specific"){
          start_pos <- transductionclassificationpatterns[[i]][[5]] *windowsize
          end_pos <- transductionclassificationpatterns[[i]][[6]] *windowsize
          match_indexes <- str_which(subset_gene_annots[,colindex], regex(paste(keywords, collapse="|"), ignore_case=T))
          matching_gene_annots <- subset_gene_annots[match_indexes,]
          start_bprange <- ifelse((start_pos-bprange < 1),1,(start_pos - bprange))
          end_bprange <- ifelse ((end_pos +bprange >(viral_subset[nrow(viral_subset),3])), (viral_subset[nrow(viral_subset),3]), (end_pos + bprange))
          matching_gene_annot_subset <- matching_gene_annots[((which(matching_gene_annots[,2] %in% c(start_bprange:end_bprange)))),]
          if (TRUE%in% (str_detect(matching_gene_annot_subset[,colindex], regex(paste(keywords, collapse="|"),ignore_case=T)))==TRUE) {
            plots[[X]] <- gene_annotation_search_and_plot(subset_gene_annots, keywords, viral_subset, classification, ref_name, colindex, start_bprange, end_bprange)
            ref_names <- c(ref_names, ref_name)
            X <- X+1
          }
        } else {
          plots[[X]] <- gene_annotation_search_and_plot(subset_gene_annots, keywords, viral_subset, classification, ref_name, colindex)
          ref_names <- c(ref_names, ref_name)
          X <- X+1
        }
      }
    }
    names(plots) <- ref_names
    print(paste(length(ref_names), "contigs had gene annotations that match one or more of the provided keywords"))
    return(plots)
  }
  else {
    for (i in  seq(1, length(transductionclassificationpatterns), 1)) {
      ref_name <- transductionclassificationpatterns[[i]][[8]]
      if (ref_name == specificcontig) {
        ref_name <- transductionclassificationpatterns[[i]][[8]]
        viral_subset <- phageread_dataset[which(phageread_dataset[,1] == ref_name),]
        subset_gene_annots <- gene_annots[which(gene_annots[,1]== ref_name),]
        classification <- transductionclassificationpatterns[[i]][[7]]
        if (genelocation== "specific" & transductionclassificationpatterns[[i]][[5]] =="NA"){
          print("Cannot use a 'specific' search on a contig with a pattern-match that spans the entire contig")
          next
        }
        colindex <- ifelse(geneorproduct=="gene", 7, 8)
        if (TRUE%in% (str_detect(subset_gene_annots[,colindex], regex(paste(keywords, collapse="|"),ignore_case=T)))==TRUE) {
          if (genelocation=="specific"){
            start_pos <- transductionclassificationpatterns[[i]][[5]] *windowsize
            end_pos <- transductionclassificationpatterns[[i]][[6]] *windowsize
            match_indexes <- str_which(subset_gene_annots[,colindex], regex(paste(keywords, collapse="|"), ignore_case=T))
            matching_gene_annots <- subset_gene_annots[match_indexes,]
            start_bprange <- ifelse((start_pos-bprange < 1),1,(start_pos - bprange))
            end_bprange <- ifelse ((end_pos +bprange >(viral_subset[nrow(viral_subset),3])), (viral_subset[nrow(viral_subset),3]), (end_pos + bprange))
            matching_gene_annot_subset <- matching_gene_annots[((which(matching_gene_annots[,2] %in% c(start_bprange:end_bprange)))),]
            if (TRUE%in% (str_detect(matching_gene_annot_subset[,colindex], regex(paste(keywords, collapse="|"),ignore_case=T)))==TRUE) {
              plot_geneannots <- gene_annotation_search_and_plot(subset_gene_annots, keywords, viral_subset, classification, ref_name, colindex)
              return(plot_geneannots)
            } else {
              print("No gene annotations match the provided keywords within the specific search range")
            }
          }else {
            plot_geneannots <- gene_annotation_search_and_plot(subset_gene_annots, keywords, viral_subset, classification, ref_name, colindex)
            return(plot_geneannots)
          }
        }
        else {
          print("No gene annotations match the provided keywords")
        }
      }
    }
  }
}
