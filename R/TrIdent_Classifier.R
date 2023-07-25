#' Classify contigs as Prophage-like, Gen/Lat/GTA, None, or HighVLPWCratio
#'
#' Classifies contigs as prophage-like, gen/lat/GTA, none, or highVLPWCratio. This function performs all the pattern-matching and summarizes the results into a list. The first item in the list is a table consisting of the summary information of all the contigs that passed through pattern-matching (i.e were not filtered out). The second item in the list is a table consisting of the summary information of all contigs that were predicted as containing a potential prophage, generalized transduction and/or low whole-community: VLP-fraction read coverage ratios. The third item in the list contains the best pattern-match information associated with each contig in the previous table. The fourth and final object in the list is a table containing the contigs that were filtered out prior to pattern_matching and the reason why.
#'
#'@param phageread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#'@param microbial_readdataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#'@param windowsize The number of base pairs to average coverage values over. Options are 100, 200, 500, 1000 and 2000 ONLY. Larger window sizes improve processing time but the resolution of read coverage patterns may be lost. Default is 1000bp.
#'
#'@export
#'
#'@examples \dontrun{
#'TrIdent_results <- TrIdent_Classifier(VLP_fracreadcov, whole_commreadcov, windowsize=1000)
#'}
  TrIdent_Classifier <- function(phageread_dataset, microbial_readdataset, windowsize = 1000){
  start_time <- Sys.time()
  windowsize <- windowsize
  viral_readdataset <- readcovdf_formatter(phageread_dataset)
  microbial_readdataset <- readcovdf_formatter(microbial_readdataset)
  print("Starting pattern-matching")
  SM_classifications_summary <- pattern_matcher(viral_readdataset, microbial_readdataset, windowsize)
  SM_classifications <- SM_classifications_summary[[1]]
  filteredoutcontigs <- SM_classifications_summary[[2]]
  print("Identifying potential transducing events")
  prophage_classifications_list <- allprophagelike_matches(SM_classifications)
  alltransduction_classifications_list <- alltransduction_events_summarylist(SM_classifications)
  classification_summary_df <- contig_classification_summary(SM_classifications)
  summary_table_withratios <- WCVF_ratio_calculator(classification_summary_df, microbial_readdataset, viral_readdataset)
  print("Determining sizes (bp) of potential transduction events")
  summary_table_matchsize <- matchsize_checker(summary_table_withratios, alltransduction_classifications_list, windowsize)
  print("Identifying highly active/abundant prophage-like elements")
  summary_table_final <- prophagelikeactivity_checker(summary_table_matchsize, prophage_classifications_list, viral_readdataset, microbial_readdataset, windowsize)
  print("Finalizing output")
  contigswith_transductioneventsandlowratios_list <- alltransduction_withlowratios_summarylist(SM_classifications, summary_table_final)
  cleaned_summary_table <- summary_table_final[which(summary_table_final[,2]=="Prophage-like"|summary_table_final[,2]=="Gen/Lat/GTA"|summary_table_final[,2]=="HighVLPWCReadCov"),]
  final_summary_list<-list(summary_table_final,cleaned_summary_table, contigswith_transductioneventsandlowratios_list, filteredoutcontigs)
  names(final_summary_list) <- c("Full_summary_table", "Cleaned_summary_table", "patternMatchInfo", "FilteredOut_contig_table")
  end_time <- Sys.time()
  print(paste("Execuion time:", end_time-start_time))
  print(paste(nrow(filteredoutcontigs), "contigs were filtered out based on short length or low read coverage"))
  print(table(final_summary_list[[1]][,2]))
  print(paste(length(which(final_summary_list[[1]][,4]=="YES")), "of the predicted prophages are highly active or abundant"))
  return(final_summary_list)
}
