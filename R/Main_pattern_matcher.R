#' Main pattern-matching function
#'
#' Creates the viral_subset, representative of one contig, that is used as input for each individual pattern-building function. After the information associated with the best match for each pattern is obtained, the pattern with the lowest mean absolute difference (score) is chosen as the prediction for the contig being assessed. Prior to the pattern-building and translating, contigs less than 30,000bp and contigs where the value 50th from the max (5000bp) is less than or equal to 10 are removed. Contigs less than 45,000bp are matched against all three prophage patterns and the 'no transduction' pattern. Contigs greater than 45,000bp but less than 100,000 are matched against all three prophage patterns, the two generalized patterns with slopes that start on the contig, and the 'no transduction' pattern. Contigs greater than 100,000bp are matched against all pattern varieties.
#'
#' @param phageread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#' @param microbialread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param windowsize The window size used to re-average read coverage datasets
#' @param minblocksize The minimum size of the prophage-like block pattern. Default is 10000 bp.
#' @param maxblocksize The maximum size of the prophage-like block pattern. Default is NA
#' @keywords internal
pattern_matcher <- function (phageread_dataset, microbialread_dataset, windowsize, minblocksize, maxblocksize) {
  refnames <- unique(phageread_dataset[,1])
  best_match_list <- list()
  filteredout_contigs <- rep(NA, length(refnames))
  reason <- rep(NA, length(refnames))
  match_scoreQC <- rep(NA, length(refnames))
  A <- 1
  B <- 1
  C <- 1
  lapply(1:length(refnames), function(p) {
    i<-refnames[[p]]
    viral_subset <- phageread_dataset[which(phageread_dataset[,1] == i),]
    if(B == floor(length(refnames)/4)){
      cat("A quarter of the way done with pattern_matching \n")
    }
    if(B == floor(length(refnames)/2)){
      cat("Half of the way done with pattern_matching \n")
    }
    if(B == floor((length(refnames)*3)/4)){
      cat("Almost done with pattern_matching! \n")
    }
    B <<- B+1
    if (viral_subset[nrow(viral_subset),3]< 30000) {
      filteredout_contigs[C] <<- i
      reason[C] <<- "Contig length too small"
      C <<- C+1
      return(NULL)
      } else if (viral_subset[(order(viral_subset[,2], decreasing=TRUE))[minblocksize/100],2] <= 10) {
      filteredout_contigs[C] <<-  i
      reason[C] <<-  "Low VLP-fraction read cov"
      C <<- C+1
      return(NULL)
    }
    viral_subset <- windowsize_func(viral_subset,windowsize)
    blocks_list<-block_builder(viral_subset, windowsize, minblocksize, maxblocksize) #combined function
    if (viral_subset[nrow(viral_subset),3]< 45000) {
      best_match_summary <- list(notransduction_pattern(viral_subset),
                                 blocks_list[[1]],
                                 blocks_list[[2]],
                                 blocks_list[[3]])
      best_match_score_summary <- c(best_match_summary[[1]][[1]],best_match_summary[[2]][[1]],
                                    best_match_summary[[3]][[1]], best_match_summary[[4]][[1]]) %>% as.numeric()
    } else if (viral_subset[nrow(viral_subset),3]> 45000 & viral_subset[nrow(viral_subset),3]< 100000){ #only do gen. pattern_matching on contigs greater than 60Kbp
      slope_list<-slope_direct_withstart(viral_subset, windowsize)
      best_match_summary <- list(notransduction_pattern(viral_subset),
                                 blocks_list[[1]],
                                 blocks_list[[2]],
                                 blocks_list[[3]],
                                 slope_list[[1]],
                                 slope_list[[2]])
      best_match_score_summary <- c(best_match_summary[[1]][[1]],best_match_summary[[2]][[1]],
                                    best_match_summary[[3]][[1]],best_match_summary[[4]][[1]],
                                    best_match_summary[[5]][[1]],best_match_summary[[6]][[1]]) %>% as.numeric()
    } else {
      slope_list<-slope_direct_withstart(viral_subset, windowsize)
      slope_list_no_start<-slope_direct(viral_subset, windowsize)
      best_match_summary <- list(notransduction_pattern(viral_subset),
                                 blocks_list[[1]],
                                 blocks_list[[2]],
                                 blocks_list[[3]],
                                 slope_list_no_start[[1]],
                                 slope_list_no_start[[2]],
                                 slope_list[[1]],
                                 slope_list[[2]])
      best_match_score_summary <- c(best_match_summary[[1]][[1]],best_match_summary[[2]][[1]],
                                    best_match_summary[[3]][[1]],best_match_summary[[4]][[1]],
                                    best_match_summary[[5]][[1]],best_match_summary[[6]][[1]],
                                    best_match_summary[[7]][[1]]) %>% as.numeric()
    }
    best_match <- best_match_summary[[which(best_match_score_summary == min(best_match_score_summary))[1]]] #may need to have a way for matches to 'tie'
    best_match_list[[A]] <<- c(best_match, i)
    match_scoreQC <- c(match_scoreQC, best_match[[1]]/mean(viral_subset$coverage))
    A <<- A+1
  })
  filteredout_contigs <- filteredout_contigs[!is.na(filteredout_contigs)]
  reason <- reason[!is.na(reason)]
  match_scoreQC <- match_scoreQC[!is.na(match_scoreQC)]
  match_scoreQC <- as.data.frame(match_scoreQC)
  colnames(match_scoreQC) <- "Match_score"
  print(ggplot(data=match_scoreQC, aes(x=Match_score))+
          geom_density()+
          labs(title="Match-score quality threshold", caption="(Lower scores are better matches)")+
          theme_bw())
  filteredout_summary_df <- cbind.data.frame(filteredout_contigs, reason)
  pattern_matching_summary <- list(best_match_list, filteredout_summary_df)
  return(pattern_matching_summary)
}
