#' Main pattern-matching function
#'
#' Creates the viral_subset, representative of one contig, that is used as input for each individual pattern-building function. After the information associated with the best match for each pattern is obtained, the pattern with the lowest mean absolute difference (score) is chosen as the prediction for the contig being assessed. Prior to the pattern-building and translating, contigs less than 30,000bp and contigs where the value 50th from the max (5000bp) is less than or equal to 10 are removed. Contigs less than 45,000bp are matched against all three prophage patterns and the 'no transduction' pattern. Contigs greater than 45,000bp but less than 100,000 are matched against all three prophage patterns, the two generalized patterns with slopes that start on the contig, and the 'no transduction' pattern. Contigs greater than 100,000bp are matched against all pattern varieties.
#'
#' @param phageread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#' @param microbialread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param windowsize The window size used to re-average read coverage datasets
#' @keywords internal
pattern_matcher <- function (phageread_dataset, microbialread_dataset, windowsize) {
  refnames <- unique(phageread_dataset[,1])
  best_match_list <- list()
  filteredout_contigs <- rep(NA, length(refnames))
  reason <- rep(NA, length(refnames))
  A <- 1
  B <- 1
  C <- 1
  for (i in refnames) {
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
    B <- B+1
    if (viral_subset[nrow(viral_subset),3]< 30000) {
      filteredout_contigs[C] <- i
      reason[C] <- "Contig length too small"
      C <- C+1
      next
      } else if (viral_subset[(order(viral_subset[,2], decreasing=TRUE))[50],2] <= 10) {
      filteredout_contigs[C] <-  i
      reason[C] <-  "Low VLP-fraction read cov"
      C <- C+1
      next
    }
    viral_subset <- windowsize_func(viral_subset,windowsize)
    viral_subset[is.nan.data.frame(viral_subset)] <- 0
    microbial_subset <- microbialread_dataset[which(microbialread_dataset[,1] == i),]
    microbial_subset <- windowsize_func(microbial_subset,windowsize)
    microbial_subset[is.nan.data.frame(microbial_subset)] <- 0
    if (viral_subset[nrow(viral_subset),3]< 45000) {
      no_transduction_best_match <- notransduction_pattern(viral_subset)
      prophage_off_left_best_match <- block_off_left_translator(viral_subset, windowsize)
      prophage_off_right_best_match <-  block_off_right_translator(viral_subset, windowsize)
      full_prophage_best_match <- full_blockpattern_builder(viral_subset, windowsize)
      best_match_summary <- list(no_transduction_best_match, prophage_off_left_best_match, prophage_off_right_best_match, full_prophage_best_match)
      best_match_score_summary <- c(no_transduction_best_match[[1]],prophage_off_left_best_match[[1]], prophage_off_right_best_match[[1]], full_prophage_best_match[[1]]) %>% as.numeric()
    } else if (viral_subset[nrow(viral_subset),3]> 45000 & viral_subset[nrow(viral_subset),3]< 100000){ #only do gen. pattern_matching on contigs greater than 60Kbp
      no_transduction_best_match <- notransduction_pattern(viral_subset)
      prophage_off_left_best_match <- block_off_left_translator(viral_subset, windowsize)
      prophage_off_right_best_match <-  block_off_right_translator(viral_subset, windowsize)
      full_prophage_best_match <- full_blockpattern_builder(viral_subset, windowsize)
      Gen_LR_wstart_best_match <- slope_LefttoRight_withstart(viral_subset, windowsize)
      Gen_RL_wstart_best_match <- slope_RighttoLeft_withstart(viral_subset, windowsize)
      best_match_summary <- list(no_transduction_best_match, prophage_off_left_best_match, prophage_off_right_best_match, full_prophage_best_match, Gen_LR_wstart_best_match, Gen_RL_wstart_best_match)
      best_match_score_summary <- c(no_transduction_best_match[[1]],prophage_off_left_best_match[[1]], prophage_off_right_best_match[[1]], full_prophage_best_match[[1]], Gen_LR_wstart_best_match[[1]], Gen_RL_wstart_best_match[[1]]) %>% as.numeric()
    } else {
      no_transduction_best_match <- notransduction_pattern(viral_subset)
      prophage_off_left_best_match <- block_off_left_translator(viral_subset, windowsize)
      prophage_off_right_best_match <-  block_off_right_translator(viral_subset, windowsize)
      full_prophage_best_match <- full_blockpattern_builder(viral_subset, windowsize)
      Gen_LR_best_match <- slope_LefttoRight(viral_subset, windowsize)
      Gen_RL_best_match <- slope_RighttoLeft(viral_subset, windowsize)
      Gen_LR_wstart_best_match <- slope_LefttoRight_withstart(viral_subset, windowsize)
      Gen_RL_wstart_best_match <- slope_RighttoLeft_withstart(viral_subset, windowsize)
      best_match_summary <- list(no_transduction_best_match, prophage_off_left_best_match, prophage_off_right_best_match, full_prophage_best_match, Gen_LR_best_match, Gen_RL_best_match, Gen_LR_wstart_best_match, Gen_RL_wstart_best_match)
      best_match_score_summary <- c(no_transduction_best_match[[1]],prophage_off_left_best_match[[1]], prophage_off_right_best_match[[1]], full_prophage_best_match[[1]], Gen_LR_best_match[[1]], Gen_RL_best_match[[1]], Gen_LR_wstart_best_match[[1]], Gen_RL_wstart_best_match[[1]]) %>% as.numeric()
    }
    best_match <- best_match_summary[[which(best_match_score_summary == min(best_match_score_summary))[1]]] #may need to have a way for matches to 'tie'
    best_match_list[[A]] <- append(best_match, i)
    A <- A+1
  }
  filteredout_contigs <- filteredout_contigs[!is.na(filteredout_contigs)]
  reason <- reason[!is.na(reason)]
  filteredout_summary_df <- cbind.data.frame(filteredout_contigs, reason)
  pattern_matching_summary <- list(best_match_list, filteredout_summary_df)
  return(pattern_matching_summary)
}
