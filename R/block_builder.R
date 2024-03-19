#' Block Builder
#'
#' Build and translate a block pattern going off the left side, right side and full length of the contig.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @param minblocksize The minimum size of the prophage-like block pattern. Default is 10000 bp.
#' @param maxblocksize The maximum size of the prophage-like block pattern. Default is NA.
#' @keywords internal
block_builder <- function (viral_subset, windowsize, minblocksize, maxblocksize) {
  max_read_cov <- max(viral_subset[,2])
  min_read_cov <- min(viral_subset[,2])
  quarter_read_cov <- abs(max_read_cov-min_read_cov)/4
  bottomtotop_read_cov <- abs(max_read_cov - (min_read_cov+quarter_read_cov))/10
  Cov_values_contig <- viral_subset[,2]
  if(min_read_cov>(max_read_cov*0.2)) {
    min_read_cov <- 0
  }
  startingcoverages <- seq((min_read_cov+quarter_read_cov),  max_read_cov, bottomtotop_read_cov)
  shape_length <- ifelse((nrow(viral_subset)-(10000/windowsize))>(maxblocksize/windowsize),maxblocksize/windowsize,nrow(viral_subset)-(10000/windowsize))
  nonshape <- nrow(viral_subset)-shape_length
  
  #full pattern
  shape_length_full <- ifelse((nrow(viral_subset)-(20000/windowsize))>(maxblocksize/windowsize), maxblocksize/windowsize, nrow(viral_subset)-(20000/windowsize))
  nonshape_full <- nrow(viral_subset)-(shape_length_full+(10000/windowsize))
  pattern_full <- c(rep(min_read_cov, 10000/windowsize), rep(startingcoverages[1], shape_length_full), rep(min_read_cov, nonshape_full))
  diff_full <- mean(abs(Cov_values_contig - pattern_full))
  start_pos_full <- (which(pattern_full == max(pattern_full))[1])
  end_pos_full <- which(pattern_full==max(pattern_full))[length(which(pattern_full==max(pattern_full)))]
  best_match_info_full <- list(diff_full, min_read_cov, startingcoverages[1], "NA", start_pos_full, end_pos_full, "NA")
  
  #right pattern
  pattern_right <- c(rep(min_read_cov, nonshape), rep(startingcoverages[1], shape_length)) 
  diff_right <- mean(abs(Cov_values_contig - pattern_right))
  start_pos_right <- (which(pattern_right == max(pattern_right))[1]) 
  best_match_info_right <- list(diff_right, min_read_cov, startingcoverages[1], "NA", start_pos_right, length(pattern_right), "NA") 
  #left pattern
  pattern_left <- c(rep(startingcoverages[1], shape_length), rep(min_read_cov, nonshape))   #Different in each function
  diff_left <- mean(abs(Cov_values_contig - pattern_left))
  end_pos_left <- (which(pattern_left == min(pattern_left))[1])-1   #The start and stop positions used for determining the pattern length are different in each function
  best_match_info_left <- list(diff_left, min_read_cov, startingcoverages[1], "NA", 1, end_pos_left, "NA")   #Recorded different for each function
  lapply(1:length(startingcoverages), function(i) {
    cov<-startingcoverages[[i]]
    pattern_full <- c(rep(min_read_cov, 10000/windowsize), rep(cov, shape_length_full), rep(min_read_cov, nonshape_full))
    pattern_left <- c(rep(cov, shape_length), rep(min_read_cov, nonshape))
    pattern_right <- c(rep(min_read_cov, nonshape), rep(cov, shape_length))
    indic_left<- 0
    indic_right<-0
    indic_full<-0
    repeat {
      #left repeats
      if(indic_left != 1){
        diff_left <- mean(abs(Cov_values_contig - pattern_left))
        end_pos_left <- (which(pattern_left == min_read_cov)[1])-1 #Different in each function
        if (diff_left < best_match_info_left[[1]]){
          best_match_info_left <<- list(diff_left, min_read_cov, cov, "NA", 1, end_pos_left, "NA") #Different in each function
        }
        pattern_left <- c(pattern_left[-c(1:(2000/windowsize))], rep(min_read_cov, (2000/windowsize))) #Different in each function
        if (length(which(pattern_left==cov)) < (minblocksize/windowsize)+1){
          indic_left <-1
        }  
      }
      #Right repeats
      if(indic_right != 1){
        diff_right <- mean(abs(Cov_values_contig - pattern_right))
        start_pos_right <- (which(pattern_right == max(pattern_right))[1])
        if (diff_right < best_match_info_right[[1]]){
          best_match_info_right <<- list(diff_right, min_read_cov, cov, "NA", start_pos_right, length(pattern_right), "NA") #Different in each function
        }
        pattern_right <- c(rep(min_read_cov,(2000/windowsize)),pattern_right[-c(((length(pattern_right))-((2000/windowsize)-1)):length(pattern_right))]) #Different in each function
        if (length(which(pattern_right==cov)) < (minblocksize/windowsize)+1){
          indic_right <- 1
        }
      }
      if(indic_full !=1){
        middle_rows_full <- which(pattern_full == cov)
        if (length(middle_rows_full) < (minblocksize/windowsize)+1){
          indic_full<-1
        }else{
          best_match_info_full <<- blockpattern_translator(viral_subset, best_match_info_full, windowsize, pattern_full)
          pattern_full <- c(pattern_full[-c(middle_rows_full[2]:middle_rows_full[(2000/windowsize)+1])],rep(min_read_cov,2000/windowsize)) #remove 2000bp at a time
        }
        
      }
      if (indic_right==1 & indic_left==1 & indic_full==1) break
    }
  })
  #left repeat
  new_pattern_coverages_left <- blockheight_optim(best_match_info_left, startingcoverages, bottomtotop_read_cov)
  lapply(1:length(new_pattern_coverages_left), function(i) {
    newcov_left<-new_pattern_coverages_left[[i]]
    pattern_left <- c(rep(newcov_left, shape_length), rep(min_read_cov, nonshape))
    repeat {
      end_pos_left <- (which(pattern_left == min(pattern_left))[1])-1 #Different in each function
      diff_left <- mean(abs(Cov_values_contig - pattern_left))
      if (diff_left < best_match_info_left[[1]]){
        best_match_info_left <<- list(diff_left, min_read_cov, newcov_left, "NA", 1, end_pos_left, "NA") #Different in each function
      }
      pattern_left <- c(pattern_left[-c(1:(2000/windowsize))], rep(min_read_cov, (2000/windowsize))) #Different in each function
      if (length(which(pattern_left==newcov_left)) < (minblocksize/windowsize)+1) break
    }
  })
  best_match_results_left <- c(best_match_info_left, "Prophage-like")
  all_vars <- ls()
  indices_to_keep <- endsWith(all_vars, "_left")
  vars_to_remove<- all_vars[indices_to_keep]
  vars_to_remove <- vars_to_remove[vars_to_remove != "best_match_results_left"]
  rm(list = vars_to_remove)
  #right repeat
  new_pattern_coverages_right <- blockheight_optim(best_match_info_right, startingcoverages, bottomtotop_read_cov)
  lapply(1:length(new_pattern_coverages_right), function(i) {
    newcov_right<-new_pattern_coverages_right[[i]]
    pattern_right <- c(rep(min_read_cov, nonshape), rep(newcov_right, shape_length)) #Different in each function
    repeat {
      diff_right <- mean(abs(Cov_values_contig - pattern_right))
      start_pos_right <- (which(pattern_right == max(pattern_right))[1]) #Different in each function
      if (diff_right < best_match_info_right[[1]]){
        best_match_info_right <<- list(diff_right, min_read_cov, newcov_right, "NA", start_pos_right, length(pattern_right), "NA")
      }
      pattern_right <- c(rep(min_read_cov,(2000/windowsize)),pattern_right[-c(((length(pattern_right))-((2000/windowsize)-1)):length(pattern_right))]) ##Different in each function
      if (length(which(pattern_right==newcov_right)) < (minblocksize/windowsize)+1) break
    }
  })
  best_match_results_right <- c(best_match_info_right, "Prophage-like")
  all_vars <- ls()
  indices_to_keep <- endsWith(all_vars, "_right")
  vars_to_remove<- all_vars[indices_to_keep]
  vars_to_remove <- vars_to_remove[vars_to_remove != "best_match_results_right"]
  rm(list = vars_to_remove)
  #full repeat
  new_pattern_coverages_full <- blockheight_optim(best_match_info_full, startingcoverages, bottomtotop_read_cov)
  lapply(1:length(new_pattern_coverages_full), function(i) {
    newcov_full<-new_pattern_coverages_full[[i]]
    pattern_full <- c(rep(min_read_cov, 10000/windowsize), rep(newcov_full, shape_length_full), rep(min_read_cov, nonshape_full))
    repeat {
      middle_rows_full <- which(pattern_full == newcov_full)
      if (length(middle_rows_full) < (minblocksize/windowsize)+1) break
      best_match_info_full <<- blockpattern_translator(viral_subset, best_match_info_full, windowsize, pattern_full)
      pattern_full <- c(pattern_full[-c(middle_rows_full[2]:middle_rows_full[(2000/windowsize)+1])],rep(min_read_cov,2000/windowsize)) #remove 2000bp at a time
    }
  })
  best_match_results_full <- c(best_match_info_full, "Prophage-like")
  return(list(best_match_results_left, best_match_results_right, best_match_results_full))
}

