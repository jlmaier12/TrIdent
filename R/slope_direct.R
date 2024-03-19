#' Sloping pattern dir pattern builder
#'
#' Build a sloping pattern that consists of a sloping line spanning the contig being assessed. The line slopes from left to right. The slope of the line is changed, but the pattern is not translated across the contig.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @keywords internal
slope_direct <- function (viral_subset, windowsize) {
  max_read_cov <- max(viral_subset[,2])
  min_read_cov <- min(viral_subset[,2])
  half_read_cov <- abs(max_read_cov-min_read_cov)/2
  bottomtotop_read_cov <- abs(max_read_cov - (min_read_cov+half_read_cov))/10
  newmax <- max_read_cov+bottomtotop_read_cov
  bottomtotop_read_cov <- abs((newmax-(min_read_cov+half_read_cov))/10)
  Cov_values_contig <- viral_subset[,2]
  
  #left to right
  cov_steps_lr <- (newmax-min_read_cov)/((nrow(viral_subset)-1))
  pattern_lr <- seq(newmax,min_read_cov, -cov_steps_lr) #different in each function
  diff_lr <- mean(abs(Cov_values_contig - pattern_lr))
  slope_lr <- (min_read_cov-newmax)/(nrow(viral_subset)-1) #different in each function
  best_match_info_lr <- list(diff_lr, min_read_cov, newmax, -cov_steps_lr, 1, length(pattern_lr), slope_lr) #different in each function
  
  #right to left
  cov_steps_rl <- (newmax-min_read_cov)/((nrow(viral_subset)-1))
  pattern_rl <- seq(min_read_cov,newmax,cov_steps_rl) #different in each function
  diff_rl <- mean(abs(Cov_values_contig - pattern_rl))
  slope_rl <- (newmax-min_read_cov)/(nrow(viral_subset)-1) #different in each function
  best_match_info_rl <- list(diff_rl, min_read_cov, newmax, cov_steps_rl, 1, length(pattern_rl), slope_rl) #different in each function
  
  lapply(seq(newmax, (min_read_cov+half_read_cov), -bottomtotop_read_cov), function(cov) {
    lr_done_indic<-0
    rl_done_indic<-0
    
    #left to right
    slope_bottom_lr <- min_read_cov
    cov_steps_lr <- (cov-slope_bottom_lr)/((nrow(viral_subset)-1))
    pattern_lr <- seq(cov,slope_bottom_lr, -cov_steps_lr) #different in each function
    slope_lr <- (slope_bottom_lr-cov)/(nrow(viral_subset)-1) #different in each function
    
    #right to left
    slope_bottom_rl <- min_read_cov
    cov_steps_rl <- (cov-slope_bottom_rl)/((nrow(viral_subset)-1))
    pattern_rl <- seq(slope_bottom_rl,cov,cov_steps_rl) #different in each function
    slope_rl <- (cov-slope_bottom_rl)/(nrow(viral_subset)-1) #different in each function
    
    if (abs(slope_lr) < 15/100000 | slope_lr > 0){
      lr_done_indinc<-1
    }else{
      repeat {
        if (diff_lr < best_match_info_lr[[1]]) {
          best_match_info_lr <<- list(diff_lr, slope_bottom_lr, cov, -cov_steps_lr, 1, length(pattern_lr), slope_lr) #different in each function
        }
        slope_bottom_lr <- slope_bottom_lr + bottomtotop_read_cov
        cov_steps_lr <- (cov-slope_bottom_lr)/((nrow(viral_subset)-1))
        pattern_lr <- seq(cov,slope_bottom_lr, -cov_steps_lr) #different in each function
        diff_lr <- mean(abs(Cov_values_contig - pattern_lr))
        slope_lr <- (slope_bottom_lr-cov)/(nrow(viral_subset)-1) #different in each function
        if (abs(slope_lr) < 15/100000 | slope_lr >0) break
      }
    } 
    if (abs(slope_rl) < 15/100000 | slope_rl < 0){
      rl_done_indic<-1
    }else{
      repeat {
        if (diff_rl < best_match_info_rl[[1]]) {
          best_match_info_rl <<- list(diff_rl, slope_bottom_rl, cov, cov_steps_rl, 1, length(pattern_rl), slope_rl) #different in each function
        }
        slope_bottom_rl <- slope_bottom_rl + bottomtotop_read_cov
        cov_steps_rl <- (cov-slope_bottom_rl)/((nrow(viral_subset)-1))
        pattern_rl <- seq(slope_bottom_rl,cov,cov_steps_rl) #different in each function
        diff_rl <- mean(abs(Cov_values_contig - pattern_rl))
        slope_rl <- (cov-slope_bottom_rl)/(nrow(viral_subset)-1) #different in each function
        if (abs(slope_rl) < 15/100000 | slope_rl < 0) break #different in each function
      }
    }
  })
  return(list(c(best_match_info_lr, "Gen/Lat/GTA"), c(best_match_info_rl, "Gen/Lat/GTA")))
}
