#' Build sloping pattern that slopes left-to-right with an initial jump-up in read coverage
#'
#' Build a sloping pattern that consists of a sloping line that starts on the contig being assessed. The pattern is used as input for the slopepattern_translator. The slope of the pattern is changed after each full translation across a contig.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @keywords internal
slope_direct_withstart <- function (viral_subset, windowsize) {
  max_read_cov <- max(viral_subset[,2])
  min_read_cov <- min(viral_subset[,2])
  half_read_cov <- abs((max_read_cov-min_read_cov))/2
  bottomtotop_read_cov <- as.numeric(abs((max_read_cov - (min_read_cov+half_read_cov)))/10)
  newmax <- max_read_cov+bottomtotop_read_cov
  Cov_values_contig <- viral_subset[,2]

  #left to right with start
  cov_steps_lrs <- (newmax-min_read_cov)/((nrow(viral_subset)-((10000/windowsize)+1)))
  pattern_lrs <- c(rep(min_read_cov,10000/windowsize), seq(newmax,min_read_cov, -cov_steps_lrs)) #different in each function
  diff_lrs <- mean(abs(Cov_values_contig - pattern_lrs))
  start_pos_lrs <- which(pattern_lrs==max(pattern_lrs)) #different in each function
  slope_lrs <- (min_read_cov-newmax)/(nrow(viral_subset)-1) #different in each function
  best_match_info_lrs <- list(diff_lrs, min_read_cov, newmax, -cov_steps_lrs, start_pos_lrs, length(pattern_lrs), slope_lrs) #different in each function
  
  #right to left with start
  cov_steps_rls <- (newmax-min_read_cov)/((nrow(viral_subset)-((10000/windowsize)+1)))
  pattern_rls <- c(seq(min_read_cov,newmax,cov_steps_rls),rep(min_read_cov,10000/windowsize)) #different in each function
  diff_rls <- mean(abs(Cov_values_contig - pattern_rls))
  end_pos_rls<- which(pattern_rls==max(pattern_rls)) #different in each function
  slope_rls <- (newmax-min_read_cov)/(nrow(viral_subset)-1) #different in each function
  best_match_info_rls <- list(diff_rls, min_read_cov, newmax, cov_steps_rls, 1, end_pos_rls, slope_rls) #different in each function

  
  lapply(seq(newmax, (min_read_cov+half_read_cov), -bottomtotop_read_cov), function(cov) {
    #indic vars
    lrs_done_indic<-0
    rls_done_indic<-0
    
    #left to right with start
    slope_bottom_lrs <- min_read_cov
    cov_steps_lrs <- (cov-slope_bottom_lrs)/((nrow(viral_subset)-((10000/windowsize)+1)))
    pattern_lrs <- c(rep(min_read_cov,10000/windowsize), seq(cov,slope_bottom_lrs, -cov_steps_lrs)) #different in each function
    slope_lrs <- (cov-slope_bottom_lrs)/(nrow(viral_subset)-((10000/windowsize)+1)) #different in each function
    step_lrs <- ((slope_bottom_lrs-cov)/10) #different in each function
    
    #right to left with start
    slope_bottom_rls <- min_read_cov
    cov_steps_rls <- (cov-slope_bottom_rls)/((nrow(viral_subset)-((10000/windowsize)+1)))
    pattern_rls <- c(seq(slope_bottom_rls,cov,cov_steps_rls),rep(min_read_cov,10000/windowsize)) #different in each function
    slope_rls <- (cov-slope_bottom_rls)/(nrow(viral_subset)-((10000/windowsize)+1)) #different in each function
    step_rls <- ((cov-slope_bottom_rls)/10) #different in each function
    
    #indic updater block
    if (abs(slope_lrs) < (15/100000) | slope_lrs >0){
      lrs_done_indic<-1
    } 
    if (abs(slope_rls) < (15/100000) | slope_rls <0){
      rls_done_indic<-1
    }
    
    if(lrs_done_indic != 1){
      repeat {
        best_match_info_lrs <- slopepattern_translator(viral_subset, best_match_info_lrs, windowsize, pattern_lrs, "lefttoright") #different in each function
        slope_bottom_lrs <- slope_bottom_lrs + step_lrs
        cov_steps_lrs <- (cov-slope_bottom_lrs)/((nrow(viral_subset)-((10000/windowsize)+1)))
        pattern_lrs <- c(rep(min_read_cov,10000/windowsize), seq(cov,slope_bottom_lrs, -cov_steps_lrs)) #different in each function
        slope_lrs <- (slope_bottom_lrs-cov)/(nrow(viral_subset)-((10000/windowsize)+1)) #different in each function
        if (abs(slope_lrs) < (15/100000) | slope_lrs >0) break #different in each function
      }
    }
    if(rls_done_indic != 1){
      repeat {
        best_match_info_rls <<- slopepattern_translator(viral_subset,best_match_info_rls, windowsize, pattern_rls, "righttoleft") #different in each function
        slope_bottom_rls <- slope_bottom_rls + step_rls
        cov_steps_rls <- (cov-slope_bottom_rls)/((nrow(viral_subset)-((10000/windowsize)+1)))
        pattern_rls <- c(seq(slope_bottom_rls,cov,cov_steps_rls),rep(min_read_cov,10000/windowsize)) #different in each function
        slope_rls <- (cov-slope_bottom_rls)/(nrow(viral_subset)-((10000/windowsize)+1)) #different in each function
        if (abs(slope_rls) < (15/100000) | slope_rls <0) break #different in each function
      }
    }
    
  })

  return(list(c(best_match_info_lrs, "Sloping"),  c(best_match_info_rls, "Sloping")))
}
