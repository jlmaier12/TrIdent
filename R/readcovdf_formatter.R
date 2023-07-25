#' Correctly formats input read coverage summary files.
#'
#' Places columns in correct order and renames columns. Clean the contig labels to remove excess informatio.
#'
#' @param read_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions
#' @keywords internal
readcovdf_formatter <- function(read_dataset) {
  column_classes <- c()
  for (i in c(1:ncol(read_dataset))) {
    column_classes <- c(column_classes, class(read_dataset[1,i]))
  }
  reformatted_readdataset <- cbind.data.frame(read_dataset[,which(column_classes == "character")],read_dataset[,which(column_classes == "numeric")], read_dataset[,which(column_classes == "integer")[1]])
  colnames(reformatted_readdataset) <- c("ref_name", "coverage", "position")
  reformatted_readdataset$ref_name <- gsub(" ", "", reformatted_readdataset$ref_name)
  reformatted_readdataset$ref_name <- gsub("length_.*", "", as.factor(reformatted_readdataset$ref_name))
  return(reformatted_readdataset)
}
