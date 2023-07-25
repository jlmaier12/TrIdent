#' .gff3 file re-formatter
#'
#' Gene annotations cleanup
#'
#' @param gff_file A gff3 file containing gene annotation information
#'
#' @export
#' @keywords internal
gff3_cleanup <- function(gff_file) {
  V9 <- NULL
  gff_df <- gff_file[-c(2,3,6,7,8)]
  colnames(gff_df)[1] <- "ref_name"
  gff_df <- separate(gff_df, col=V9, into=c('ID', 'annotations'), sep=';', extra="merge")
  inference <- str_extract(gff_df[,5], "(?<=inference=)[\\s\\S]*") %>% str_extract(".+?(?=;)")
  EC_number <- str_extract(gff_df[,5], "(?<=eC_number=)[\\s\\S]*") %>% str_extract(".+?(?=;)")
  gene <- str_extract(gff_df[,5], "(?<=gene=)[\\s\\S]*") %>% str_extract(".+?(?=;)")
  product <- str_extract(gff_df[,5], "(?<=;product=)[\\s\\S]*")
  cleaned_gff_df <- cbind.data.frame(gff_df[,c(1:4)], inference, EC_number, gene, product)
  cleaned_gff_df[,1] <- gsub("_length.*", "", gff_df[,1])
  cleaned_gff_df <- cleaned_gff_df[which(str_starts(cleaned_gff_df[,1],"NODE.*")==TRUE),]
  return(cleaned_gff_df)
}
