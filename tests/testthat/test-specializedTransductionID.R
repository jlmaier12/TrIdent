test_that("specializedTransductionID", {
  data("VLPFractionSamplePileup")
  data("WholeCommunitySamplePileup")
  TrIdent_results <- TrIdentClassifier(VLPpileup=VLPFractionSamplePileup, WCpileup=WholeCommunitySamplePileup)
  temp_Specialized_transduction <- specializedTransductionID(VLPpileup=VLPFractionSamplePileup, TrIdentResults=TrIdent_results)
  copied_Specialized_transduction <- Specialized_transduction
  temp_plot<-temp_Specialized_transduction$Plots$plot
  copied_plot <- copied_Specialized_transduction$Plots$plot
  expect_equal(temp_Specialized_transduction$Summary_table, copied_Specialized_transduction$Summary_table)
  expect_equal(temp_plot, copied_plot)

  temp_Specialized_transduction_n44 <- specializedTransductionID(VLPpileup=VLPFractionSamplePileup, TrIdentResults=TrIdent_results, specificContig="NODE_44")
  copied_Specialized_transduction_NODE44 <- Specialized_transduction_NODE44
  expect_equal(temp_Specialized_transduction_n44$Plots$plot, copied_Specialized_transduction_NODE44$Plots$plot)
})

#test_that("SpecializedTransduction_ID_Node_44", {
#  TrIdent_results <- TrIdentClassifier(VLP_pileup=VLPFractionSampleData, WC_pileup=WholeCommunitySamplePileup, windowsize=1000, minblocksize=10000, maxblocksize=Inf, cleanup=TRUE)
#  temp_Specialized_transduction_n44 <- SpecializedTransduction_ID(VLP_pileup=VLPFractionSampleData, transductionclassifications=TrIdent_results, specificcontig="NODE_44", noreadcov=500, spectranslength=2000, cleanup=TRUE)
#  copied_Specialized_transduction_NODE44 <- Specialized_transduction_NODE44
#  expect_equal(temp_Specialized_transduction_n44$Plots$plot, copied_Specialized_transduction_NODE44$Plots$plot)
#})
