context("SpecializedTransduction_ID")

test_that("SpecializedTransduction_ID", {
  TrIdent_results <- TrIdent_Classifier(VLP_pileup=VLPFraction_sampledata, WC_pileup=WholeCommunity_sampledata, windowsize=1000, minblocksize=10000, maxblocksize=Inf, cleanup=TRUE)
  temp_Specialized_transduction <- SpecializedTransduction_ID(VLP_pileup=VLPFraction_sampledata, transductionclassifications=TrIdent_results, noreadcov=500, spectranslength=2000, cleanup=TRUE)
  

  temp_Specialized_transduction$Plots<-NULL
  copied_Specialized_transduction <-Specialized_transduction
  copied_Specialized_transduction$Plots<-NULL
  expect_equal(temp_Specialized_transduction, copied_Specialized_transduction)
})

test_that("SpecializedTransduction_ID_Node_44", {
  TrIdent_results <- TrIdent_Classifier(VLP_pileup=VLPFraction_sampledata, WC_pileup=WholeCommunity_sampledata, windowsize=1000, minblocksize=10000, maxblocksize=Inf, cleanup=TRUE)
  temp_Specialized_transduction_n44 <- SpecializedTransduction_ID(VLP_pileup=VLPFraction_sampledata, transductionclassifications=TrIdent_results, specificcontig="NODE_44", noreadcov=500, spectranslength=2000, cleanup=TRUE)
  temp_Specialized_transduction_n44<-NULL
  copied_Specialized_transduction_NODE44 <- Specialized_transduction_NODE44
  copied_Specialized_transduction_NODE44<-NULL
  
  expect_equal(temp_Specialized_transduction_n44, copied_Specialized_transduction_NODE44)
})
