context("Plot_TrIdentPatternMatches")

test_that("Plot_TrIdentPatternMatches", {
  TrIdent_results <- TrIdent_Classifier(VLP_pileup=VLPFraction_sampledata, WC_pileup=WholeCommunity_sampledata, windowsize=1000, minblocksize=10000, maxblocksize=Inf, cleanup=TRUE)
  temp_patternmatcher <- Plot_TrIdentPatternMatches(VLP_pileup=VLPFraction_sampledata, WC_pileup=WholeCommunity_sampledata, transductionclassifications=TrIdent_results, cleanup=TRUE)
  expect_equal(temp_patternmatcher, TrIdent_patternmatches)
})
