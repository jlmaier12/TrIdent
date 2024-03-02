context("Trident_classifer")

test_that("default_Trident_classifer", {
  #devtools::load_all()
  temp_tri_results <- TrIdent_Classifier(VLP_pileup=VLPFraction_sampledata, WC_pileup=WholeCommunity_sampledata, windowsize=1000, minblocksize=10000, maxblocksize=Inf, cleanup=TRUE)
  expect_equal(temp_tri_results , TrIdent_results)
})
