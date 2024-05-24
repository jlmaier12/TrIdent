context("Trident_classifer")

test_that("default_Trident_classifer", {
  temp_tri_results <- TrIdent_Classifier(VLP_pileup=VLPFraction_sampledata, WC_pileup=WholeCommunity_sampledata)
  expect_equal(temp_tri_results , default_trident_results)
})
