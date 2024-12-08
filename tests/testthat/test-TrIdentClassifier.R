
test_that("default_TridentClassifer", {
  data("VLPFractionSamplePileup")
  data("WholeCommunitySamplePileup")
  temp_tri_results <- TrIdentClassifier(VLPpileup=VLPFractionSamplePileup,
                                        WCpileup=WholeCommunitySamplePileup)
  expect_equal(temp_tri_results , trident_results_v1)
})
