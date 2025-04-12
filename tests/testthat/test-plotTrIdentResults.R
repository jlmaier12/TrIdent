# test_that("plotTrIdentResults", {
#   data("VLPFractionSamplePileup")
#   data("WholeCommunitySamplePileup")
#   TrIdent_results <- TrIdentClassifier(
#     VLPpileup = VLPFractionSamplePileup,
#     WCpileup = WholeCommunitySamplePileup
#   )
#   temp_patternmatcher <- plotTrIdentResults(
#     VLPpileup = VLPFractionSamplePileup,
#     WCpileup = WholeCommunitySamplePileup,
#     TrIdentResults = TrIdent_results
#   )
#   copied_TrIdent_patternmatches <- TrIdentPlots$plot
#   temp_patternmatcher <- temp_patternmatcher$plot
#   expect_equal(copied_TrIdent_patternmatches, temp_patternmatcher)
# })
