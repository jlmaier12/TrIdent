# test_that("specializedTransductionID", {
#   data("VLPFractionSamplePileup")
#   data("WholeCommunitySamplePileup")
#   TrIdent_results <- TrIdentClassifier(
#     VLPpileup = VLPFractionSamplePileup,
#     WCpileup = WholeCommunitySamplePileup
#   )
#   temp_Specialized_transduction <- specializedTransductionID(
#     VLPpileup = VLPFractionSamplePileup,
#     TrIdentResults = TrIdent_results
#   )
#   copied_Specialized_transduction <- SpecTransduc
#   temp_plot <- temp_Specialized_transduction$Plots$plot
#   copied_plot <- copied_Specialized_transduction$Plots$plot
#   expect_equal(
#     temp_Specialized_transduction$Summary_table,
#     copied_Specialized_transduction$Summary_table
#   )
#   expect_equal(temp_plot, copied_plot)
# 
#   temp_Specialized_transduction_n62 <- specializedTransductionID(
#     VLPpileup = VLPFractionSamplePileup,
#     TrIdentResults = TrIdent_results,
#     specificContig = "NODE_62"
#   )
#   copied_Specialized_transduction_NODE62 <- SpecTransduc_node62
#   expect_equal(
#     temp_Specialized_transduction_n62$Plots$plot,
#     copied_Specialized_transduction_NODE62$Plots$plot
#   )
# })
