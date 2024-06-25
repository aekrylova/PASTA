data("polyA_small")

polyA_small <- CalcPolyAResiduals(polyA_small)

test_that("Residuals are 0 when cells have 0 counts for that gene", {
  cells.0 <- colnames(polyA_small)[polyA_small[['polyA']]@counts["11-118409548-118409847",] + polyA_small[['polyA']]@counts["11-118408896-118409195",] == 0]
  expect_equal(sum(polyA_small[['polyA']]@scale.data["11-118409548-118409847",cells.0]), 0)
})

