data("polyA_small")

pct <- percentage.usage(polyA_small, features = rownames(polyA_small),
                        cells = head(colnames(polyA_small), 20))
names(pct) <- rownames(polyA_small)
test_that("percentages within same gene sum to 1", {
  expect_equal(as.numeric(pct["11-118408896-118409195"] + pct["11-118409548-118409847"]), 1)
})
