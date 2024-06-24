data("polyA_small")

pct <- percentage.usage(polyA_small, features = rownames(polyA_small),
                        cells = head(colnames(polyA_small), 20))
names(pct) <- rownames(polyA_small)
test_that("percentages within same gene sum to 1", {
  expect_equal(as.numeric(pct["11-118408896-118409195"] + pct["11-118409548-118409847"]), 1)
})

test_that("FindDifferentialPolyA throws error if no residuals present", {
expect_error(FindDifferentialPolyA(polyA_small, ident.1 =  "A"),
             regexp = "Run CalcPolyAResiduals prior to FindPolyASites")
})


polyA_small <- CalcPolyAResiduals(polyA_small)
Idents(polyA_small) <- polyA_small$group
m <- FindDifferentialPolyA(polyA_small, ident.1 =  "A", ident.2="B")

test_that("Percentages within a gene sum to 1", {
  expect_equal(sum(subset(m, symbol=="UBC")$percent.1), 1)
  expect_equal(sum(subset(m, symbol=="UBC")$percent.2), 1)
  expect_equal(sum(subset(m, symbol=="ARPC1B")$percent.1), 1)
  expect_equal(sum(subset(m, symbol=="ARPC1B")$percent.2), 1)
})
