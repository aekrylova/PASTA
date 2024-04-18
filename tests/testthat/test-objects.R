
data("polyA_small")

test_that("Subsetting cells works", {
    x <- subset(polyA_small, cells = colnames(polyA_small)[1:10])
    expect_equal(Cells(x), colnames(x = polyA_small)[1:10])
})

test_that("Subsetting features works", {
  x <- subset(polyA_small, features = rownames(polyA_small)[1:10])
  expect_equal(rownames(x), rownames(x = polyA_small)[1:10])
})


test_that("Merging objects works", {
  obj.list <- SplitObject(polyA_small, split.by = "group")
  m <- merge(obj.list[[1]], obj.list[[2]])
  expect_equal(ncol(m), ncol(polyA_small))
})
