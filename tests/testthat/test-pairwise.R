test_that("pairwiseLR(example2) gives known one-marker LRs", {
  res = pairwiseLR(example2, verbose = FALSE)

  expected = matrix(c(
    1, 0, NA,
    1, 5, NA,
    NA, NA, 5),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("V1", "V2", "V3"), c("M1", "M2", "M3")))

  expect_equal(res$LRmatrix, expected)
})

test_that("pairwiseLR limit reduction keeps only supported pairings plus star", {
  p = pairwiseLR(example2, limit = 2, verbose = FALSE)$pairings

  expect_setequal(p$V1, "*")
  expect_setequal(p$V2, c("*", "M2"))
  expect_setequal(p$V3, c("*", "M3"))
})

test_that("findUndisputed handles threshold-dependent cascading matches", {
  res5 = findUndisputed(example2, threshold = 5, verbose = FALSE)

  expect_equal(
    res5$summary,
    data.frame(
      Family = "F2",
      Missing = "M3",
      Sample = "V3",
      LR = 5,
      Conclusion = "Undisputed",
      Comment = "Step 1"
    )
  )

  res2 = findUndisputed(example2, threshold = 2, verbose = FALSE)

  expect_equal(
    res2$summary,
    data.frame(
      Family = c("F1", "F2", "F1"),
      Missing = c("M2", "M3", "M1"),
      Sample = c("V2", "V3", "V1"),
      LR = c(5, 5, 10),
      Conclusion = "Undisputed",
      Comment = c("Step 1", "Step 1", "Step 2")
    )
  )
})


test_that("findExcluded removes incompatible pairings without dropping compatible individuals", {
  res = findExcluded(example2, maxIncomp = 0, verbose = FALSE)

  expect_equal(res$exclusionMatrix["V1", "M2"], 1)
  expect_false("M2" %in% res$dviReduced$pairings$V1)
  expect_equal(res$excluded$sample, character())
  expect_equal(res$excluded$missing, character())
})
