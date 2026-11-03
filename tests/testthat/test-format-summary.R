test_that("fixed identifications and unfilled rows share the same ordering", {
  dvi = setPairing(example2, c(V1 = "M1"), prune = FALSE, verbose = FALSE)$dviReduced
  dvi$missing = rev(dvi$missing)
  dvi$pm = rev(dvi$pm)
  am = formatSummary(list(dvi$identifications), "AM", dvi = dvi)
  pm = formatSummary(list(dvi$identifications), "PM", dvi = dvi)

  expect_equal(am$Family, c("F1", "F1", "F2"))
  expect_equal(am$Missing, c("M1", "M2", "M3"))
  expect_equal(pm$Sample, c("V1", "V2", "V3"))
  expect_equal(am$Sample, c("V1", NA, NA))
  expect_equal(pm$Missing, c("M1", NA, NA))
})

test_that("identifications from a removed family are not placed last", {
  dvi = setPairing(example2, c(V2 = "M2", V1 = "M1"), verbose = FALSE)$dviReduced
  am = formatSummary(list(dvi$identifications), "AM", dvi = dvi)
  pm = formatSummary(list(dvi$identifications), "PM", dvi = dvi)

  expect_equal(names(dvi$am), "F2")
  expect_equal(am$Family, c("F1", "F1", "F2"))
  expect_equal(am$Missing, c("M1", "M2", "M3"))
  expect_equal(pm$Sample, c("V1", "V2", "V3"))
})

test_that("dviSolve combines fixed and computed identifications in order", {
  dvi = setPairing(example2, c(V1 = "M1"), prune = FALSE, verbose = FALSE)$dviReduced
  res = dviSolve(dvi, threshold = 2, verbose = FALSE)

  expect_equal(res$AM$Family, c("F1", "F1", "F2"))
  expect_equal(res$AM$Missing, c("M1", "M2", "M3"))
  expect_equal(res$PM$Sample, c("V1", "V2", "V3"))
  expect_equal(which(res$AM$Conclusion == "Provided"), 1L)
  expect_equal(which(res$PM$Conclusion == "Provided"), 1L)
})

test_that("dviSolve retains all identifications when no active samples remain", {
  dvi = setPairing(example2, c(V3 = "M3", V2 = "M2", V1 = "M1"),
                   prune = FALSE, verbose = FALSE)$dviReduced
  res = dviSolve(dvi, detailedOutput = TRUE, verbose = FALSE)

  expect_equal(res$AM$Missing, c("M1", "M2", "M3"))
  expect_equal(res$PM$Sample, c("V1", "V2", "V3"))
  expect_equal(res$AM$Conclusion, rep("Provided", 3))
  expect_equal(res$PM$Conclusion, rep("Provided", 3))
  expect_null(res$LRmatrix)
  expect_length(res$jointTables, 0)
})
