test_that("generatePairings respects sex unless ignoreSex is TRUE", {
  p1 = generatePairings(example2)
  p2 = generatePairings(example2, ignoreSex = TRUE)

  expect_setequal(p1$V1, c("*", "M1", "M2"))
  expect_setequal(p1$V2, c("*", "M1", "M2"))
  expect_setequal(p1$V3, c("*", "M3"))
  expect_setequal(p2$V1, c("*", "M1", "M2", "M3"))
  expect_setequal(p2$V3, c("*", "M1", "M2", "M3"))
})

test_that("subsetDVI removes stale missing persons from pairings", {
  x = subsetDVI(example2, pm = c("V1", "V3"), missing = c("M1", "M3"), verbose = FALSE)

  expect_equal(names(x$pm), c("V1", "V3"))
  expect_equal(x$missing, c("M1", "M3"))
  expect_false("M2" %in% unlist(x$pairings, use.names = FALSE))
})

test_that("excludePairing affects only the selected victim and is used in dviJoint", {
  x = excludePairing(example2, victim = "V1", missing = "M1")

  expect_false("M1" %in% x$pairings$V1)
  expect_equal(x$pairings$V2, example2$pairings$V2)

  joint = dviJoint(x, verbose = FALSE, progress = FALSE)
  expect_false(any(joint$V1 == "M1"))
})

test_that("dviSolve detailed output stores first-stage matrices", {
  res = dviSolve(example2, threshold = 5, detailedOutput = TRUE, verbose = FALSE)

  expect_true(all(c("AM", "PM", "LRmatrix", "exclusionMatrix", "jointTables") %in% names(res)))
  expect_equal(res$LRmatrix, pairwiseLR(example2, verbose = FALSE)$LRmatrix)
  expect_equal(res$exclusionMatrix, exclusionMatrix(example2))
  expect_length(res$jointTables, 1)
})
