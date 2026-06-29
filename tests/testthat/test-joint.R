test_that("jointDVI(example2) reproduces table", {
  res = suppressMessages(jointDVI(example2, verbose = FALSE))

  expected = data.frame(
    V1 = c("M1", "M1", "*", "M1", "*", "*", "*", "M1", "*", "*"),
    V2 = c("M2", "M2", "M2", "*", "M1", "M2", "*", "*", "M1", "*"),
    V3 = c("M3", "*", "M3", "M3", "M3", "*", "M3", "*", "*", "*"))

  expect_equal(res[c("V1", "V2", "V3")], expected)
  expect_equal(round(res$loglik, 5), c(-16.11810, -17.72753, -18.42068, -20.03012,
                                       -20.03012, -20.03012, -20.03012, -21.63956,
                                       -21.63956, -21.63956))
  expect_equal(res$LR, c(250, 50, 25, 5, 5, 5, 5, 1, 1, 1))
  expect_equal(res$posterior, res$LR/sum(res$LR))
})

test_that("Bmarginal gives documented posterior probabilities", {
  joint = suppressMessages(jointDVI(example2, verbose = FALSE))

  expected = matrix(c(
    0.87931034, 0,          0,         0.12068966,
    0.01724138, 0.94827586, 0,         0.03448276,
    0,          0,          0.83333333,0.16666667),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("V1", "V2", "V3"), c("M1", "M2", "M3", "*")))

  expect_equal(Bmarginal(joint, example2$missing), expected, tolerance = 1e-8)

  prior = c(1, rep(0, nrow(joint) - 1))
  expected = matrix(c(
    1,0,0,0,
    0,1,0,0,
    0,0,1,0),
    nrow = 3, byrow = TRUE, dimnames = dimnames(expected))
  
  B = Bmarginal(joint, example2$missing, prior = prior)
  expect_equal(B, expected)
})

test_that("dviJoint ranks example2 and keeps PM/AM orientations consistent", {
  res = dviJoint(example2, verbose = FALSE, progress = FALSE)

  expect_equal(unlist(res[1, c("V1", "V2", "V3")], use.names = FALSE), c("M1", "M2", "M3"))
  expect_equal(unlist(res[1, c("M1", "M2", "M3")], use.names = FALSE), c("V1", "V2", "V3"))
  expect_equal(res$LR[1], 1)

  finite_ll = res$loglik[is.finite(res$loglik)]
  expect_true(all(diff(finite_ll) <= sqrt(.Machine$double.eps)))

  assig = res[c("V1", "V2", "V3")]
  expect_equal(swapOrientation(swapOrientation(assig)), assig)
})
