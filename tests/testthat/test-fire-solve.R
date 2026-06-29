test_that("dviSolve handles fire example with detailed output", {
  r = dviSolve(fire, detailedOutput = TRUE, verbose = FALSE)

  expect_named(r, c("AM", "PM", "LRmatrix", "exclusionMatrix", "jointTables"))

  expect_equal(
    r$AM,
    data.frame(
      Family = c("F1", "F1", "F1"),
      Missing = c("M1", "M2", "M3"),
      Sample = c("V1", "V2/V3", "V2/V3"),
      LR = c(121.347585, NA, NA),
      GLR = c(70582.39, 2947149.10, 2947149.10),
      Conclusion = c("Match (GLR)", "Symmetric match", "Symmetric match"),
      Comment = c("Joint: {M1,M2,M3}",
                  "Full siblings: {M2, M3} = {V2, V3}",
                  "Full siblings: {M2, M3} = {V2, V3}")
    ),
    tolerance = 1e-4
  )

  expect_equal(
    r$PM,
    data.frame(
      Sample = c("V1", "V2", "V3"),
      Missing = c("M1", "M2/M3", "M2/M3"),
      Family = c("F1", "F1", "F1"),
      LR = c(121.347585, 36.563503, 1.312181),
      GLR = c(70582.39, 2947149.10, 2947149.10),
      Conclusion = c("Match (GLR)", "Symmetric match", "Symmetric match"),
      Comment = c("Joint: {M1,M2,M3}",
                  "Full siblings: {V2, V3} = {M2, M3}",
                  "Full siblings: {V2, V3} = {M2, M3}")
    ),
    tolerance = 1e-4
  )

  expect_equal(
    r$LRmatrix,
    matrix(
      c(121.347585, 0.5826857, 0.002019597,
        398.582218, 36.563503, 1.312181,
        398.582218, 36.563503, 1.312181),
      nrow = 3,
      dimnames = list(c("V1", "V2", "V3"), c("M1", "M2", "M3"))
    ),
    tolerance = 1e-6
  )

  expect_equal(
    r$exclusionMatrix,
    matrix(0, nrow = 3, ncol = 3,
           dimnames = list(c("V1", "V2", "V3"), c("M1", "M2", "M3")))
  )

  j = r$jointTables$F1
  expect_equal(nrow(j), 34)
  expect_equal(names(j), c("V1", "V2", "V3", "loglik", "LR", "M1", "M2", "M3"))

  expect_equal(
    j[1:6, ],
    data.frame(
      V1 = c("M1", "M1", "*", "*", "M2", "M3"),
      V2 = c("M2", "M3", "M2", "M3", "M1", "M1"),
      V3 = c("M3", "M2", "M3", "M2", "M3", "M2"),
      loglik = c(-257.7327, -257.7327, -268.8973, -268.8973, -272.6291, -272.6291),
      LR = c(1, 1, 70582.39, 70582.39, 2947149.10, 2947149.10),
      M1 = c("V1", "V1", "*", "*", "V2", "V2"),
      M2 = c("V2", "V3", "V2", "V3", "V1", "V3"),
      M3 = c("V3", "V2", "V3", "V2", "V3", "V1")
    ),
    tolerance = 1e-4
  )
})
