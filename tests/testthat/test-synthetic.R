# Helper function to create a simple trio DVI dataset with 1 missing and 1 victim.
trio_dvi = function(vgeno = "1/1", p = 0.1, parents = c("1/2", "1/2"), sex = 1) {
  afr = c("1" = p, "2" = 1 - p)
  pm = singleton("V", sex = sex) |> addMarker(geno = vgeno, afreq = afr, name = "m")
  am = nuclearPed(father = "fa", mother = "mo", children = "M", sex = sex) |>
    addMarker(geno = c(parents, NA), afreq = afr, name = "m")
  dviData(pm = pm, am = am, missing = "M")
}


test_that("a one-child trio has the analytically known pairwise LR", {
  dvi = trio_dvi(vgeno = "1/1", p = 0.1)
  lr = pairwiseLR(dvi, verbose = FALSE)$LRmatrix

  expect_equal(lr["V", "M"], 25, tolerance = 1e-12)
})

test_that("dviSolve identifies a simple trio match above threshold", {
  dvi = trio_dvi(vgeno = "1/1", p = 0.1)
  res = dviSolve(dvi, threshold = 10, verbose = FALSE)

  expect_equal(res$AM$Missing, "M")
  expect_equal(res$AM$Sample, "V")
  expect_equal(res$AM$Conclusion, "Undisputed")
  expect_equal(res$AM$LR, 25)
  expect_equal(res$PM$Sample, "V")
  expect_equal(res$PM$Conclusion, "Undisputed")
})

test_that("dviSolve reports no match when the best LR is below one", {
  dvi = trio_dvi(vgeno = "2/2", p = 0.1)
  expected_lr = 1/(4 * 0.9^2)

  expect_equal(pairwiseLR(dvi, verbose = FALSE)$LRmatrix["V", "M"], expected_lr, tolerance = 1e-12)

  res = dviSolve(dvi, threshold = 10, verbose = FALSE)
  expect_equal(res$AM$Conclusion, "No match")
  expect_equal(res$PM$Conclusion, "No match")
  expect_equal(res$AM$LR, expected_lr, tolerance = 1e-12)
})

test_that("Mendelian incompatibility is excluded before the LR workflow", {
  dvi = trio_dvi(vgeno = "2/2", p = 0.5, parents = c("1/1", "1/1"))
  excl = findExcluded(dvi, maxIncomp = 0, verbose = FALSE)

  expect_equal(excl$exclusionMatrix["V", "M"], 1L)
  expect_equal(excl$summary$PM$Sample, "V")
  expect_equal(excl$summary$AM$Missing, "M")

  res = dviSolve(dvi, threshold = 10, maxIncomp = 0, verbose = FALSE)
  expect_equal(res$PM$Conclusion, "Excluded")
  expect_equal(res$AM$Conclusion, "Excluded")
})

test_that("setPairing transfers victim data before reducing away the matched missing person", {
  dvi = trio_dvi(vgeno = "1/1", p = 0.1)
  res = setPairing(dvi, match = c(V = "M"), verbose = FALSE)

  expect_equal(res$summary$Missing, "M")
  expect_equal(res$summary$Sample, "V")
  expect_length(res$dviReduced$pm, 0)
  expect_length(res$dviReduced$missing, 0)
})

test_that("setPairing works when a whole AM family is removed", {
  r = setPairing(example2, match = c(V3 = "M3"), verbose = FALSE)

  expect_equal(
    r$summary,
    data.frame(
      Family = "F2",
      Missing = "M3",
      Sample = "V3",
      Conclusion = "Provided",
      Comment = ""
    )
  )

  expect_setequal(labels(r$dviReduced$pm), c("V1", "V2"))
  expect_equal(r$dviReduced$missing, c("M1", "M2"))
  expect_equal(names(r$dviReduced$am), "F1")
})
