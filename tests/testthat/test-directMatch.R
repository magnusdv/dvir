test_that("directMatch handles dropout", {
  p = 0.1
  q = 1 - p
  d = 0.2
  s = 1 - d
  afr = c("1" = p, "2" = q)

  pm = singletons(c("V1", "V2", "V3", "V4")) |>
    addMarker(V1 = "1/2", V2 = "1/2", V3 = "1/1", V4 = "2/2",
              afreq = afr, name = "M")

  # Equal heterozygotes: dropout cancels from the LR
  expect_equal(directMatch(pm[[1]], pm[[2]]), 1/(2 * p * q))
  expect_equal(directMatch(pm[[1]], pm[[2]], dropout = d), 1/(2 * p * q))

  # Heterozygote vs homozygote: exclusion without dropout
  expect_equal(directMatch(pm[[1]], pm[[3]]), 0)

  likAA = p^2 * (1 - d^2) + 2 * p * q * d * s
  expect_equal(directMatch(pm[[1]], pm[[3]], dropout = d), d * s/likAA)

  # Opposite homozygotes are possible through opposite dropout
  likBB = q^2 * (1 - d^2) + 2 * p * q * d * s
  expect_equal(directMatch(pm[[3]], pm[[4]], dropout = d),
               2 * p * q * (d * s)^2/(likAA * likBB))
})