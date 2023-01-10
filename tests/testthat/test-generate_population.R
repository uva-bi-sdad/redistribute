test_that("loci misspecifications and verbose works", {
  expect_identical(
    substring(capture.output(generate_population(
      2,
      attraction_loci = 3, random_regions = FALSE,
      cost_loci = 3, size_loci = 3, verbose = TRUE
    ), type = "message"), 3),
    c(
      "regions: sequence along `N`, with housing unit capacities",
      "preparing 2 households",
      "too many attraction loci; setting `attraction_loci` to 1",
      "calculating region similarities",
      "rescaling region similarities",
      "too many cost loci; setting `cost_loci` to 1",
      "too many size loci; setting `size_loci` to 1",
      "drawing renting status",
      "drawing race baserates",
      "generating individuals"
    )
  )
})

skip_if(!grepl("R_LIBS", getwd(), fixed = TRUE), "not generating larger population")

test_that("IDs and regions are as requested", {
  n_households <- 1e4
  regions <- seq_len(n_households)
  pop <- generate_population(n_households)
  expect_identical(pop$households$household, seq_len(n_households))
  expect_true(all(pop$households$region %in% seq_len(n_households)))
  expect_identical(
    as.integer(pop$individuals$household),
    rep(pop$households$household, pop$households$size)
  )
})
