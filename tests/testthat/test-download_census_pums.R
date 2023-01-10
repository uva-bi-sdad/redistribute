test_that("files are downloaded and loaded", {
  dir <- "../../census_download_tests"
  if (grepl("R_LIBS", getwd(), fixed = TRUE)) dir.create(dir, FALSE, TRUE)
  if (!dir.exists(dir)) dir <- "../../../../census_download_tests"
  skip_if_not(dir.exists(dir), "not downloading data")
  data <- download_census_pums(dir, "ak", calculate_error = TRUE)
  expect_identical(names(data), c(
    "year", "state", "dictionary", "household", "person",
    "household_error", "person_error", "crosswalk"
  ))
  subset <- download_census_pums(dir, "ak", geoids = c("02013", "02020", "02090000100"))
  expect_identical(unique(subset$household$PUMA), c("00101", "00102", "00400", "00300"))
})
