test_that("example works", {
  dir <- "../../census_download_tests"
  if (grepl("R_LIBS", getwd(), fixed = TRUE)) dir.create(dir, FALSE, TRUE)
  if (!dir.exists(dir)) dir <- "../../../../census_download_tests"
  skip_if_not(dir.exists(dir), "not downloading data")
  options(tigris_use_cache = TRUE)
  tracts <- tidycensus::get_acs(
    year = 2021,
    state = "51",
    county = "013",
    geography = "tract",
    output = "wide",
    variables = c(total = "B01001_001"),
    geometry = TRUE,
    cache_table = TRUE
  )
  parcel_file <- paste0(dir, "/parcels.rds")
  if (!file.exists(parcel_file)) {
    parcels <- sf::st_read(paste0(
      "https://arlgis.arlingtonva.us/arcgis/rest/",
      "services/Open_Data/od_MHUD_Polygons/",
      "FeatureServer/0/query?where=1=1&outFields=*&outSR=4326&f=json"
    ))
    saveRDS(parcels, parcel_file, compress = "xz")
  }
  parcels <- readRDS(parcel_file)
  map_tr_to_parcel <- redistribute(
    tracts, parcels,
    target_id = "OBJECTID", return_map = TRUE
  )
  pums <- download_census_pums(dir, "51", geoids = tracts$GEOID)

  wrapper <- redistribute_parcel_pums_adj(
    tracts[, -2], parcels, pums,
    map = map_tr_to_parcel, target_id = "OBJECTID"
  )

  residents <- parcels$Total_Units * tapply(
    table(pums$person$SERIALNO)[pums$household$SERIALNO],
    pums$household$BLD %in% c("02", "03"),
    mean,
    na.rm = TRUE
  )[as.character(parcels$Unit_Type != "MULTI")]
  manual <- redistribute(
    tracts[, -2], parcels, map_tr_to_parcel,
    target_id = "OBJECTID", weight = residents
  )

  expect_identical(wrapper$totalE, manual$totalE)
})
