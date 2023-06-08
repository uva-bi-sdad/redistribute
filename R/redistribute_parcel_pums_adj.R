#' Method: Proportional, resident-estimate-adjusted parcel unit counts
#'
#' A wrapper around the \code{\link{redistribute}} function, that
#' takes PUMS data as additional input, and uses it to calculates an adjusted weight for
#' redistribution to parcel data.
#'
#' It is assumed that initial weights are unit counts, and these are to be adjusted
#' based on PUMS household data.
#'
#' @param source,target Source and target of \code{\link{redistribute}}.
#' @param households PUMS household data. Can be a list with \code{household} and \code{person}
#' entries, like that returned from \code{\link{download_census_pums}}.
#' @param target_total Column name in \code{target} (e.g., parcel data; or equivalent vector)
#' containing the unit count
#' @param target_indicator Column name of a logical variable in \code{target} (or equivalent vector)
#' to use to assign adjustment values, where \code{TRUE} indicates a single family home.
#' If the values of the column are characters, anything not matching \code{"MULTI"} will be \code{TRUE}.
#' @param households_indicator Same a \code{target_indicator}, but in \code{households}, or
#' a vector (or column name) of BLD codes, where either \code{"02"} or \code{"03"} would be \code{TRUE}.
#' @param households_size A vector of household sizes, or a column in \code{households} containing
#' such values. If this is not specified or found, \code{person} will be used to calculate
#' household size, based on \code{hosuehold_id} and \code{person_household_id}.
#' @param households_id A vector of household IDs aligning with \code{households},
#' or a column name in \code{households} containing such IDs.
#' Only used to calculate \code{households_size} if necessary.
#' @param person PUMS person data. Only used to calculate \code{households_size} if necessary.
#' @param person_household_id A vector of household IDs aligning with \code{person},
#' or a column name in \code{person} containing such IDs.
#' Only used to calculate \code{households_size} if necessary.
#' @param ... Additional arguments to pass to \code{\link{redistribute}}. You will likely need
#' to specify a \code{map}, and potentially \code{source_id} and/or \code{target_id}.
#' @examples
#' \dontrun{
#' if (require("tidycensus")) {
#'   # download source, target, and household data
#'   tracts <- tidycensus::get_acs(
#'     year = 2021,
#'     state = "51",
#'     county = "013",
#'     geography = "tract",
#'     output = "wide",
#'     variables = c(total = "B01001_001"),
#'     geometry = TRUE
#'   )
#'   parcels <- sf::st_read(paste0(
#'     "https://arlgis.arlingtonva.us/arcgis/rest/",
#'     "services/Open_Data/od_MHUD_Polygons/",
#'     "FeatureServer/0/query?where=1=1&outFields=*&outSR=4326&f=json"
#'   ))
#'   pums <- download_census_pums(tempdir(), "51", geoids = tracts$GEOID)
#'
#'   # calculate map between source and target
#'   map_tr_to_parcel <- redistribute(
#'     tracts, parcels,
#'     target_id = "OBJECTID", return_map = TRUE
#'   )
#'
#'   # redistribute tract-level summaries to parcels,
#'   # using a resident-estimate-adjusted weight
#'   parcels_filled <- redistribute_parcel_pums_adj(
#'     tracts[, -2], parcels, pums,
#'     target_total = "Total_Units", map = map_tr_to_parcel, target_id = "OBJECTID"
#'   )
#'
#'   # this can also be calculated with the underlying function
#'
#'   ## calculate resident estimates
#'   household_size <- tapply(
#'     table(pums$person$SERIALNO)[pums$household$SERIALNO],
#'     pums$household$BLD %in% c("02", "03"),
#'     mean,
#'     na.rm = TRUE
#'   )
#'   residents <- parcels$Total * household_size[
#'     as.character(parcels$Unit_Type != "MULTI")
#'   ]
#'
#'   ## supply as weights
#'   parcels_filled_manual <- redistribute(
#'     tracts[, -2], parcels, map_tr_to_parcel,
#'     target_id = "OBJECTID", weight = residents
#'   )
#'
#'   # either way, you could now use the result to make estimates
#'   # at higher-level geographies by redistributing from the
#'   # parcel-level estimates to the new target
#' }
#' }
#' @returns Result of \code{\link{redistribute}}. These are assumed to be parcel-level
#' estimates, which could then be aggregated to make higher-level estimates.
#' @export

redistribute_parcel_pums_adj <- function(
    source, target, households, target_total = "Total_Units", target_indicator = "Unit_Type",
    households_size = "size", households_indicator = "BLD", households_id = "SERIALNO",
    person = NULL, person_household_id = "SERIALNO", ...) {
  if (is.null(dim(households)) && !is.null(households$household)) {
    if (is.null(person)) person <- households$person
    households <- households$household
  }
  hn <- nrow(households)
  hcols <- colnames(households)
  tn <- nrow(target)
  tcols <- colnames(target)

  # resolve main vectors
  if (is.character(target_total) && length(target_total) == 1) {
    if (!target_total %in% tcols) cli_abort("{.arg target_total} ({target_total}) is not in {.arg target}")
    target_total <- target[[target_total]]
  } else if (length(target_total) != tn) cli_abort("{.arg target_total} is not the same length as {.code nrow(target)}")
  if (is.character(target_indicator) && length(target_indicator) == 1) {
    if (!target_indicator %in% tcols) cli_abort("{.arg target_indicator} ({target_indicator}) is not in {.arg target}")
    target_indicator <- target[[target_indicator]]
  } else if (length(target_indicator) != tn) cli_abort("{.arg target_indicator} is not the same length as {.code nrow(target)}")
  if (is.character(households_indicator) && length(households_indicator) == 1) {
    if (!households_indicator %in% hcols) {
      cli_abort("{.arg households_indicator} ({households_indicator}) is not in {.arg households}")
    }
    households_indicator <- households[[households_indicator]]
  } else if (length(households_indicator) != hn) {
    cli_abort("{.arg households_indicator} is not the same length as {.code nrow(households)}")
  }
  if (!is.logical(target_indicator)) {
    target_indicator <- if (is.numeric(target_indicator)) {
      target_indicator == 1
    } else {
      !grepl("^multi", target_indicator, TRUE)
    }
  }
  if (!is.logical(households_indicator)) {
    households_indicator <- if (is.numeric(households_indicator)) {
      households_indicator == 1
    } else {
      if (all(c("02", "03") %in% households_indicator)) {
        grepl("^0[23]$", households_indicator)
      } else {
        !grepl("^multi", households_indicator, TRUE)
      }
    }
  }

  # resolve size
  size <- NULL
  if (length(households_size) != 1 && length(households_size) != hn) {
    cli_abort("{.arg households_size} is not the same length as {.code nrow(households)}")
  }
  if (length(households_id) != 1 && length(households_id) != hn) {
    cli_abort("{.arg households_id} is not the same length as {.code nrow(households)}")
  }
  if (is.character(households_size) && length(households_size) == 1) {
    if (households_size %in% hcols) {
      size <- households[[households_size]]
    } else {
      if (length(households_id) == 1 && households_id %in% hcols) households_id <- households[[households_id]]
      if (!is.null(person) && length(person_household_id) == 1 && person_household_id %in% colnames(person)) {
        person_household_id <- person[[person_household_id]]
      } else if (length(person_household_id) == 1) {
        cli_abort("{.arg person_household_id} ({person_household_id}) is not in {.arg person}")
      }
      if (!any(households_id %in% person_household_id)) {
        cli_abort("no {.arg households_id}s are present in {.arg person_household_id}")
      }
      size <- table(person_household_id)[households_id]
    }
  }
  if (is.null(size)) cli_abort("failed to resolve {.arg households_size}")

  residents <- target_total * tapply(
    size, households_indicator, mean,
    na.rm = TRUE
  )[as.character(target_indicator)]
  redistribute(source, target, weight = residents, ...)
}
