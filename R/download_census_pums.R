#' Download U.S. Census Microdata
#'
#' Download and load U.S. census American Community Survey (ACS) Public Use Microdata Samples (PUMS):
#' (census.gov/programs-surveys/acs/microdata.html)[https://www.census.gov/programs-surveys/acs/microdata.html]
#'
#'
#' @param dir Directory in which to save the file(s).
#' @param state Postal or FIPS code of the state.
#' @param year 4 digit year, between 2005 and the most recent year.
#' @param level A character indicating whether to get the person- or household-level sample. Defaults to both.
#' @param one_year Logical; if \code{FALSE}, will get the 5-year estimates rather than the 1-year file. If
#' not specified, will fall back to the 5-year file if the 1-year file is not available.
#' @param calculate_error Logical; if \code{TRUE}, will calculate standard errors for each variable using
#' the Successive Difference Replicate Formula from the
#' [PUMS handbook](https://www.census.gov/programs-surveys/acs/library/handbooks/pums.html).
#' @param crosswalk Logical; if \code{FALSE}, will not retrieve the PUMA relationship files for associating
#' Census tracts with PUM areas. This will be treated as \code{TRUE} if
#' \code{geoids} is specified.
#' @param geoids A vector of county, tract, or block group GEOIDs within the specified state to
#' select PUM areas by; defaults to all areas. If specified, \code{crosswalk} will be treated as \code{TRUE}.
#' @param verbose Logical; if \code{FALSE}, will not print status messages.
#' @examples
#' \dontrun{
#' download_census_pums(".", "va")
#' }
#' @return A list with entries for \code{year} and \code{state} (as specified),
#' \code{dictionary} (containing the data dictionary for that year), \code{household} and/or \code{person}
#' (with survey data), and optionally \code{household_error} and/or \code{person_error}
#' (if \code{calculate_error} is \code{TRUE}), and \code{crosswalk} (if \code{crosswalk} is \code{TRUE} or
#' \code{geoids} is specified).
#' @export

download_census_pums <- function(dir, state, year = 2021, level = "both", one_year = TRUE, calculate_error = FALSE,
                                 crosswalk = TRUE, geoids = NULL, verbose = TRUE) {
  us_fips <- list(
    post = c(
      "us", "al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc", "fl", "ga", "hi", "id", "il", "in", "ia", "ks", "ky",
      "la", "me", "md", "ma", "mi", "mn", "ms", "mo", "mt", "ne", "nv", "nh", "nj", "nm", "ny", "nc", "nd", "oh", "ok",
      "or", "pa", "ri", "sc", "sd", "tn", "tx", "ut", "vt", "va", "wa", "wv", "wi", "wy", "pr"
    ),
    fips = c(
      "us", 1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
      33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56, 72
    )
  )
  state <- gsub("^0|[^a-z0-9]", "", tolower(state[[1]]))
  fid <- which(us_fips$fips == state)
  if (!length(fid) && nchar(state) == 2) fid <- which(us_fips$post == state)
  if (!length(fid) || is.na(fid)) cli_abort("failed to recognize {.arg state}")
  state <- us_fips$post[fid]
  level <- tolower(substring(level, 1, 1))
  if (!level %in% c("h", "p", "b")) {
    cli_abort(
      '{.arg level} should be one of {.field "houshold"}, {.field "person"}, or {.field "both"}.'
    )
  }
  summary <- if (one_year) "1-Year" else "5-Year"
  urls <- paste0(
    "https://www2.census.gov/programs-surveys/acs/data/pums/", year,
    "/", summary, "/csv_", c("h", "p"), state, ".zip"
  )
  if (!is.null(dir)) {
    dir <- tdir <- paste0(normalizePath(dir, "/", FALSE), "/")
    if (!dir.exists(dir)) dir.create(dir, FALSE, TRUE)
  } else {
    tdir <- paste0(tempdir(), "/")
  }
  folder <- paste0(tdir, "pums_", state, "/", year, "/")
  dir.create(folder, FALSE, TRUE)
  res <- list(year = year, state = state, dictionary = NULL, household = NULL, person = NULL)
  dict_json <- paste0(folder, paste0("PUMS_Data_Dictionary_", year, ".json"))
  if (file.exists(dict_json)) {
    res$dictionary <- read_json(dict_json, simplifyVector = TRUE)
  } else {
    dict_file <- sub(".json", ".csv", dict_json, fixed = TRUE)
    if (!file.exists(dict_file)) {
      url <- paste0("https://www2.census.gov/programs-surveys/acs/tech_docs/pums/data_dict/", basename(dict_file))
      status <- curl_fetch_disk(url, dict_file)
      if (status$status_code != 200) cli_abort("failed to retrieve data dictionary from {.url {url}}")
    }
    dict <- readLines(dict_file, warn = FALSE)
    res$dictionary <- lapply(
      split(dict, sub(",.*$", "", sub("^[^,]+,", "", dict, perl = TRUE), perl = TRUE)),
      function(e) {
        h <- suppressWarnings(scan(text = e[[1]], what = "", sep = ",", quiet = TRUE))
        map <- as.data.frame(do.call(rbind, lapply(
          e[-1], function(l) suppressWarnings(scan(text = l, what = "", sep = ",", quiet = TRUE))[5:7]
        )))
        colnames(map) <- c("from", "to", "value")
        list(
          description = h[5],
          type = if (h[3] == "C") "character" else "integer",
          length = as.numeric(h[4]),
          map = map
        )
      }
    )
    write_json(res$dictionary, dict_json, auto_unbox = TRUE)
  }
  files <- paste0(folder, "/", c("h", "p"), summary, ".csv.xz")
  if (level != "b") {
    lid <- (level == "p") + 1
    urls <- urls[lid]
    files <- files[lid]
  }
  can_change <- missing(one_year)
  for (i in seq_along(urls)) {
    url <- urls[i]
    clevel <- if (grepl("csv_h", url, fixed = TRUE)) "household" else "person"
    if (!file.exists(files[i])) {
      file_zip <- paste0(folder, basename(url))
      if (!file.exists(file_zip)) {
        if (verbose) cli_alert_info("retrieving {clevel} sample: {.url {url}}")
        status <- curl_fetch_disk(url, file_zip)
        if (status$status_code != 200 && can_change) {
          new_summary <- if (one_year) "5-Year" else "1-Year"
          urls <- sub(summary, new_summary, urls, fixed = TRUE)
          files <- sub(summary, new_summary, files, fixed = TRUE)
          summary <- new_summary
          url <- urls[i]
          status <- curl_fetch_disk(url, file_zip)
        }
        can_change <- FALSE
        if (status$status_code != 200) {
          unlink(file_zip)
          cli_abort(paste0(
            "failed to download {.url {url}}: ", strsplit(rawToChar(status$headers), "[\r\n]")[[1]][1]
          ))
        }
      }
      rfile <- paste0(folder, grep(".csv", unzip(file_zip, list = TRUE)$Name, fixed = TRUE, value = TRUE))
      unzip(file_zip, exdir = folder)
      cols <- strsplit(readLines(rfile, 1), ",")[[1]]
      types <- paste(vapply(cols, function(n) {
        type <- res$dictionary[[n]]$type
        if (is.null(type)) if (grepl("WGTP", n, fixed = TRUE)) "i" else "c" else substring(type, 1, 1)
      }, ""), collapse = "")
      res[[clevel]] <- vroom(rfile, col_types = types)
      vroom_write(res[[clevel]], files[i], ",")
      on.exit(unlink(rfile), add = TRUE)
    } else {
      if (verbose) cli_alert_info("loading {clevel} sample: {.field {basename(files[i])}}")
      cols <- strsplit(readLines(files[i], 1), ",")[[1]]
      types <- paste(vapply(cols, function(n) {
        type <- res$dictionary[[n]]$type
        if (is.null(type)) if (grepl("WGTP", n, fixed = TRUE)) "i" else "c" else substring(type, 1, 1)
      }, ""), collapse = "")
      res[[clevel]] <- vroom(files[i], col_types = types)
    }
    weight <- if (clevel == "person") "PWGTP" else "WGTP"
    if (calculate_error && paste0(weight, "1") %in% cols) {
      if (verbose) cli_alert_info("calculating {clevel} errors")
      vars <- cols[!grepl("^(?:RT|SERIAL|DIV|PUMA|REG|ST|ADJ)|WGTP", cols) & cols %in% names(res$dictionary)]
      names(vars) <- vars
      wv <- res[[clevel]][, weight, drop = TRUE]
      w <- as.matrix(res[[clevel]][, grep("WGTP\\d+", cols)])
      res[[paste0(clevel, "_error")]] <- as.data.frame(do.call(cbind, lapply(vars, function(name) {
        v <- res[[clevel]][, name, drop = TRUE]
        map <- res$dictionary[[name]]
        if (map$type == "integer") {
          sqrt(.05 * rowSums((w - v * wv)^2))
        } else {
          levels <- map$map$from
          names(levels) <- paste0(name, "_", levels)
          vapply(levels, function(l) {
            sqrt(.05 * rowSums((w - (v == l) * wv)^2))
          }, numeric(length(v)))
        }
      })))
    }
  }
  if (crosswalk || !is.null(geoids)) {
    decade <- paste0(substring(year, 1, 3), "0")
    url <- if (decade == "2020") {
      "https://www2.census.gov/geo/docs/maps-data/data/rel2020/2020_Census_Tract_to_2020_PUMA.txt"
    } else {
      "https://www2.census.gov/geo/docs/maps-data/data/rel/2010_Census_Tract_to_2010_PUMA.txt"
    }
    crosswalk_file <- paste0(folder, basename(url))
    if (!file.exists(crosswalk_file)) {
      if (verbose) cli_alert_info("retrieving crosswalk {.url {url}}")
      status <- curl_fetch_disk(url, crosswalk_file)
      if (status$status_code != 200) cli_warn("failed to retrieve crosswalk from {.url {url}}")
    }
    if (file.exists(crosswalk_file)) {
      if (verbose) cli_alert_info("loading crosswalk {.field {basename(crosswalk_file)}}")
      res$crosswalk <- vroom(crosswalk_file, col_types = "cccc")
      if (!is.null(geoids)) {
        geoids <- substring(geoids, 1, 11)
        su <- do.call(paste0, res$crosswalk[, 1:2]) %in% geoids | do.call(paste0, res$crosswalk[, 1:3]) %in% geoids
        if (any(su)) {
          pids <- unique(res$crosswalk$PUMA5CE[su])
          if (verbose) {
            cli_alert_info(
              "subsetting to {length(pids)} of {length(unique(res$crosswalk$PUMA5CE))} PUM areas"
            )
          }
          for (set in c("household", "person")) {
            su <- formatC(as.numeric(res[[set]]$PUMA), width = 5, flag = 0) %in% pids
            res[[set]] <- res[[set]][su, ]
            set_error <- paste0(set, "_error")
            if (set_error %in% names(res)) res[[set_error]] <- res[[set_error]][su, ]
          }
        } else {
          cli_warn("none of the specified {.arg geoids} were found in the crosswalk")
        }
      }
    }
  }
  res
}
