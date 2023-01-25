#' Redistribute Data
#'
#' Distribute data from a source frame to a target frame.
#'
#' @param source A matrix-like object you want to distribute from; usually this will be
#' the real or more complete dataset, and is often at a lower resolution / higher level.
#' @param target A matrix-like object you want to distribute to: usually this will be
#' the dataset you want but isn't available, and is often at a higher resolution / lower level.
#' @param map A list with entries named with \code{source} IDs (or aligning with those IDs),
#' containing vectors of associated \code{target} IDs (or indices of those IDs). Entries
#' can also be numeric vectors with IDs as names, which will be used to weigh the relationship.
#' If IDs are related by substrings (the first characters of \code{target} IDs are \code{source} IDs),
#' then a map can be automatically generated from them. If \code{source} and \code{target}
#' contain \code{sf} geometries, a map will be made with \code{\link[sf]{st_intersects}}
#' (\code{st_intersects(source, target)}). If an intersects map is made, and \code{source}
#' is being aggregated to \code{target}, and map entries contain multiple target IDs,
#' those entries will be weighted by their proportion of overlap with the source area.
#' @param source_id,target_id Name of a column in \code{source} / \code{target},
#' or a vector containing IDs. For \code{source}, this will default to the first column. For
#' \code{target}, columns will be searched through for one that appears to relate to the
#' source IDs, falling back to the first column.
#' @param weight Name of a column, or a vector containing weights (or single value to apply to all cases),
#' which apply to \code{target} when disaggregating, and \code{source} when aggregating.
#' Defaults to unit weights (all weights are 1).
#' @param source_variable,source_value If \code{source} is tall (with variables spread across
#' rows rather than columns), specifies names of columns in \code{source} containing variable names
#' and values for conversion.
#' @param aggregate Logical; if specified, will determine whether to aggregate or disaggregate
#' from \code{source} to \code{target}. Otherwise, this will be \code{TRUE} if there are more
#' \code{source} observations than \code{target} observations.
#' @param weight_agg_method Means of aggregating \code{weight}, in the case that target IDs contain duplicates.
#' Options are \code{"sum"}, \code{"average"}, or \code{"auto"} (default; which will sum if \code{weight}
#' is integer-like, and average otherwise).
#' @param outFile Path to a CSV file in which to save results.
#' @param overwrite Logical; if \code{TRUE}, will overwrite an existing \code{outFile}.
#' @param make_intersect_map Logical; if \code{TRUE}, will opt to calculate an intersect-based map
#' rather than an ID-based map, if both seem possible. If specified as \code{FALSE}, will
#' never calculate an intersect-based map.
#' @param overlaps If specified and not \code{TRUE} or \code{"keep"}, will assign \code{target}
#' entities that are mapped to multiple \code{source} entities to a single source entity. The value
#' determines how entities with the same weight should be assigned, between \code{"first"} (default),
#' \code{"last"}, and \code{"random"}.
#' @param use_all Logical; if \code{TRUE} (default), will redistribute map weights so they sum to 1.
#' Otherwise, entities may be partially weighted.
#' @param return_geometry Logical; if \code{FALSE}, will not set the returned \code{data.frame}'s
#' geometry to that of \code{target}, if it exists.
#' @param return_map Logical; if \code{TRUE}, will only return the map, without performing the
#' redistribution. Useful if you want to inspect an automatically created map, or use it in a later call.
#' @param verbose Logical; if \code{TRUE}, will show status messages.
#' @examples
#' # minimal example
#' source <- c(a = 1, b = 2)
#' target <- 1:5
#' (redistribute(source, target, verbose = TRUE))
#'
#' # multi-entity example
#' source <- data.frame(id = c("a", "b"), cat = c("aaa", "bbb"), num = c(1, 2))
#' target <- data.frame(
#'   id = sample(paste0(c("a", "b"), rep(1:5, 2))),
#'   population = sample.int(1e5, 10)
#' )
#' (redistribute(source, target, verbose = TRUE))
#' @returns A \code{data.frame} with a row for each \code{target_ids} (identified by the first column,
#' \code{id}), and a column for each variable from \code{source}.
#' @importFrom cli cli_abort cli_bullets cli_alert_info cli_alert_warning cli_warn
#' cli_progress_step cli_progress_update cli_progress_done
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom utils unzip
#' @importFrom sf st_intersects st_intersection st_geometry st_geometry<- st_crs st_crs<- st_geometry_type
#' st_coordinates st_centroid st_boundary st_cast st_polygon st_union
#' @importFrom s2 s2_area
#' @importFrom lingmatch lma_simets
#' @importFrom jsonlite read_json write_json
#' @importFrom curl curl_fetch_disk
#' @importFrom vroom vroom vroom_write
#' @importFrom stats rbeta rbinom rnorm rpois
#' @importFrom Matrix Matrix drop0
#' @useDynLib redistribute, .registration = TRUE
#' @export

redistribute <- function(source, target = NULL, map = list(), source_id = "GEOID",
                         target_id = source_id, weight = NULL, source_variable = NULL, source_value = NULL,
                         aggregate = NULL, weight_agg_method = "auto", outFile = NULL, overwrite = FALSE,
                         make_intersect_map = FALSE, overlaps = "keep", use_all = TRUE, return_geometry = TRUE,
                         return_map = FALSE, verbose = FALSE) {
  if (!overwrite && !is.null(outFile) && file.exists(outFile)) {
    cli_abort("{.arg outFile} already exists; use {.code overwrite = TRUE} to overwrite it")
  }
  can_intersects <- missing(make_intersect_map) || make_intersect_map
  source_sf <- can_intersects && inherits(source, "sf")
  target_sf <- can_intersects && inherits(target, "sf")
  intersect_map <- FALSE
  if (length(dim(source)) != 2) source <- t(source)
  if (is.null(colnames(source))) colnames(source) <- paste0("V", seq_len(ncol(source)))
  if (length(source_id) > 1) {
    if (verbose) cli_alert_info("source IDs: {.arg source_id} vector")
    sid <- source_id
  } else if (source_id %in% colnames(source)) {
    if (verbose) cli_alert_info("source IDs: {.field {source_id}} column of {.arg source}")
    sid <- source[, source_id, drop = TRUE]
    source <- source[, colnames(source) != source_id]
  } else if (nrow(source) == 1) {
    if (verbose) cli_alert_info("source IDs: {.field 1}")
    sid <- if (missing(source_id)) "1" else source_id
  } else {
    if (source_sf && target_sf) {
      if (verbose) cli_alert_info("source IDs: sequence, assuming map from geometries")
      intersect_map <- TRUE
      sid <- seq_len(nrow(source))
    } else {
      if (verbose) cli_alert_info("source IDs: {.field {colnames(source)[1]}} column of {.arg source}")
      sid <- source[, 1, drop = TRUE]
      source <- source[, -1, drop = FALSE]
    }
  }
  sid <- as.character(sid)
  if (length(sid) != nrow(source)) cli_abort("{.arg source_id} is not the same length as {.field nrow(source)}")
  if (anyDuplicated(sid)) {
    if (is.null(source_variable)) {
      if ("variable" %in% colnames(source)) source_variable <- "variable"
      if ("measure" %in% colnames(source)) source_variable <- "measure"
    }
    if (is.null(source_value)) {
      if ("value" %in% colnames(source)) source_value <- "value"
    }
  }
  if (!is.null(source_variable) || !is.null(source_value)) {
    source_numeric <- vapply(seq_len(ncol(source)), function(col) is.numeric(source[, col, drop = TRUE]), TRUE)
    if (is.null(source_value) && any(source_numeric)) source_value <- colnames(source)[which(source_numeric)[1]]
    if (is.null(source_variable) && any(!source_numeric)) {
      source_variable <- colnames(source)[
        which.min(vapply(which(!source_numeric), function(col) length(unique(source[, col])), 0))
      ]
    }
    if (verbose) {
      cli_alert_info(paste(
        "converting {.arg source} format, breaking",
        if (length(source_value) == 1) "{.field {source_value}}" else "{.arg source_value}",
        "values across",
        if (length(source_variable) == 1) "{.field {source_variable}}" else "{.arg source_variable}",
        "columns"
      ))
    }
    usid <- unique(sid)
    nr <- nrow(source)
    if (length(source_value) != nr) source_value <- source[, source_value, drop = TRUE]
    if (length(source_variable) != nr) source_variable <- source[, source_variable, drop = TRUE]
    ss <- split(data.frame(sid, source_value), source_variable)
    source <- do.call(cbind, lapply(ss, function(vd) {
      vd <- vd[!duplicated(vd[[1]]), ]
      rownames(vd) <- vd[, 1]
      vd[usid, 2, drop = TRUE]
    }))
    colnames(source) <- names(ss)
    source <- as.data.frame(source)
    sid <- usid
  }
  if (length(dim(target)) == 2 && is.null(colnames(target))) colnames(target) <- paste0("V", seq_len(ncol(target)))
  if (length(target_id) > 1) {
    if (verbose) cli_alert_info("target IDs: {.arg target_id} vector")
    tid <- target_id
  } else if (is.null(target)) {
    if (verbose) cli_alert_info("target IDs: {.field 1}")
    tid <- 1
  } else if (!is.null(target_id) && length(dim(target)) == 2 && target_id %in% colnames(target)) {
    if (verbose) cli_alert_info("target IDs: {.field {target_id}} column of {.arg target}")
    tid <- target[, target_id, drop = TRUE]
    target <- target[, colnames(target) != target_id]
  } else {
    if (length(dim(target)) != 2) {
      if (is.null(weight) && is.numeric(target) && !is.null(names(target))) {
        if (verbose) cli_alert_info("target IDs: names of {.arg target} vector; {.arg weight} set to {.arg target}")
        weight <- unname(target)
        tid <- names(target)
      } else {
        if (verbose) cli_alert_info("target IDs: {.arg target} vector")
        tid <- target
      }
    } else {
      if (can_intersects && target_sf && source_sf) {
        if (verbose) cli_alert_info("target IDs: sequence, assuming map from geometries")
        intersect_map <- TRUE
        tid <- seq_len(nrow(target))
      } else {
        idi <- 1
        if (length(map)) {
          tids <- unique(unlist(map, use.names = FALSE))
          for (i in seq_len(ncol(target))) {
            if (any(target[, i] %in% tids)) {
              idi <- i
              break
            }
          }
        } else if (all(nchar(sid) == nchar(sid[1]))) {
          sl <- nchar(sid[1])
          for (i in seq_len(ncol(target))) {
            if (any(substr(target[, i], 1, sl) %in% sid)) {
              idi <- i
              break
            }
          }
        }
        if (verbose) cli_alert_info("target IDs: {.field {colnames(target)[idi]}} column of {.arg target}")
        tid <- target[, idi, drop = TRUE]
        target <- target[, -idi, drop = FALSE]
      }
    }
  }
  tid <- as.character(tid)
  if (!is.logical(aggregate)) aggregate <- length(sid) > length(tid)
  if (make_intersect_map && source_sf && target_sf) intersect_map <- TRUE
  if (length(map)) {
    if (verbose) cli_alert_info("map: provided list")
    intersect_map <- FALSE
    if (is.null(names(map))) {
      if (length(map) != length(sid)) cli_abort("{.arg map} has no names, and is not the same length as source IDs")
      names(map) <- sid
    } else if (!any(sid %in% names(map))) cli_abort("no source IDs were present in the ID map")
  } else {
    if (nrow(source) == 1) {
      if (verbose) cli_alert_info("map: all target IDs for single source")
      intersect_map <- FALSE
      map[[sid]] <- tid
    } else if (is.null(target) || length(tid) == 1) {
      if (verbose) cli_alert_info("map: all source IDs to a single target")
      intersect_map <- FALSE
      map <- as.list(structure(rep(tid, length(sid)), names = sid))
    } else if (aggregate && !intersect_map && all(nchar(tid) == nchar(tid[1])) &&
      nchar(sid[1]) > nchar(tid[1]) && any(substring(sid, 1, nchar(tid[1])) %in% tid)) {
      if (verbose) cli_alert_info("map: first {.field {nchar(tid[1])}} character{?s} of source IDs")
      intersect_map <- FALSE
      map <- as.list(structure(substr(sid, 1, nchar(tid[1])), names = sid))
    } else if (!intersect_map && all(nchar(sid) == nchar(sid[1])) && nchar(tid[1]) > nchar(sid[1]) &&
      any(substring(tid, 1, nchar(sid[1])) %in% sid)) {
      if (verbose) cli_alert_info("map: first {.field {nchar(sid[1])}} character{?s} of target IDs")
      map <- split(tid, substr(tid, 1, nchar(sid[1])))
    } else if (can_intersects && source_sf && target_sf) {
      if (verbose) cli_alert_info("map: intersections between geometries")
      intersect_map <- TRUE
      op <- options(sf_use_s2 = FALSE)
      on.exit(options(sf_use_s2 = op[[1]]))
      if (st_crs(source) != st_crs(target)) st_crs(target) <- st_crs(source)
      map <- tryCatch(suppressMessages(st_intersects(st_geometry(source), st_geometry(target))), error = function(e) NULL)
      if (is.null(map)) {
        cli_abort(c(
          x = "{.fn st_intersects} failed between {.arg source} and {.arg target}",
          i = "if you can do this or similar separetely, you can pass the result as {.arg map}"
        ))
      }
      names(map) <- sid
    } else {
      cli_abort("no map was provided, and could not make one from IDs")
    }
  }
  if (aggregate && length(map) == length(tid) && length(map) != length(sid)) {
    ids <- if (is.null(names(map))) seq_along(map) else names(map)
    names(map) <- NULL
    child_counts <- vapply(map, length, 0)
    map <- as.list(unlist(map))
    names(map) <- rep(ids, child_counts)
  }
  if (intersect_map) {
    if (st_crs(source) != st_crs(target)) st_crs(target) <- st_crs(source)
    source_geom <- st_geometry(source)
    names(source_geom) <- sid
  }
  if ((is.logical(overlaps) && overlaps) || grepl("^[Kk]", overlaps)) {
    dodedupe <- FALSE
  } else {
    tiebreak <- c(f = "first", l = "last", r = "random")[tolower(substring(overlaps, 1, 1))]
    dodedupe <- !is.na(tiebreak)
    deduper <- function(e, tm = tiebreak) {
      su <- which(e == max(e))
      if (length(su) > 1) {
        su <- switch(tm,
          first = su[1],
          last = su[length(su)],
          random = sample(su, 1)
        )
      }
      e[su]
    }
  }
  any_partial <- FALSE
  map <- map[sid]
  mls <- vapply(map, length, 0)
  polys <- intersect_map && any(grepl("POLY", st_geometry_type(
    if (aggregate) source_geom else st_geometry(target)
  ), fixed = TRUE))
  if (polys && any(mls > 1)) {
    map <- lapply(sid, function(id) {
      e <- map[[id]]
      if (length(e)) {
        res <- if (is.null(names(e))) {
          w <- rep(1, length(e))
          if (intersect_map && length(e) > 1) {
            reg <- st_geometry(target)[e]
            totals <- s2_area(reg)
            su <- totals > 0
            if (any(su)) {
              part <- suppressMessages(st_intersection(reg[su], source_geom[id]))
              w[su] <- s2_area(part) / (if (aggregate) s2_area(source_geom[[id]]) else totals[su])
            }
          }
          names(w) <- if (is.integer(e)) tid[e] else e
          w
        } else {
          e
        }
        if (dodedupe) {
          deduper(res)
        } else {
          if (!aggregate && any(res < 1)) any_partial <<- TRUE
          res
        }
      } else {
        numeric()
      }
    })
  } else {
    map <- lapply(map, function(e) {
      if (is.null(names(e))) {
        w <- rep(1, length(e))
        names(w) <- if (is.integer(e)) tid[e] else e
        e <- w
      }
      if (!aggregate && any(e < 1)) any_partial <<- TRUE
      e
    })
  }
  names(map) <- sid
  if (return_map) {
    if (verbose) cli_alert_info("returning map")
    return(map)
  }
  mtid <- unlist(lapply(map, names), use.names = FALSE)
  nout <- length(if (aggregate) sid else tid)
  w <- if (length(weight) > 1) {
    if (verbose) cli_alert_info("weights: {.arg weight} vector")
    if (length(weight) != nout) {
      cli_abort("{.arg weight} is not the same length as {.arg {if (aggregate) 'source' else 'target'}} IDs")
    }
    weight
  } else if (aggregate && !is.null(weight) && weight %in% colnames(source)) {
    if (verbose) cli_alert_info("weights: {.field {weight}} column of {.arg source}")
    w <- source[, weight, drop = TRUE]
    source <- source[, colnames(weight) != weight, drop = FALSE]
    w
  } else if (!is.null(weight) && weight %in% colnames(target)) {
    if (verbose) cli_alert_info("weights: {.field {weight}} column of {.arg target}")
    target[, weight, drop = TRUE]
  } else {
    weight <- if (!is.numeric(weight)) 1 else weight
    if (verbose) cli_alert_info("weights: {.field {weight}}")
    rep(weight, nout)
  }
  w <- as.numeric(w)
  realign <- FALSE
  su <- tid %in% mtid
  if (!all(su)) {
    if (verbose) cli_alert_info("some {.arg target_id}s were dropped because they were not present in {.arg map}")
    realign <- TRUE
    otid <- tid
    tid <- tid[su]
    if (!aggregate) w <- w[su]
  }
  if (!aggregate && anyDuplicated(tid)) {
    if (!realign) {
      realign <- TRUE
      otid <- tid
    }
    agger <- if (weight_agg_method == "auto") {
      if (all(w == as.integer(w))) "summing" else "averaging"
    } else if (grepl("^[Ss]", weight_agg_method)) "summing" else "averaging"
    if (verbose) cli_alert_info("{.arg target_id} contains duplicates, so {agger} {.arg weight}s")
    w <- tapply(w, tid, if (agger == "summing") sum else mean)
    tid <- names(w)
    w <- unname(w)
  }
  if (inherits(source, "sf")) st_geometry(source) <- NULL
  source_numeric <- vapply(seq_len(ncol(source)), function(col) is.numeric(source[, col, drop = TRUE]), TRUE)
  if (!all(source_numeric)) {
    non_numeric <- which(!source_numeric)
    level_map <- list()
    for (col in non_numeric) {
      x <- as.factor(source[, col, drop = TRUE])
      source[, col] <- as.numeric(x)
      level_map[[as.character(col)]] <- list(
        index = col,
        levels = levels(x)
      )
    }
  }
  if (verbose) {
    cli_alert_info("{if (aggregate) 'aggregating' else 'disaggregating'} {.field {length(source_numeric)}} variable{?s}:")
    var_groups <- split(colnames(source), !source_numeric)
    names(var_groups) <- c(
      "TRUE" = "(char; {sum(!source_numeric)})", "FALSE" = "(numb; {sum(source_numeric)})"
    )[names(var_groups)]
    cli_bullets(structure(
      vapply(names(var_groups), function(type) {
        paste(type, if (length(var_groups[[type]]) < 10) paste(var_groups[[type]], collapse = ", "))
      }, ""),
      names = rep("*", length(var_groups))
    ))
  }
  method <- as.integer(source_numeric)
  if (any_partial) aggregate <- TRUE
  if (aggregate) {
    method <- method + 10
    if (any_partial) method[method == 11] <- 12
  }
  balance <- any_partial && (!is.logical(use_all) || use_all)
  res <- process_distribute(as.matrix(source), method, tid, w, map, aggregate, balance)
  res <- as.data.frame(res)
  colnames(res) <- colnames(source)
  if (!all(source_numeric)) {
    if (verbose) cli_alert_info("re-converting categorical levels")
    for (l in level_map) {
      x <- res[, l$index]
      x[x == 0] <- NA
      res[, l$index] <- l$levels[x]
    }
  }
  res <- cbind(id = tid, res)
  if (realign) {
    if (verbose) cli_alert_info("realigning with original target IDs")
    rownames(res) <- tid
    res <- res[otid, ]
    res$id <- otid
    rownames(res) <- NULL
  }
  if (!is.null(outFile)) {
    if (verbose) cli_alert_info("writing results to {.file {outFile}}")
    dir.create(dirname(outFile), FALSE, TRUE)
    if (!grepl("\\.\\w", outFile)) outFile <- paste0(outFile, ".csv")
    vroom_write(res, outFile, ",")
  }
  if (return_geometry && target_sf) st_geometry(res) <- st_geometry(target)
  invisible(res)
}
