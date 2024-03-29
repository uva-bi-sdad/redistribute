#' Redistribute Data
#'
#' Distribute data from a source frame to a target frame.
#'
#' @param source A matrix-like object you want to distribute from; usually this will be
#' the real or more complete dataset, and is often at a lower resolution / higher level.
#' @param target A matrix-like object you want to distribute to: usually this will be
#' the dataset you want but isn't available, and is often at a higher resolution / lower level
#' (for disaggregation). Can also be a single number, representing the number of initial
#' characters of \code{source} IDs to derive target IDs from (useful for aggregating up nested groups).
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
#' @param rescale Logical; if \code{FALSE}, will not adjust target values after redistribution such that they
#' match source totals.
#' @param drop_extra_sources Logical; if \code{TRUE}, will remove any source rows that are not mapped
#' to any target rows. Useful when inputting a source with regions outside of the target area,
#' especially when \code{rescale} is \code{TRUE}.
#' @param default_value Value to set to any unmapped target ID.
#' @param outFile Path to a CSV file in which to save results.
#' @param overwrite Logical; if \code{TRUE}, will overwrite an existing \code{outFile}.
#' @param make_intersect_map Logical; if \code{TRUE}, will opt to calculate an intersect-based map
#' rather than an ID-based map, if both seem possible. If specified as \code{FALSE}, will
#' never calculate an intersect-based map.
#' @param fill_targets Logical; if \code{TRUE}, will make new \code{target} rows for any
#' un-mapped \code{source} row.
#' @param overlaps If specified and not \code{TRUE} or \code{"keep"} (default), will assign \code{target}
#' entities that are mapped to multiple \code{source} entities to a single source entity. The value
#' determines how entities with the same weight should be assigned, between \code{"first"},
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
#' source <- data.frame(a = 1, b = 2)
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
#' @importFrom sf st_intersects st_intersection st_geometry st_geometry<- st_crs st_geometry_type
#' st_coordinates st_centroid st_boundary st_cast st_polygon st_union st_transform st_buffer st_as_sf
#' st_is_valid st_make_valid st_sfc st_point
#' @importFrom s2 s2_area s2_is_valid
#' @importFrom lingmatch lma_simets
#' @importFrom jsonlite read_json write_json
#' @importFrom curl curl_fetch_disk
#' @importFrom vroom vroom vroom_write
#' @importFrom stats rbeta rbinom rnorm rpois
#' @importFrom Matrix Matrix drop0
#' @useDynLib redistribute, .registration = TRUE
#' @export

redistribute <- function(source, target = NULL, map = list(), source_id = "GEOID", target_id = source_id,
                         weight = NULL, source_variable = NULL, source_value = NULL, aggregate = NULL,
                         weight_agg_method = "auto", rescale = TRUE, drop_extra_sources = FALSE,
                         default_value = NA, outFile = NULL, overwrite = FALSE, make_intersect_map = FALSE,
                         fill_targets = FALSE, overlaps = "keep", use_all = TRUE, return_geometry = TRUE,
                         return_map = FALSE, verbose = FALSE) {
  if (!overwrite && !is.null(outFile) && file.exists(outFile)) {
    cli_abort("{.arg outFile} already exists; use {.code overwrite = TRUE} to overwrite it")
  }
  can_intersects <- missing(make_intersect_map) || make_intersect_map
  source_sf <- can_intersects && inherits(source, c("sfc", "sf"))
  target_sf <- can_intersects && inherits(target, c("sfc", "sf"))
  can_intersects <- source_sf && target_sf
  if (!can_intersects && make_intersect_map) {
    if (!source_sf && is.data.frame(source) && any(vapply(source, inherits, TRUE, "sfc"))) {
      cli_warn("{.arg source} is not an sf object, but contains a simple features column, so converting it")
      source <- st_as_sf(source)
      source_sf <- TRUE
    }
    if (!target_sf && is.data.frame(target) && any(vapply(target, inherits, TRUE, "sfc"))) {
      cli_warn("{.arg target} is not an sf object, but contains a simple features column, so converting it")
      target <- st_as_sf(target)
      target_sf <- TRUE
    }
    can_intersects <- source_sf && target_sf
    if (!can_intersects) {
      cli_abort(
        "{.arg make_intersect_map} was set to {.field TRUE}, but {.arg source} and {.arg target} are not both sf objects"
      )
    }
  }
  intersect_map <- FALSE
  if (length(dim(source)) != 2) source <- data.frame(source)
  if (is.null(colnames(source))) colnames(source) <- paste0("V", seq_len(ncol(source)))
  if (length(source_id) > 1) {
    if (verbose) cli_alert_info("source IDs: {.arg source_id} vector")
    sid <- source_id
  } else if (source_id %in% colnames(source)) {
    if (verbose) cli_alert_info("source IDs: {.field {source_id}} column of {.arg source}")
    sid <- source[, source_id, drop = TRUE]
    source <- source[, colnames(source) != source_id, drop = FALSE]
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
  if (length(target) == 1 && is.numeric(target)) {
    if (verbose) cli_alert_info("target IDs: first {.field {target}} characters of {.arg source} IDs")
    tid <- target <- unique(substring(sid, 1, target))
  } else if (length(target_id) > 1) {
    if (verbose) cli_alert_info("target IDs: {.arg target_id} vector")
    tid <- target_id
  } else if (is.null(target)) {
    if (verbose) cli_alert_info("target IDs: {.field 1}")
    tid <- 1
  } else if (!is.null(target_id) && length(dim(target)) == 2 && target_id %in% colnames(target)) {
    if (verbose) cli_alert_info("target IDs: {.field {target_id}} column of {.arg target}")
    tid <- target[, target_id, drop = TRUE]
    target <- target[, colnames(target) != target_id, drop = FALSE]
  } else {
    if (length(dim(target)) != 2) {
      if (is.null(weight) && (target_sf || is.numeric(target)) && !is.null(names(target))) {
        if (target_sf) {
          if (verbose) cli_alert_info("target IDs: names of {.arg target} vector")
        } else {
          if (verbose) cli_alert_info("target IDs: names of {.arg target} vector; {.arg weight} set to {.arg target}")
          weight <- unname(target)
        }
        tid <- names(target)
      } else {
        if (target_sf) {
          if (verbose) cli_alert_info("target IDs: sequence")
          tid <- seq_along(target)
        } else {
          if (verbose) cli_alert_info("target IDs: {.arg target} vector")
          tid <- target
        }
      }
    } else {
      if (can_intersects) {
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
        } else if (length(sid) > 1 && (sl <- nchar(sid[1])) > 1 && all(nchar(sid) == sl)) {
          for (i in seq_len(ncol(target))) {
            if (any(substr(target[, i, drop = TRUE], 1, sl) %in% sid)) {
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
    } else if (!any(sid %in% names(map))) {
      if (!any(tid %in% names(map))) cli_abort("no source or target IDs were in map names")
      mw <- unlist(unname(map))
      if (is.null(names(mw))) mw <- structure(numeric(length(mw)) + 1, names = mw)
      onames <- names(mw)
      if (!any(sid %in% onames)) {
        if (all(!grepl("[^0-9]", onames))) {
          onames <- as.integer(onames)
          if (!anyNA(onames) && min(onames) > 0 && max(onames) <= length(sid)) {
            if (verbose) cli_alert_info("inverting map, with entry values mapped to source IDs")
            names(mw) <- rep(names(map), vapply(map, length, 0))
            map <- split(mw, sid[onames])
          } else {
            cli_abort("map appeared to refer to target indices, but they are out of range")
          }
        } else {
          cli_abort("no source IDs were present in the ID map")
        }
      } else {
        if (verbose) cli_alert_info("inverting map")
        names(mw) <- rep(names(map), vapply(map, length, 0))
        map <- split(mw, onames)
      }
    }
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
      any(substring(tid, 1, nchar(sid[1])) %in% sid) &&
      (!can_intersects || !all(substring(tid, 1, nchar(sid[1])) == substring(tid[1], 1, nchar(sid[1]))))) {
      if (verbose) cli_alert_info("map: first {.field {nchar(sid[1])}} character{?s} of target IDs")
      map <- split(tid, substr(tid, 1, nchar(sid[1])))
    } else if (can_intersects) {
      if (verbose) {
        cli_alert_info("map: intersections between geometries")
        cli_progress_step("initial mapping...", msg_done = "done initial mapping")
      }
      intersect_map <- TRUE
      op <- options(sf_use_s2 = FALSE)
      on.exit(options(sf_use_s2 = op[[1]]))
      if (st_crs(source) != st_crs(target)) {
        if (length(tid) > length(sid)) {
          source <- st_transform(source, st_crs(target))
        } else {
          target <- st_transform(target, st_crs(source))
        }
      }
      su <- !st_is_valid(st_geometry(source))
      if (any(su)) {
        st_geometry(source)[su] <- st_make_valid(st_geometry(source)[su])
        su <- !st_is_valid(st_geometry(source))
        if (any(su)) {
          cli_warn("some {.arg source} geometries were not valid, so removed unrepairable ones")
          source <- source[!su, ]
          sid <- sid[!su]
        } else {
          cli_warn("some {.arg source} geometries were not valid, so repaired them")
        }
      }
      su <- !st_is_valid(st_geometry(target))
      if (any(su)) {
        st_geometry(target)[su] <- st_make_valid(st_geometry(target)[su])
        su <- !st_is_valid(st_geometry(target))
        if (any(su)) {
          cli_warn("some {.arg target} geometries were not valid, so removed unrepairable ones")
          target <- target[!su, ]
          tid <- tid[!su]
        } else {
          cli_warn("some {.arg target} geometries were not valid, so repaired them")
        }
      }
      map <- tryCatch(suppressMessages(st_intersects(st_geometry(source), st_geometry(target))), error = function(e) NULL)
      if (is.null(map)) {
        cli_abort(c(
          x = "{.fn st_intersects} failed between {.arg source} and {.arg target}",
          i = "if you can do this or similar separetely, you can pass the result as {.arg map}"
        ))
      }
      names(map) <- sid
    } else {
      cli_abort("no map was provided, and it could not be made from IDs")
    }
  }
  if (aggregate && length(map) == length(tid) && length(map) != length(sid)) {
    ids <- if (is.null(names(map))) seq_along(map) else names(map)
    names(map) <- NULL
    child_counts <- vapply(map, length, 0)
    map <- as.list(unlist(map))
    names(map) <- rep(ids, child_counts)
  }
  if ((is.logical(overlaps) && overlaps) || grepl("^[Kk]", overlaps)) {
    dodedupe <- FALSE
  } else {
    tiebreak <- c(f = "first", l = "last", r = "random")[tolower(substring(overlaps, 1, 1))]
    dodedupe <- !is.na(tiebreak)
    deduper <- switch(tiebreak,
      first = function(e, agg = TRUE) {
        su <- which(e == max(e))
        if (length(su) > 1) su <- su[1]
        if (agg) e[su] else su
      },
      last = function(e, agg = TRUE) {
        su <- which(e == max(e))
        if (length(su) > 1) su <- su[length(su)]
        if (agg) e[su] else su
      },
      random = function(e, agg = TRUE) {
        su <- which(e == max(e))
        if (length(su) > 1) su <- sample(su, 1)
        if (agg) e[su] else su
      }
    )
  }
  any_partial <- FALSE
  map <- map[sid]
  mls <- vapply(map, length, 0)
  if (!fill_targets && drop_extra_sources) {
    su <- vapply(map, length, 0) != 0
    if (!all(su)) {
      if (verbose) cli_alert_info("removing {sum(!su)} {.arg source}{?s} with no mapped {.arg target}s")
      source <- source[su, , drop = FALSE]
      sid <- sid[su]
      map <- map[su]
    }
  }
  if (intersect_map) {
    source_geom <- st_geometry(source)
    names(source_geom) <- sid
  }
  polys <- intersect_map && any(grepl("POLY", st_geometry_type(
    if (aggregate) source_geom else st_geometry(target)
  ), fixed = TRUE))
  if (polys && any(mls > 1)) {
    if (verbose) {
      oi <- 0
      cli_progress_step(
        "calculating map intersections ({oi}/{length(map)})",
        msg_done = "calculated map intersections ({length(map)})"
      )
    }
    map <- lapply(seq_along(map), function(i) {
      if (verbose) {
        oi <<- i
        cli_progress_update(.envir = parent.frame(2))
      }
      e <- map[[i]]
      if (length(e)) {
        res <- if (is.null(names(e))) {
          w <- rep(1, length(e))
          if (intersect_map && length(e) > 1) {
            reg <- st_geometry(target)[e]
            totals <- s2_area(reg)
            su <- totals > 0
            if (any(su)) {
              part <- suppressMessages(st_intersection(reg[su], source_geom[sid[i]], model = "closed"))
              pv <- numeric(length(part))
              tsu <- s2_is_valid(part)
              pv[tsu] <- s2_area(part[tsu])
              w[su] <- pv / (if (aggregate) s2_area(source_geom[[sid[i]]]) else totals[su])
            }
          }
          names(w) <- if (is.integer(e)) tid[e] else e
          w[w > 0]
        } else {
          e
        }
        if (aggregate && dodedupe) {
          if (length(res) > 1) {
            deduper(res)
          } else {
            res
          }
        } else {
          if (!aggregate && any(res < 1)) any_partial <<- TRUE
          res
        }
      } else {
        numeric()
      }
    })
    if (verbose) cli_progress_done()
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
  targets <- unlist(unname(map))
  if (!any_partial) any_partial <- anyDuplicated(names(targets))
  if (dodedupe && any_partial) {
    if (verbose) cli_progress_step("assigning each target to a single source")
    any_partial <- FALSE
    esids <- factor(rep(sid, vapply(map, length, 0)), unique(sid))
    if (anyDuplicated(names(targets))) {
      for (ct in unique(names(targets))) {
        su <- which(names(targets) == ct)
        if (length(su) > 1) {
          sel <- deduper(targets[su], FALSE)
          targets[su[-sel]] <- 0
          targets[su[sel]] <- 1
        }
      }
      su <- which(targets != 0)
      map <- split(targets[su], esids[su])
    } else {
      targets[] <- 1L
      map <- split(targets, esids)
    }
    if (verbose) cli_progress_done()
  } else {
    names(map) <- sid
  }
  n_filled <- 0
  if (fill_targets) {
    su <- vapply(map, length, 0) == 0
    if (any(su)) {
      missed_sources <- names(which(su))
      ftids <- paste0("filled_", missed_sources)
      n_filled <- length(ftids)
      if (verbose) cli_alert_info("adding {n_filled} target ID{?s} to cover all source rows")
      tid <- c(tid, ftids)
      map[missed_sources] <- split(structure(rep(1, n_filled), names = ftids), missed_sources)
    }
  }
  if (verbose && intersect_map) cli_progress_done()
  if (return_map) {
    if (verbose) cli_alert_info("returning map")
    return(map)
  }
  mtid <- names(unlist(unname(map)))
  nout <- length(if (aggregate) sid else tid)
  w <- if (length(weight) > 1) {
    if (verbose) cli_alert_info("weights: {.arg weight} vector")
    if (length(weight) != nout && (!n_filled || length(weight) + n_filled != nout)) {
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
  if (n_filled) w <- c(w, rep(1, n_filled))
  if (anyNA(w)) {
    cli_warn("weights contained NAs, which were set to 0")
    w[is.na(w)] <- 0
  }
  realign <- FALSE
  su <- tid %in% mtid
  if (!all(su)) {
    if (verbose) {
      cli_alert_info(paste(
        "{.field {sum(!su)}} of {.field {length(su)}} target IDs",
        "were dropped because they were not present in {.arg map}"
      ))
    }
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
  if (source_sf) st_geometry(source) <- NULL
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
    cli_alert_info(paste(
      "redistributing {.field {length(source_numeric)}} variable{?s}",
      "from {.field {length(sid)}} source{?s} to {.field {length(tid)}} target{?s}:"
    ))
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
    type <- if (aggregate) "aggregating" else "disaggregating"
    cli_progress_step(paste0(type, "..."), msg_done = paste("done", type))
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
  if (verbose) cli_progress_done()
  if (!all(source_numeric)) {
    if (verbose) cli_alert_info("re-converting categorical levels")
    for (l in level_map) {
      x <- res[, l$index]
      x[x == 0] <- NA
      res[, l$index] <- l$levels[x]
    }
  }
  if (rescale) {
    if (verbose) {
      status <- "checking totals"
      cli_progress_step("{status}")
    }
    source_totals <- vapply(
      seq_along(source_numeric),
      function(i) if (source_numeric[i]) sum(source[, i], na.rm = TRUE) else 0,
      0
    )
    res_totals <- vapply(
      seq_along(source_numeric),
      function(i) if (source_numeric[i]) sum(res[, i], na.rm = TRUE) else 0,
      0
    )
    su <- which(source_totals != res_totals)
    if (length(su)) {
      if (verbose) {
        status <- paste0("rescaling ", length(su), " variable", if (length(su) == 1) "" else "s")
        cli_progress_update()
      }
      res_totals[res_totals == 0] <- 1
      for (i in su) res[, i] <- res[, i] / res_totals[i] * source_totals[i]
      if (verbose) {
        status <- sub("rescaling", "rescaled", status, fixed = TRUE)
        cli_progress_done()
      }
    } else if (verbose) {
      status <- "totals are aligned"
      cli_progress_done()
    }
  }
  res <- cbind(id = tid, res)
  if (realign) {
    if (verbose) cli_alert_info("realigning with original target IDs")
    rownames(res) <- tid
    res <- res[otid, ]
    res$id <- otid
    if (!is.na(default_value)) res[!res$id %in% tid, -1] <- default_value
    rownames(res) <- NULL
  }
  if (!is.null(outFile)) {
    if (verbose) cli_progress_step("writing results to {.file {outFile}}")
    dir.create(dirname(outFile), FALSE, TRUE)
    if (!grepl("\\.\\w", outFile)) outFile <- paste0(outFile, ".csv")
    vroom_write(res, outFile, ",")
    if (verbose) cli_progress_done()
  }
  if (return_geometry && target_sf) {
    st_geometry(res) <- if (n_filled) {
      c(st_geometry(target), if (source_sf) {
        structure(source_geom, names = sid)[missed_sources]
      } else {
        st_sfc(lapply(seq_len(n_filled), function(i) st_point()), crs = st_crs(target))
      })
    } else {
      st_geometry(target)
    }
  }
  invisible(res)
}
