#' Redistribute Data
#'
#' Distribute data from a source frame to a target frame.
#'
#' @param source A matrix-like object you want to distribute from; usually this will be
#' the real or more complete dataset, and is often at a lower resolution / higher level.
#' @param target A matrix-like object you want to distribute to: usually this will be
#' the dataset you want but isn't available, and is often at a higher resolution / lower level.
#' @param map A list with entries named with \code{source} IDs (or aligning with those IDs),
#' containing vectors of associated \code{target} IDs (or indices of those IDs). If IDs
#' are related by substrings (the first characters of \code{target} IDs are \code{source} IDs),
#' then a map can be automatically generated from them. If \code{source} and \code{target}
#' contain \code{sf} geometries, a map will be made with \code{\link[sf]{st_intersects}}
#' (\code{st_intersects(source, target)}).
#' @param source_id,target_id Name of a column in \code{source} / \code{target},
#' or a vector containing IDs. For \code{source}, this will default to the first column. For
#' \code{target}, columns will be searched through for one that appears to relate to the
#' source IDs, falling back to the first column.
#' @param target_weight Name of a column in \code{target}, or a vector containing target weights.
#' Defaults to unit weights (a weight of 1 for each \code{target} entry).
#' @param source_variable,source_value If \code{source} is tall (with variables spread across
#' rows rather than columns), specifies names of columns in \code{source} containing variable names
#' and values for conversion.
#' @param outFile Path to a CSV file in which to save results.
#' @param overwrite Logical; if \code{TRUE}, will overwrite an existing \code{outFile}.
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
#' @returns A \code{data.frame} the same dimensions as \code{target}, with an initial column
#' (\code{id}) containing \code{target_ids}, and additional columns from \code{source}.
#' @importFrom cli cli_abort cli_bullets cli_alert_info
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom utils write.csv
#' @importFrom sf st_intersects st_geometry<-
#' @useDynLib redistribute, .registration = TRUE
#' @export

redistribute <- function(source, target, map = list(), source_id = "GEOID", target_id = source_id,
                         target_weight = NULL, source_variable = NULL, source_value = NULL,
                         outFile = NULL, overwrite = FALSE, verbose = FALSE) {
  if (!overwrite && !is.null(outFile) && file.exists(outFile)) {
    cli_abort("{.arg outFile} already exists; use {.code overwrite = TRUE} to overwrite it")
  }
  source_sf <- inherits(source, "sf")
  target_sf <- inherits(target, "sf")
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
  tall <- !is.null(source_variable) || !is.null(source_value)
  if (tall || anyDuplicated(sid)) {
    source_numeric <- vapply(seq_len(ncol(source)), function(col) is.numeric(source[, col, drop = TRUE]), TRUE)
    if (is.null(source_value) && any(source_numeric)) source_value <- colnames(source)[which(source_numeric)[1]]
    if (is.null(source_variable) && any(!source_numeric)) {
      source_variable <- colnames(source)[
        which.min(vapply(which(!source_numeric), function(col) length(unique(source[, col])), 0))
      ]
    }
    if (is.null(source_variable) || is.null(source_value)) {
      cli_abort("{.arg source_id} contains duplicates, and {.arg source} does not appear to be tall")
    } else {
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
  }
  if (length(dim(target)) == 2 && is.null(colnames(target))) colnames(target) <- paste0("V", seq_len(ncol(target)))
  if (length(target_id) > 1) {
    if (verbose) cli_alert_info("target IDs: {.arg target_id} vector")
    tid <- target_id
  } else if (!is.null(target_id) && length(dim(target)) == 2 && target_id %in% colnames(target)) {
    if (verbose) cli_alert_info("target IDs: {.field {target_id}} column of {.arg target}")
    tid <- target[, target_id, drop = TRUE]
    target <- target[, colnames(target) != target_id]
  } else {
    if (length(dim(target)) != 2) {
      if (is.null(target_weight) && is.numeric(target) && !is.null(names(target))) {
        if (verbose) cli_alert_info("target IDs: names of {.arg target} vector; {.arg target_weight} set to {.arg target}")
        target_weight <- unname(target)
        tid <- names(target)
      } else {
        if (verbose) cli_alert_info("target IDs: {.arg target} vector")
        tid <- target
      }
    } else {
      if (target_sf && source_sf) {
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
  if (length(map)) {
    if (verbose) cli_alert_info("map: provided list")
    if (is.null(names(map))) {
      if (length(map) != length(sid)) cli_abort("{.arg map} has no names, and is not the same length as source IDs")
      names(map) <- sid
    } else if (!any(sid %in% names(map))) cli_abort("no source IDs were present in the ID map")
  } else {
    if (nrow(source) == 1) {
      if (verbose) cli_alert_info("map: all target IDs for single source")
      map[[sid]] <- tid
    } else if (!intersect_map && all(nchar(sid) == nchar(sid[1])) && nchar(tid[1]) > nchar(sid[1])) {
      if (verbose) cli_alert_info("map: first {.field {nchar(sid[1])}} character{?s} of target IDs")
      map <- split(tid, substr(tid, 1, nchar(sid[1])))
    } else if (source_sf && target_sf) {
      if (verbose) cli_alert_info("map: intersections between geometries")
      op <- options(sf_use_s2 = FALSE)
      on.exit(options(sf_use_s2 = op[[1]]))
      map <- tryCatch(st_intersects(source, target), error = function(e) NULL)
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
  map <- lapply(sid, function(id) {
    e <- map[[id]]
    if (length(e)) {
      if (is.integer(e)) tid[e] else e
    } else {
      character()
    }
  })
  tw <- if (length(target_weight) > 1) {
    if (verbose) cli_alert_info("target weights: {.arg target_weight} vector")
    if (length(target_weight) != length(tid)) cli_abort("{.arg target_weight} is not the same length as {.arg target_id}")
    target_weight
  } else if (!is.null(target_weight) && length(dim(target)) == 2 && target_weight %in% colnames(target)) {
    if (verbose) cli_alert_info("target weights: {.field {target_weight}} column of {.arg target}")
    target[, target_weight, drop = TRUE]
  } else if (length(dim(target)) != 2) {
    if (verbose) cli_alert_info("target weights: {.field 1}")
    rep(1, length(tid))
  } else {
    su <- vapply(seq_len(ncol(target)), function(col) is.numeric(target[, col, drop = TRUE]), TRUE)
    if (any(su)) {
      target_weight <- colnames(target)[which(su)[1]]
      if (verbose) cli_alert_info("target weights: {.field {target_weight}} column of {.arg target}")
      target[, target_weight, drop = TRUE]
    } else {
      if (verbose) cli_alert_info("target weights: {.field 1}")
      rep(1, length(tid))
    }
  }
  tw <- as.numeric(tw)
  if (anyDuplicated(tid)) {
    agger <- if (all(tw == as.integer(tw))) "summing" else "averaging"
    if (verbose) cli_alert_info("{.arg target_id} contains duplicates, so {agger} {.arg target_weight}s")
    tw <- tapply(tw, tid, if (agger == "summing") sum else mean)
    tid <- names(tw)
    tw <- unname(tw)
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
    cli_alert_info("redistributing {.field {length(source_numeric)}} variable{?s}:")
    var_groups <- split(colnames(source), !source_numeric)
    names(var_groups) <- c("TRUE" = "(char)", "FALSE" = "(numb)")[names(var_groups)]
    cli_bullets(structure(
      vapply(names(var_groups), function(type) paste(type, paste(var_groups[[type]], collapse = ", ")), ""),
      names = rep("*", length(var_groups))
    ))
  }
  res <- process_distribute(as.matrix(source), as.integer(source_numeric), tid, tw, map)
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
  if (!is.null(outFile)) {
    if (verbose) cli_alert_info("writing results to {.file {outFile}}")
    dir.create(dirname(outFile), FALSE, TRUE)
    if (!grepl("\\.\\w", outFile)) outFile <- paste0(outFile, ".csv")
    write.csv(res, outFile, row.names = FALSE)
  }
  invisible(res)
}
