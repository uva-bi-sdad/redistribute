#' Redistribute Data
#'
#' Distribute data from a source frame to a target frame.
#'
#' @param source A matrix-like object you want to distribute from; usually this will be
#' the real or more complete dataset, and is often at a lower resolution / higher level.
#' @param target A matrix-like object you want to distribute to: usually this will be
#' the dataset you want but isn't available, and is often at a higher resolution / lower level.
#' @param map A list with entries named with \code{source} IDs, containing vectors of associated
#' \code{target} IDs.
#' @param source_id,target_id Name of a column in \code{source} / \code{target},
#' or a vector containing IDs.
#' @param target_weight Name of a column in \code{target}, or a vector containing target weights.
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
#' @useDynLib redistribute, .registration = TRUE
#' @export

redistribute <- function(source, target, map = list(), source_id = "GEOID", target_id = source_id,
                         target_weight = "population", outFile = NULL, overwrite = FALSE, verbose = FALSE) {
  if (!overwrite && !is.null(outFile) && file.exists(outFile)) {
    cli_abort("{.arg outFile} already exists; use {.code overwrite = TRUE} to overwrite it")
  }
  if (is.null(dim(source))) source <- t(source)
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
    if (verbose) cli_alert_info("source IDs: {.field {colnames(source)[1]}} column of {.arg source}")
    sid <- source[, 1, drop = TRUE]
    source <- source[, -1, drop = FALSE]
  }
  sid <- as.character(sid)
  vars <- colnames(source)
  if (length(sid) != nrow(source)) cli_abort("{.arg source_id} is not the same length as {.field nrow(source)}")
  if (anyDuplicated(sid)) cli_abort("{.arg source_id} contains duplicates")
  if (!is.null(dim(target)) && is.null(colnames(target))) colnames(target) <- paste0("V", seq_len(ncol(target)))
  if (length(target_id) > 1) {
    if (verbose) cli_alert_info("target IDs: {.arg target_id} vector")
    tid <- target_id
  } else if (!is.null(dim(target)) && target_id %in% colnames(target)) {
    if (verbose) cli_alert_info("target IDs: {.field {target_id}} column of {.arg target}")
    tid <- target[, target_id, drop = TRUE]
    target <- target[, colnames(target) != target_id]
  } else {
    if (is.null(dim(target))) {
      if (verbose) cli_alert_info("target IDs: {.arg target} vector")
      tid <- target
    } else {
      if (verbose) cli_alert_info("target IDs: {.field {colnames(target)[1]}} column of {.arg target}")
      tid <- target[, 1, drop = TRUE]
      target <- target[, -1, drop = FALSE]
    }
  }
  tid <- as.character(tid)
  if (length(map)) {
    if (verbose) cli_alert_info("map: provided list")
    if (!any(sid %in% names(map))) cli_abort("no source IDs were present in the ID map")
  } else {
    if (nrow(source) == 1) {
      if (verbose) cli_alert_info("map: all target IDs for single source")
      map[[sid]] <- tid
    } else if (all(nchar(sid) == nchar(sid[1])) && nchar(tid[1]) > nchar(sid[1])) {
      if (verbose) cli_alert_info("map: first {.field {nchar(sid[1])}} character{?s} of target IDs")
      map <- split(tid, substr(tid, 1, nchar(sid[1])))
    } else {
      cli_abort("no map was provided, and could not make one from IDs")
    }
  }
  map <- lapply(sid, function(id) {
    e <- map[[id]]
    if (length(e)) e else character()
  })
  tw <- if (length(target_weight) > 1) {
    if (verbose) cli_alert_info("target weights: provided vector")
    if (length(target_weight) != length(tid)) cli_abort("{.arg target_weight} is not the same length as {.arg target_id}")
    target_weight
  } else if (!is.null(dim(target)) && target_weight %in% colnames(target)) {
    if (verbose) cli_alert_info("target weights: {.field {target_weight}} column of {.arg target}")
    target[, target_weight, drop = TRUE]
  } else if (is.null(dim(target))) {
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
    cli_bullets(structure(
      paste(ifelse(source_numeric, "(numb) ", "(char)"), colnames(source)),
      names = rep("*", ncol(source))
    ))
  }
  res <- process_distribute(as.matrix(source), as.integer(source_numeric), tid, tw, map)
  res <- as.data.frame(res)
  colnames(res) <- vars
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
