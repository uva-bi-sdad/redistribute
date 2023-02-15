#' Make Higher-Order Polygons
#'
#' Randomly combine given polygons into contiguous, larger polygons.
#'
#' @param polygons An \code{sf} object with polygons to be reformed.
#' @param n Number of polygons to aim to combine for each new polygon.
#' @param strict_n Logical; if \code{TRUE} (default), \code{n} represents the number of
#' intersecting polygons to select, depending on availability (with a minimum of 2).
#' If \code{FALSE}, \code{n} represents the number of nearest polygons (by centroid)
#' to use when calculating a box to use when selecting polygons (with a minimum of 1).
#' @param verbose Logical; if \code{FALSE}, will not print status messages.
#' @returns A list with entries for \code{new} (an \code{sf} object containing the new polygons)
#' and \code{map} (a list mapping between the old and new polygon indices).
#' @export

make_parent_polygons <- function(polygons, n = 2, strict_n = TRUE, verbose = TRUE) {
  if (!inherits(polygons, "sf")) cli_abort("{.arg polygons} must be an sf object")
  polygons <- st_geometry(polygons)
  polygons <- st_cast(st_boundary(polygons), "POLYGON", do_split = FALSE)
  n_poly <- length(polygons)
  n <- max(if (strict_n) 2 else 1, min(n, n_poly - 1))
  sels <- seq_len(n) + 1L
  centers <- st_coordinates(suppressWarnings(st_centroid(polygons)))
  rownames(centers) <- seq_len(nrow(centers))
  map <- new <- NULL
  accounted <- integer(n_poly)
  if (verbose) {
    cli_progress_step(
      "mapping child polygons ({sum(accounted)}/{n_poly})",
      msg_done = paste0(
        "created map from {.field {length(map)}} parent polygons to ",
        "{.field {n_poly}} child polygons"
      ),
      spinner = TRUE
    )
  }
  for (i in seq_len(n_poly)) {
    if (verbose) cli_progress_update()
    available <- which(accounted == 0L)
    if (!length(available)) break
    start <- if (length(available) == 1) available else sample(available, 1)
    if (strict_n) {
      set <- unique(c(start, available[st_intersects(
        polygons[start], polygons[available]
      )[[1]]]))
      sn <- length(set)
      if (sn < n) {
        while (length(set) < n && length(available)) {
          set <- unique(c(set, available[st_intersects(
            polygons[set], polygons[available]
          )[[1]]]))
          if (length(set) > n) {
            set <- set[seq_len(n)]
            break
          } else {
            su <- !available %in% set
            if (all(su)) {
              break
            } else {
              available <- available[su]
            }
          }
        }
      } else {
        set <- set[seq_len(n)]
      }
    } else {
      nearest <- as.integer(names(sort(
        lma_simets(centers, centers[start, ], metric = "euclidean"), TRUE
      )[sels]))
      adj <- max(abs(centers[nearest, ] - rep(unname(centers[start, ]), each = n)))
      set <- unique(c(start, available[st_intersects(
        st_polygon(list(matrix(
          rep(centers[start, ], each = 5) - adj * c(1, -1, -1, 1, 1, 1, 1, -1, -1, 1), 5
        ))), polygons[available]
      )[[1]]]))
    }
    accounted[set] <- 1L
    map <- c(list(set), map)
  }
  i <- 0
  map <- data.frame(
    super = rep(seq_along(map), vapply(map, length, 0)),
    sub = unlist(map, use.names = FALSE)
  )
  map <- map[!duplicated(map$sub), ]
  map <- split(map$sub, map$super)
  if (verbose) {
    cli_progress_step(
      "creating new polygons ({i}/{length(map)})",
      msg_done = "created {.field {length(map)}} new polygons", spinner = TRUE
    )
  }
  new <- do.call(c, lapply(map, function(l) {
    if (verbose) {
      i <<- i + 1
      cli_progress_update(.envir = parent.frame(2))
    }
    st_union(polygons[l])
  }))
  list(new = new, map = map)
}
