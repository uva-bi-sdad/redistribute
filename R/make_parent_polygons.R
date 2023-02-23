#' Make Higher-Order Polygons
#'
#' Randomly combine given polygons into contiguous, larger polygons.
#'
#' @param polygons An \code{sf} object with polygons to be reformed.
#' @param ids A vector of IDs the same length as \code{polygons}, or the column name
#' in \code{polygons} to use as IDs.
#' @param n Number of polygons to aim to combine for each new polygon.
#' @param strict_n Logical; if \code{TRUE} (default), \code{n} represents the number of
#' intersecting polygons to select, depending on availability (with a minimum of 2).
#' If \code{FALSE}, \code{n} represents the number of nearest polygons (by centroid)
#' to use when calculating a box to use when selecting polygons (with a minimum of 1).
#' @param n_as_min Logical; if \code{TRUE}, will merge any parents with fewer than \code{n}
#' children with a random neighbor. Otherwise (and by default), parents may have fewer than \code{n}
#' children. Applies if \code{strict_n} is \code{TRUE}.
#' @param buffer_dist Distance around each initial shape or set of shapes, used to define
#' neighboring shapes. Applies if \code{strict_n} is \code{TRUE}
#' @param min_overlap Minimal area of overlap between potential neighboring shapes and
#' the buffered target shape, used to define neighboring shapes.
#' Applies if \code{strict_n} is \code{TRUE}
#' @param verbose Logical; if \code{FALSE}, will not print status messages.
#' @returns A list with entries for \code{new} (an \code{sf} object containing the new polygons)
#' and \code{map} (a list mapping between the old and new polygon indices).
#' @export

make_parent_polygons <- function(
    polygons, ids = NULL, n = 2, strict_n = TRUE, n_as_min = FALSE, buffer_dist = 5e-5,
    min_overlap = NULL, verbose = TRUE) {
  if (!inherits(polygons, "sf")) cli_abort("{.arg polygons} must be an sf object")
  if (length(ids) == 1 && ids %in% colnames(polygons)) ids <- polygons[[ids]]
  polygons <- st_geometry(polygons)
  polygons <- st_cast(st_boundary(polygons), "POLYGON", do_split = FALSE)
  if (strict_n && is.null(NULL)) min_overlap <- min(s2_area(polygons)) * .0035
  n_poly <- length(polygons)
  if (!is.null(ids)) {
    if (length(ids) != n_poly) cli_abort("{.arg ids} do not align with {.arg polygons}")
    names(polygons) <- ids
  }
  n <- max(if (strict_n) 2 else 1, min(n, n_poly - 1))
  sels <- seq_len(n) + 1L
  centers <- st_coordinates(suppressWarnings(st_centroid(polygons)))
  rownames(centers) <- seq_len(nrow(centers))
  map <- new <- NULL
  accounted <- integer(n_poly)
  if (verbose) {
    cli_progress_step(
      "mapping child polygons ({.field {sum(accounted)}}/{.field {n_poly}})",
      msg_done = paste0(
        "created initial map from {.field {length(map)}} parent polygons to ",
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
      sg <- st_buffer(polygons[[start]], buffer_dist, 1)
      ag <- polygons[available]
      pre <- st_intersects(sg, ag)[[1]]
      if (length(pre)) {
        pre <- pre[vapply(ag[pre], function(t) {
          s2_area(st_intersection(t, sg, model = "closed"))
        }, 0) > min_overlap]
      }
      set <- unique(c(start, available[pre]))
      sn <- length(set)
      if (sn < n) {
        available <- available[!available %in% set]
        while (length(set) < n && length(available)) {
          sg <- st_buffer(st_boundary(polygons[set])[[1]], buffer_dist, 1)
          ag <- polygons[available]
          pre <- st_intersects(sg, ag)[[1]]
          if (length(pre)) {
            pre <- pre[vapply(ag[pre], function(t) {
              s2_area(st_intersection(t, sg, model = "closed"))
            }, 0) > min_overlap]
          }
          set <- unique(c(set, available[pre]))
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
  map <- data.frame(
    super = rep(seq_along(map), vapply(map, length, 0)),
    sub = unlist(map, use.names = FALSE)
  )
  map <- map[!duplicated(map$sub), ]
  map <- split(map$sub, map$super)
  if (n_as_min) {
    su <- which(vapply(map, length, 0) < n)
    if (length(su)) {
      i <- 0
      defunct <- NULL
      if (verbose) {
        cli_progress_step(
          "merging in {.field {length(defunct)}}/{.field {length(su)}} parent{?/s} with too few children",
          msg_done = "merged in {.field {length(defunct)}} of {.field {length(su)}} small parent{?/s}",
          spinner = TRUE
        )
      }
      tmap <- data.frame(
        super = rep(seq_along(map), vapply(map, length, 0)),
        sub = unlist(map, use.names = FALSE)
      )
      available <- seq_along(polygons)
      for (i in su) {
        available <- available[!available %in% map[[i]]]
        sg <- st_buffer(st_boundary(polygons[map[[i]]])[[1]], buffer_dist, 1)
        ag <- polygons[available]
        pre <- st_intersects(sg, ag)[[1]]
        if (length(pre)) {
          pre <- pre[vapply(ag[pre], function(t) {
            s2_area(st_intersection(t, sg, model = "closed"))
          }, 0) > min_overlap]
        }
        if (length(pre)) {
          ns <- available[if (length(pre) == 1) pre else sample(pre, 1)]
          np <- tmap[[1]][tmap[[2]] == ns]
          map[[np]] <- c(map[[np]], map[[i]])
          available <- available[!available %in% map[[np]]]
          defunct <- c(defunct, i)
        }
        if (verbose) cli_progress_update()
      }
      if (length(defunct)) map <- map[-defunct]
      if (verbose) cli_progress_done()
    }
    map <- data.frame(
      super = rep(seq_along(map), vapply(map, length, 0)),
      sub = unlist(map, use.names = FALSE)
    )
    map <- map[!duplicated(map$sub), ]
    map <- split(map$sub, map$super)
  }
  if (!is.null(ids)) map <- lapply(map, function(l) ids[l])
  i <- 0
  if (verbose) {
    cli_progress_step(
      "creating new polygons ({.field {i}}/{.field {length(map)}})",
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
