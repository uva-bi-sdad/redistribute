library(sf)

test_that("basic example works", {
  square <- matrix(c(0, 1, 1, 0, 0, 0, 0, 1, 1, 0), 5)
  adj_up <- matrix(c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1), 5)
  adj_right <- matrix(c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0), 5)
  original <- st_as_sf(data.frame(
    geometry = do.call(c, lapply(0:8, function(col) {
      st_sfc(lapply(0:4, function(r) {
        st_polygon(list(round(square + adj_right * .5 + adj_up * r + adj_right * col, 2) * 1e-4))
      }))
    }))
  ))
  parents <- NULL
  expect_message(parents <- make_parent_polygons(original, n = 3), "mapping")
  expect_true(all(vapply(parents$map, length, 0) < 4))
  expect_identical(st_bbox(original), st_bbox(parents$new))
  parents <- make_parent_polygons(original, strict_n = FALSE)
  expect_identical(sort(unlist(parents$map, use.names = FALSE)), seq_len(nrow(original)))
})
