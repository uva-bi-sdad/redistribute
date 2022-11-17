test_that("minimal example works", {
  expect_identical(
    redistribute(1:3, 1:5), data.frame(id = as.character(1:5), V1 = 1 / 5, V2 = 2 / 5, V3 = 3 / 5)
  )
})

test_that("from data, with weights works", {
  source <- data.frame(a = 1, b = 2, c = 3)
  target <- data.frame(id = 1:5, weight = c(1, 1, 2, 2, 5))
  w <- target$weight / sum(target$weight)
  expect_equal(
    redistribute(source, target),
    data.frame(id = as.character(1:5), a = 1 * w, b = 2 * w, c = 3 * w)
  )
})

test_that("multiple source IDs, with categorical variables works", {
  source <- data.frame(id = c("a", "b"), cat = c("aaa", "bbb"), num = c(1, 2))
  target <- data.frame(id = sample(paste0(c("a", "b"), rep(1:5, 2))), population = sample.int(1e5, 10))
  res <- redistribute(source, target)
  map <- split(target$id, substr(target$id, 1, 1))
  mres <- do.call(rbind, lapply(names(map), function(l) {
    source_row <- source[source$id == l, ]
    tar <- target[target$id %in% map[[l]], ]
    data.frame(id = tar$id, cat = source_row$cat, num = source_row$num * tar$population / sum(tar$population))
  }))
  rownames(mres) <- mres$id
  mres <- mres[target$id, ]
  rownames(mres) <- NULL
  expect_equal(redistribute(source, target), mres)
})

test_that("tall works", {
  source <- data.frame(id = 1:5, a = rnorm(5), b = rnorm(5))
  source_tall <- data.frame(
    id = rep(source$id, 2),
    variable = rep(c("a", "b"), each = 5),
    value = c(source$a, source$b)
  )
  target <- data.frame(id = paste0(rep(1:5, 6), 1:40), population = sample.int(1e5, 40))
  expect_equal(redistribute(source_tall, target), redistribute(source, target))
})

test_that("intersect map work", {
  source <- sf::st_as_sf(data.frame(id = c("a", "b"), a = rnorm(2), b = rnorm(2), c = rnorm(2), geometry = c(
    "POLYGON ((.5 .5, 1 -1, 1 1, .5 .5))", "POLYGON ((0 0, .1 -.1, .1 .1, 0 0))"
  )), wkt = "geometry")
  target <- data.frame(id = sample(paste0(c("a", "b"), rep(1:5, 2))), population = sample.int(1e5, 10))
  sf::st_geometry(target) <- sf::st_geometry(source)[grepl("^b", target$id) + 1]
  map <- sf::st_intersects(source, target)
  expect_identical(substring(capture.output(
    res <- redistribute(source[, -1], target, verbose = TRUE),
    type = "message"
  ), 3), c(
    "source IDs: sequence, assuming map from geometries",
    "target IDs: sequence, assuming map from geometries",
    "map: intersections between geometries",
    "weights: population column of `target`",
    "disaggregating 3 variables:",
    "(numb) a, b, c"
  ))
  expect_identical(redistribute(source, target, source_id = "id")[, -1], res[, -1])
  expect_identical(redistribute(source, target, map = map, source_id = "id")[, -1], res[, -1])
})

test_that("aggregation works", {
  source_original <- data.frame(id = c("a", "b"), cat = c("aaa", "bbb"), num = c(1, 2))
  target_original <- data.frame(id = sample(paste0(c("a", "b"), rep(1:5, 2))), population = sample.int(1e5, 10))
  source <- redistribute(source_original, target_original)
  expect_equal(redistribute(source, target = c("a", "b")), source_original)
})
