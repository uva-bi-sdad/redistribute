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
    redistribute(source, target, weight = "weight"),
    data.frame(id = as.character(1:5), a = 1 * w, b = 2 * w, c = 3 * w)
  )
})

test_that("multiple source IDs, with categorical variables works", {
  source <- data.frame(id = c("a", "b"), cat = c("aaa", "bbb"), num = c(1, 2))
  target <- data.frame(id = sample(paste0(c("a", "b"), rep(1:5, 2))), population = sample.int(1e5, 10))
  res <- redistribute(source, target, weight = "population")
  map <- split(target$id, substr(target$id, 1, 1))
  mres <- do.call(rbind, lapply(names(map), function(l) {
    source_row <- source[source$id == l, ]
    tar <- target[target$id %in% map[[l]], ]
    data.frame(id = tar$id, cat = source_row$cat, num = source_row$num * tar$population / sum(tar$population))
  }))
  rownames(mres) <- mres$id
  mres <- mres[target$id, ]
  rownames(mres) <- NULL
  expect_equal(redistribute(source, target, weight = "population"), mres)
})

test_that("tall works", {
  source <- data.frame(id = 1:5, a = rnorm(5), b = rnorm(5))
  source_tall <- data.frame(
    id = rep(source$id, 2),
    variable = rep(c("a", "b"), each = 5),
    value = c(source$a, source$b)
  )
  target <- data.frame(id = paste0(rep(1:5, 6), 1:40), population = sample.int(1e5, 40))
  expect_equal(
    redistribute(source_tall, target, weight = "population"),
    redistribute(source, target, weight = "population")
  )
})

test_that("aggregation works", {
  source_original <- data.frame(id = c("a", "b"), cat = c("aaa", "bbb"), num = c(1, 2))
  target_original <- data.frame(id = sample(paste0(c("a", "b"), rep(1:5, 2))), population = sample.int(1e5, 10))
  source <- redistribute(source_original, target_original, weight = "population")
  expect_equal(redistribute(source, target = c("a", "b")), source_original)
})

test_that("intersect map work", {
  library(sf)
  square <- matrix(c(0, 1, 1, 0, 0, 0, 0, 1, 1, 0), 5)
  adj_up <- matrix(c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1), 5)
  adj_right <- matrix(c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0), 5)
  source <- st_as_sf(data.frame(
    id = c("a", "b"), a = c(225, 2250), b = rnorm(2) * 1000, c = rnorm(2) * 1000,
    geometry = st_sfc(st_polygon(list(square * 5e-4)), st_polygon(list(square * 5e-4 + adj_right * 5e-4)))
  ))
  vars <- colnames(source)[2:4]
  target <- st_as_sf(data.frame(
    id = paste0(rep(c("a", "ab", "b"), c(20, 5, 20)), c(1:20, 1:5, 1:20)), population = 1000, # sample.int(1e5, 45),
    geometry = do.call(c, lapply(0:8, function(col) {
      st_sfc(lapply(0:4, function(r) {
        st_polygon(list(round(square + adj_right * .5 + adj_up * r + adj_right * col, 2) * 1e-4))
      }))
    }))
  ))

  ## disaggregation
  map <- redistribute(
    source, target,
    source_id = "id", target_id = "id", return_map = TRUE, make_intersect_map = TRUE
  )
  mres <- target
  for (v in vars) mres[, v] <- 0
  mres <- mres[, c("id", vars, "population")]
  rownames(mres) <- mres$id
  for (id in names(map)) {
    r <- st_drop_geometry(source[source$id == id, ])
    e <- map[[id]]
    tids <- names(e)
    w <- mres[tids, "population", drop = TRUE] * e
    for (v in vars) {
      mres[tids, v] <- mres[tids, v, drop = TRUE] + r[[v]] * w / sum(w)
    }
  }
  rownames(mres) <- NULL
  expect_identical(sub(" \\[.*$", "", substring(capture.output(
    res <- redistribute(
      source, target,
      weight = "population", source_id = "id", target_id = "id", map = map, verbose = TRUE
    ),
    type = "message"
  ), 3)), c(
    "source IDs: id column of `source`",
    "target IDs: id column of `target`",
    "map: provided list",
    "weights: population column of `target`",
    "redistributing 3 variables from 2 sources to 45 targets:",
    "(numb; 3) a, b, c",
    "disaggregating...",
    "done disaggregating",
    ""
  ))
  expect_equal(res[, 1:4, drop = TRUE], mres[, 1:4, drop = TRUE])
  expect_identical(redistribute(source[, -1], target, weight = "population")[, -1], res[, -1])

  ## aggregation
  map_weight <- c(1, .5)[grepl("^ab", target$id) + 1]
  mares <- source
  mares[, 2:4] <- 0
  for (v in vars) {
    for (id in source$id) {
      su <- grepl(paste0("^(?:", id, "|ab)"), target$id)
      mares[mares$id == id, v] <- sum(mres[su, v, drop = TRUE] * map_weight[su])
    }
  }
  agg_map <- redistribute(
    target[, -2], source,
    source_id = "id", target_id = "id", return_map = TRUE, make_intersect_map = TRUE
  )
  expect_equal(redistribute(mres[, 1:4], source, map = agg_map, source_id = "id", target_id = "id"), mares)
  expect_equal(redistribute(mres[, 1:4], source, source_id = "id", target_id = "id", make_intersect_map = TRUE), mares)
})
