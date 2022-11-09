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
