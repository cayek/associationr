library(associationr)
library(datasets)

context("bootstrap function testing")

test_that("test boot.lm on iris datasets", {
  iris.boot <- boot.lm(Sepal.Length ~ Petal.Length, iris)
  expect_is(iris.boot, "function")
})
