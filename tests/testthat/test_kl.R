set.seed(123)
x <- rprobMat(10,c(2,3,4,5))
y <- x[4,]

test_that("KL-divergence works", {
  expect_equal(kl(x, y), apply(x, 1, function(z) kl(z,y)))
  expect_equal(kl(x, y)[4], 0)
})

