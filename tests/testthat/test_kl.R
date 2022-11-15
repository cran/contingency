set.seed(123)
x <- rprobMat(10,c(2,3,4,5))
y <- x[4,]

test_that("KL-divergence works", {
  expect_equal(kl(x, y), apply(x, 1, function(z) kl(z,y)))
  expect_equal(kl(x, y)[4], 0)
})

x0 <- x
x0[] <- 1
for (i in seq_along(tdim(x))) {
  x0 <- x0*margin2(x, i)
  if (i == 3) x1 <- x0[,,,,1]
}

test_that("multiinformation works", {
  expect_equal(multiInf(x), entropy(x/x0))
  expect_equal(multiInf(x, 1:3), entropy(margin(x, 1:3)/x1))
})
