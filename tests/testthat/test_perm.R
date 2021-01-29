x <- rprobMat(10,c(2,3,4))

test_that("aperm works as expected", {
  expect_equal(aperm(x[1,], c(2,3,1)), aperm(x, c(2,3,1))[1,])
  expect_equal(aperm(x, c(1,2,3)), x)
})
