set.seed(123)
x <- rprobMat(10,c(2,3,4))

test_that("rprobMat creates tables object of distributions correctly", {
  expect_equal(dim(x), c(10, 24))
  expect_equal(tdim(x), c(2,3,4))
  expect_equal(apply(x,1,sum), rep(1,10))
})

test_that("[ operator for tables works correctly", {
  expect_equal(dim(x[1,]), c(2,3,4))
  expect_equal(dim(x[1,,keep=TRUE]), c(1,24))
  expect_equal(dim(x[,1,3,2:3]), c(10,2))  
  expect_equal(tdim(x[,1,3,2:3]), c(2))  
  expect_equal(dim(x[,1,3,2:3,drop=FALSE]), c(10,2))  
  expect_equal(tdim(x[,1,3,2:3,drop=FALSE]), c(1,1,2))  
})

