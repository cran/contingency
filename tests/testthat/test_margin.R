set.seed(123)
x <- rprobMat(10,c(2,3,4,2))

test_that("conditioning works for object of class 'tables'", {
  expect_equal(conditional2(x,3:4,1:2), x/c(margin(x,1:2)))
  expect_equal(conditional(x,3:4,1:2,condition.value=list(1:2,1)), 
               aperm((x/c(margin(x,1:2)))[,1:2,1,,,drop=FALSE],c(3,4,1,2)))
  #expect_equal(c(intervention(x,3,1:2)), c(t((apply(as.array(x), 1, intervention, 3, 1:2)))))
})

test_that("intervention works for object of class 'tables'", {
  expect_equal(intervention(x,3,1:2), c(margin(x,1:2))*conditional2(x,4,1:3))
  expect_equal(c(intervention(x,3,1:2)), c(t((apply(as.array(x), 1, intervention, 3, 1:2)))))
})

