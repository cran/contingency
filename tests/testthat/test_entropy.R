set.seed(123)
x <- rprobMat(10,c(2,3,4,5))
#x[c(3,19,26,112)] = 0
x2 <- as.array(x)
x3 <- rprobdist(c(2,3,4,5))

test_that("entropy works", {
  expect_equal(entropy(x), apply(x, 1, entropy))
  expect_equal(entropy(x[1,]), -sum(x[1,]*log(x[1,])))
})

test_that("mutual information works", {
  expect_equal(mutualInf(x,1,2), entropy(x,1)+entropy(x,2)-entropy(x,1:2))
  expect_equal(mutualInf(x,1,2,3), entropy(x,c(1,3))+entropy(x,2:3)-entropy(x,1:3)-entropy(x,3))
  expect_equal(mutualInf(x,1,2), apply(x2, 1, mutualInf, 1, 2))
  expect_equal(mutualInf(x,1,2,3), apply(x2, 1, mutualInf, 1, 2, 3))  
  expect_equal(mutualInf(x3,1,2), entropy(x3,1)+entropy(x3,2)-entropy(x3,1:2))
  expect_equal(mutualInf(x3,1,2,3), entropy(x3,c(1,3))+entropy(x3,2:3)-entropy(x3,1:3)-entropy(x3,3))
})

test_that("interaction information works", {
  expect_equal(interactionInf(x, 1, 2, 3), mutualInf(x,1,2,3) - mutualInf(x,1,2))
  expect_equal(interactionInf(x3, 1, 2, 3), mutualInf(x3,1,2,3) - mutualInf(x3,1,2))
  expect_equal(interactionInf(x, 1, 2, 3), apply(x2, 1, interactionInf,1,2,3))
  expect_equal(interactionInf(x, 1, 2, 3:4), mutualInf(x,1,2,3:4) - mutualInf(x,1,2)) 
  expect_equal(interactionInf(x3, 1, 2, 3:4), mutualInf(x3,1,2,3:4) - mutualInf(x3,1,2))
})
