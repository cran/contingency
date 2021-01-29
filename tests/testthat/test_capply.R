set.seed(127)
x <- rprobMat(100, dim=c(2,3,4), alpha=2)
x2 <- capply(x, margin, c(1,3))

test_that("capply works on tables", {
  expect_equal(c(margin(x, c(1,3))), c(x2))
  expect_equal(capply(x, entropy), apply(as.array(x), 1, entropy))
  expect_equal(capply(x, mutualInf, m1=1, m2=3, condition=2), 
               apply(as.array(x), 1, mutualInf, m1=1,m2=3,condition=2))
})
