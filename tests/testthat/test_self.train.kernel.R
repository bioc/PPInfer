

test_that("self.train.kernel()", {
  data(litG)
  litG <- graph_from_graphnel(litG)
  sg <- decompose(litG, min.vertices=50)
  sg <- sg[[1]]
  K <- net.kernel(sg)
  y <- rep(NA, length(V(sg)))
  y[1:10] <- 1
  y[11:20] <- 0
  y <- factor(y)
  
  expect_is(self.train.kernel(K, y), "factor")
  expect_is(self.train.kernel(K, y, type = 'probabilities'), "matrix")
  expect_is(self.train.kernel(K, y, type = 'votes'), "matrix")
  expect_is(self.train.kernel(K, y, type = 'decision'), "matrix")
  
  expect_false(identical(levels(self.train.kernel(K, y)), unique(y)))
  expect_true(identical(levels(self.train.kernel(K, y)), levels(y)))
  
  expect_length(self.train.kernel(K, y), length(V(sg)))
  expect_error(self.train.kernel(matrix(rnorm(9),3,3), factor(c(1,0,NA,0,1))))
  expect_error(self.train.kernel(matrix(rnorm(20),5,4), factor(c(1,0,NA,0,1))))
})

