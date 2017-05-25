

test_that("net.kernel()", {
  data(litG)
  litG <- graph_from_graphnel(litG)
  sg <- decompose(litG, min.vertices = 50)
  sg <- sg[[1]]
  K <- net.kernel(sg)
  expect_is(K, "matrix")
  
  n <- dim(K)[1]
  expect_equal(n, length(diag(K)))
})

