

test_that("net.infer()", {
  data(litG)
  litG <- graph_from_graphnel(litG)
  sg <- decompose(litG, min.vertices = 50)
  sg <- sg[[1]]
  K <- net.kernel(sg)
  expect_error(net.infer(names(V(sg))[1], K, top=20))
  expect_error(net.infer(names(V(sg)), K, top=20))
  
  Net.infer <- net.infer(names(V(sg))[1:10], K, top=20)
  expect_is(Net.infer, 'list')
  expect_null(Net.infer$CVerror)
  
  expect_equivalent(Net.infer, net.infer(names(V(sg))[1:10], K, top=20))
  expect_equal(Net.infer, net.infer(names(V(sg))[1:10], K, top=20))
  expect_identical(Net.infer, net.infer(names(V(sg))[1:10], K, top=20))
  
  expect_equal(Net.infer$list, names(V(sg))[1:10])
  expect_true(identical(intersect(net.infer(names(V(sg))[1:10], K, top=10)$top, names(V(sg))[1:10]), character(0)))
  expect_gt(length(Net.infer$top), length(net.infer(names(V(sg))[1:10], K, top=10)$top))
})

