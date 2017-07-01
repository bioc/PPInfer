

test_that("net.infer.ST()", {
  data(litG)
  litG <- graph_from_graphnel(litG)
  sg <- decompose(litG, min.vertices = 50)
  sg <- sg[[1]]
  K <- net.kernel(sg)
  expect_error(net.infer.ST(names(V(sg))[1], K, top=20))
  expect_error(net.infer.ST(names(V(sg)), K, top=20))
  
  Net.infer.ST <- net.infer.ST(names(V(sg))[1:10], K, top=20)
  expect_is(Net.infer.ST, 'list')
  
  expect_equivalent(Net.infer.ST, net.infer.ST(names(V(sg))[1:10], K, top=20))
  expect_equal(Net.infer.ST, net.infer.ST(names(V(sg))[1:10], K, top=20))
  expect_identical(Net.infer.ST, net.infer.ST(names(V(sg))[1:10], K, top=20))
  
  expect_equal(Net.infer.ST$list, names(V(sg))[1:10])
  expect_true(identical(intersect(net.infer.ST(names(V(sg))[1:10], K, top=10)$top, names(V(sg))[1:10]), character(0)))
  expect_gt(length(Net.infer.ST$top), length(net.infer.ST(names(V(sg))[1:10], K, top=10)$top))
})

