

test_that("ORA()", {
  data(examplePathways)
  data(exampleRanks)
  geneNames <- names(exampleRanks)
  set.seed(1)
  gene.id <- sample(geneNames, 100)
  
  expect_error(ORA(examplePathways))
  expect_error(ORA(gene.id))

  
  result.ORA <- ORA(examplePathways, gene.id)
  expect_is(result.ORA$pvalue, 'numeric')
  
  expect_equivalent(result.ORA, ORA(examplePathways, gene.id))
  expect_equal(result.ORA, ORA(examplePathways, gene.id))
  expect_identical(result.ORA, ORA(examplePathways, gene.id))
})

