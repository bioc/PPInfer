

test_that("ORA.dotplot()", {
  data(examplePathways)
  data(exampleRanks)
  geneNames <- names(exampleRanks)
  set.seed(1)
  gene.id <- sample(geneNames, 100)
  result.ORA <- ORA(examplePathways, gene.id)
  
  plot.ORA <- tempfile("plot", fileext = ".png")
  png(filename = plot.ORA, width = 2000, height = 1600, res = 300)
  ORA.dotplot(result.ORA, category = "Category", size = "Size",
              count = "Count", pvalue = "pvalue", sort = "pvalue")
  dev.off()
  
  expect_error(ORA.dotplot(result.ORA, size = "Size", count = "Count",
                           pvalue = "pvalue", sort = "pvalue"))
  expect_error(ORA.dotplot(result.ORA, category = "Category", count = "Count",
                           pvalue = "pvalue", sort = "pvalue"))
  expect_error(ORA.dotplot(result.ORA, category = "Category", size = "Size",
                           pvalue = "pvalue", sort = "pvalue"))
  
  result.table <- ORA.dotplot(result.ORA, category = "Category", size = "Size",
                              count = "Count", pvalue = "pvalue", sort = "pvalue",
                              plot = FALSE)
  expect_is(result.table$Pvalue, 'numeric')
  
  expect_equivalent(result.table, ORA.dotplot(result.ORA, category = "Category", size = "Size",
                                              count = "Count", pvalue = "pvalue", sort = "pvalue",
                                              plot = FALSE))
  expect_equal(result.table, ORA.dotplot(result.ORA, category = "Category", size = "Size",
                                         count = "Count", pvalue = "pvalue", sort = "pvalue",
                                         plot = FALSE))
  expect_identical(result.table, ORA.dotplot(result.ORA, category = "Category", size = "Size",
                                             count = "Count", pvalue = "pvalue", sort = "pvalue",
                                             plot = FALSE))
})

