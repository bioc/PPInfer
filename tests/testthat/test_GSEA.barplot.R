

test_that("GSEA.barplot()", {
  data(examplePathways)
  data(exampleRanks)
  set.seed(1)
  result.GSEA <- fgsea(examplePathways, exampleRanks, nperm = 1000)
  
  plot.GSEA <- tempfile("plot", fileext = ".png")
  png(filename = plot.GSEA, width = 2000, height = 1600, res = 300)
  GSEA.barplot(result.GSEA, category = 'pathway', score = 'NES', pvalue = 'pval',
               sort = 'NES', decreasing = TRUE)
  dev.off()
  
  expect_error(GSEA.barplot(result.GSEA, score = 'NES', pvalue = 'pval',
                            sort = 'NES', decreasing = TRUE))
  expect_error(GSEA.barplot(result.GSEA, category = 'pathway', pvalue = 'pval',
                            sort = 'NES', decreasing = TRUE))
  expect_error(GSEA.barplot(result.GSEA, category = 'pathway', score = 'NES',
                            sort = 'NES', decreasing = TRUE))
  
  result.table <- GSEA.barplot(result.GSEA, category = 'pathway', score = 'NES',
                               pvalue = 'pval', sort = 'NES', decreasing = TRUE,
                               plot = FALSE)
  expect_is(result.table$pval, 'numeric')
  
  expect_equivalent(result.table, GSEA.barplot(result.GSEA, category = 'pathway', score = 'NES',
                                               pvalue = 'pval', sort = 'NES', decreasing = TRUE,
                                               plot = FALSE))
  expect_equal(result.table, GSEA.barplot(result.GSEA, category = 'pathway', score = 'NES',
                                        pvalue = 'pval', sort = 'NES', decreasing = TRUE,
                                        plot = FALSE))
  expect_identical(result.table, GSEA.barplot(result.GSEA, category = 'pathway', score = 'NES',
                                            pvalue = 'pval', sort = 'NES', decreasing = TRUE,
                                            plot = FALSE))
})

