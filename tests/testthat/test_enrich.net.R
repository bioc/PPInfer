

test_that("enrich.net()", {
  data(examplePathways)
  data(exampleRanks)
  set.seed(1)
  result.GSEA <- fgsea(examplePathways, exampleRanks, nperm = 1000)
  
  plot.enrich.net <- tempfile("plot", fileext = ".png")
  png(filename = plot.enrich.net, width = 2000, height = 1600, res = 300)
  enrich.net(result.GSEA, examplePathways, node.id = 'pathway',
             pvalue = 'pval', edge.cutoff = 0.6, degree.cutoff = 1,
             n = 50, vertex.color = 'red', vertex.label.cex = 0.75,
             show.legend = FALSE, layout=igraph::layout.kamada.kawai)
  dev.off()
})

