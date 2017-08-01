

test_that("ppi.infer.mouse() gives correct errors", {
  string.db.10090 <- STRINGdb$new(version='10', species = 10090,
                                  score_threshold = 999)
  string.db.10090.graph <- string.db.10090$get_graph()
  K.10090 <- net.kernel(string.db.10090.graph)
  rownames(K.10090) <- substring(rownames(K.10090), 7)
  colnames(K.10090) <- substring(colnames(K.10090), 7)
  target <- colnames(K.10090)
  
  expect_is(ppi.infer.mouse(target[1:100], K.10090,
                            input="ensembl_peptide_id"), 'list')
  expect_error(ppi.infer.mouse(target[1], K.10090,
                               input="ensembl_peptide_id"))
  expect_error(ppi.infer.mouse(target, K.10090,
                               input="ensembl_peptide_id"))
  expect_error(ppi.infer.mouse(target[1:100], K.10090,
                               input="mgi_symbol"))
  
  expect_is(ppi.infer.mouse(target[1:100], K.10090,
                            classifier = net.infer.ST,
                            input="ensembl_peptide_id"), 'list')
  expect_error(ppi.infer.mouse(target[1], K.10090,
                               classifier = net.infer.ST,
                               input="ensembl_peptide_id"))
  expect_error(ppi.infer.mouse(target, K.10090,
                               classifier = net.infer.ST,
                               input="ensembl_peptide_id"))
  expect_error(ppi.infer.mouse(target[1:100], K.10090,
                               classifier = net.infer.ST,
                               input="mgi_symbol"))
})

