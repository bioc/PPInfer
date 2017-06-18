

test_that("ppi.infer.human()", {
  string.db.9606 <- STRINGdb$new(version='10', species = 9606,
                                 score_threshold = 999)
  string.db.9606.graph <- string.db.9606$get_graph()
  K.9606 <- net.kernel(string.db.9606.graph)
  rownames(K.9606) <- substring(rownames(K.9606), 6)
  colnames(K.9606) <- substring(colnames(K.9606), 6)
  target <- colnames(K.9606)
  
  expect_is(ppi.infer.human(target[1:100], K.9606,
                            input="ensembl_peptide_id"), 'list')
  expect_error(ppi.infer.human(target[1], K.9606,
                               input="ensembl_peptide_id"))
  expect_error(ppi.infer.human(target, K.9606,
                               input="ensembl_peptide_id"))
  expect_error(ppi.infer.human(target[1:100], K.9606,
                               input="hgnc_symbol"))
  
  expect_is(ppi.infer.human(target[1:100], K.9606,
                            classifier = net.infer.ST,
                            input="ensembl_peptide_id"), 'list')
  expect_error(ppi.infer.human(target[1], K.9606,
                               classifier = net.infer.ST,
                               input="ensembl_peptide_id"))
  expect_error(ppi.infer.human(target, K.9606,
                               classifier = net.infer.ST,
                               input="ensembl_peptide_id"))
  expect_error(ppi.infer.human(target[1:100], K.9606,
                               classifier = net.infer.ST,
                               input="hgnc_symbol"))
})

