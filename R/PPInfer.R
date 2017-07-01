




# kernel matrix
net.kernel <- function (g, decay = 0.5) 
{
  L <- as.matrix(laplacian_matrix(g, normalized = TRUE))
  colname <- colnames(L)
  rowname <- rownames(L)
  I <- diag(1, dim(L))
  K <- chol2inv(chol(I + decay * L))
  colnames(K) <- colname
  rownames(K) <- rowname
  return(K)
}



# inference for network
net.infer <- function(target, kernel, top=NULL, cross=0,
                      C1 = 1, nu = 0.2, epsilon = 0.1, cache1 = 40,
                      tol1 = 0.001, shrinking1 = TRUE, C2 = 1,
                      cache2 = 40, tol2 = 0.001, shrinking2 = TRUE)
{
  # find new list
  node <- rownames(kernel)
  n <- length(node)
  if(length(target)<=1) stop("size of list is too small")
  index <- match(target,node)!='NA'
  new.target <- target[!is.na(index)]
  n1 <- length(na.omit(index))
  if(n-2*n1<=0) stop("size of list is too large")
  
  
  #### OCSVM
  
  # train
  new.index <- !is.na(match(node,new.target))
  train.kernel <- kernel[new.index,new.index]
  model <- ksvm(train.kernel, type="one-svc",kernel="matrix",
                C = C1, nu = nu, epsilon = epsilon, cache = cache1,
                tol = tol1, shrinking = shrinking1)
  # predict
  test.kernel <- as.kernelMatrix(kernel[,new.index][,SVindex(model), drop=FALSE]) 
  pred.decision <- predict(model, test.kernel,type='decision')
  
  
  # assign class 1
  y <- ifelse(new.index==TRUE,1,'NA')
  
  pred <- cbind(node,pred.decision,y)
  sort.pred <- pred[order(as.numeric(pred[,2])),]
  sort.pred2 <- subset(sort.pred,sort.pred[,3]=='NA')
  new.sort.pred2 <- sort.pred2[order(as.numeric(sort.pred2[,2])),]
  target0 <- as.vector(new.sort.pred2[1:n1,1])
  
  # assign pseudo class 0
  new.index <- !is.na(match(node,target0))
  y <- ifelse(new.index==TRUE,0,y)
  
  
  
  #### SVM
  
  # train
  model <- ksvm(kernel[y!="NA", y!="NA"], as.numeric(y[y!="NA"]),
                type="C-svc", C=C2, kernel="matrix", cross=cross,
                cache = cache2, tol = tol2, shrinking = shrinking2)

  
  # predict
  test.kernel <- as.kernelMatrix(kernel[, y!="NA"][,SVindex(model), drop=FALSE]) 
  pred.decision <- predict(model, test.kernel,type='decision')
  
  
  
  
  pred <- cbind(node,pred.decision,y)
  colnames(pred)[2] <- 'decision values'
  
  sort.pred <- pred[order(as.numeric(pred[,2]),decreasing = TRUE),]
  sort.pred2 <- subset(sort.pred,sort.pred[,3]=='NA')
  new.sort.pred2 <- sort.pred2[order(as.numeric(sort.pred2[,2]),
                                     decreasing = TRUE),]
  
  
  if( is.null(top) ) {top=n-2*n1}  
  if( (n-2*n1)<top ) {top=n-2*n1}
  pred.svm <- list()
  pred.svm$list <- new.target
  pred.svm$error <- error(model)
  if (cross > 0) {pred.svm$CVerror=cross(model)}
  pred.svm$top <- as.vector(new.sort.pred2[1:top,1])
  pred.svm$score <- as.numeric(new.sort.pred2[1:top,2])
  
  pred.svm
}





# self training
self.train.kernel <- function (K, y, type = 'response', C = 1, 
                               cache = 40, tol = 0.001,
                               shrinking = TRUE, thrConf = 0.9,
                               maxIts = 10, percFull = 1, verbose = FALSE) 
{
  if(dim(K)[1]!=dim(K)[2])
    stop(" kernel matrix not square!")
  
  learner = "ksvm"
  learner.pars = list(type = "C-svc", prob.model = TRUE,
                      kernel = "matrix", C = C, cache = cache,
                      tol = tol, shrinking = shrinking)
  pred <- function(m,d)
  {
    p <- predict(m, as.kernelMatrix(d), type="prob")
    data.frame(cl = colnames(p)[apply(p, 1, which.max)],
               p = apply(p, 1, max))
  }
  pred.pars = list()
  
  data <- data.frame(y, K)
  N <- nrow(data)
  it <- 0
  sup <- which(!is.na(y))
  repeat
  {
    it <- it + 1
    model <- do.call(learner, c(list(K[sup, sup], data[sup, 1]), 
                                learner.pars))
    test.kernel <- as.kernelMatrix(K[,sup][, SVindex(model), drop = FALSE])
    pred.type <- predict(model, test.kernel, type = type)
    probPreds <- do.call(pred, c(list(model, test.kernel), pred.pars))
    new <- which(probPreds[, 2] > thrConf)
    if (verbose) 
      cat("IT.", it, "\t nr. added exs. =", length(new), "\n")
    if (length(new))
    {
      index <- setdiff(new, sup)
      sup <- unique(c(sup, new))
      data[index, 1] <- probPreds[index, 1]
    }
    else break
    if (it == maxIts || length(sup)/N >= percFull) 
      break
  }
  return(pred.type)
}





# inference for network with self training
net.infer.ST <- function(target, kernel, top = NULL, C1 = 1, nu = 0.2, 
                         epsilon = 0.1, cache1 = 40, tol1 = 0.001, shrinking1 = TRUE, 
                         C2 = 1, cache2 = 40, tol2 = 0.001, shrinking2 = TRUE,
                         thrConf = 0.9, maxIts = 10, percFull = 1, verbose = FALSE) 
{
  node <- rownames(kernel)
  n <- length(node)
  if (length(target) <= 1) 
    stop("size of list is too small")
  index <- match(target, node) != "NA"
  new.target <- target[!is.na(index)]
  n1 <- length(na.omit(index))
  if (n - 2 * n1 <= 0) 
    stop("size of list is too large")
  
  
  #### OCSVM
  new.index <- !is.na(match(node, new.target))
  train.kernel <- kernel[new.index, new.index]
  model <- ksvm(train.kernel, type = "one-svc", kernel = "matrix", 
                C = C1, nu = nu, epsilon = epsilon, cache = cache1, tol = tol1, 
                shrinking = shrinking1)
  test.kernel <- as.kernelMatrix(kernel[, new.index][, SVindex(model), 
                                                     drop = FALSE])
  pred.decision <- predict(model, test.kernel, type = "decision")
  y <- ifelse(new.index == TRUE, 1, "NA")
  pred <- cbind(node, pred.decision, y)
  sort.pred <- pred[order(as.numeric(pred[, 2])), ]
  sort.pred2 <- subset(sort.pred, sort.pred[, 3] == "NA")
  new.sort.pred2 <- sort.pred2[order(as.numeric(sort.pred2[, 
                                                           2])), ]
  target0 <- as.vector(new.sort.pred2[1:n1, 1])
  new.index <- !is.na(match(node, target0))
  y <- ifelse(new.index == TRUE, 0, y)
  
  # self training
  pred.decision  <- self.train.kernel(kernel, factor(y, levels=c(0,1)), type = 'decision',
                                      C = C2, cache = cache2, tol = tol2, shrinking = shrinking2,
                                      thrConf = thrConf, maxIts = maxIts,
                                      percFull = percFull, verbose = verbose)
  
  pred <- cbind(node, pred.decision, y)
  colnames(pred)[2] <- "decision values"
  sort.pred <- pred[order(as.numeric(pred[, 2]), decreasing = TRUE), 
                    ]
  sort.pred2 <- subset(sort.pred, sort.pred[, 3] == "NA")
  new.sort.pred2 <- sort.pred2[order(as.numeric(sort.pred2[, 
                                                           2]), decreasing = TRUE), ]
  if (is.null(top)) {
    top = n - 2 * n1
  }
  if ((n - 2 * n1) < top) {
    top = n - 2 * n1
  }
  pred.svm <- list()
  pred.svm$list <- new.target
  pred.svm$error <- error(model)
  
  pred.svm$top <- as.vector(new.sort.pred2[1:top, 1])
  pred.svm$score <- as.numeric(new.sort.pred2[1:top, 2])
  pred.svm
}






### inference for human
ppi.infer.human <- function (target, kernel, top = 10, classifier = net.infer,
                             input = "hgnc_symbol", output = "hgnc_symbol", ...) 
{
  ensembl <- useMart("ensembl")
  human.ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
  
  # input
  new.list <- getBM(attributes = c("ensembl_peptide_id", input), 
                    filters = input, values = target, mart = human.ensembl)[, 1]
  new.list <- na.omit(new.list)
  
  # main
  ppi.pred.9606 <- classifier(new.list, kernel, ...)
  
  # output
  protein_score <- data.frame(ppi.pred.9606$top, ppi.pred.9606$score)
  format.protein <- getBM(attributes = c("ensembl_peptide_id", output),
                          filters = "ensembl_peptide_id", values = ppi.pred.9606$top, 
                          mart = human.ensembl)
  index <- match(protein_score[, 1], format.protein[, 1])
  new.protein_score <- cbind(protein_score, format.protein[index, 2])
  new.protein_score[new.protein_score == ""] <- NA
  na.omit.new.protein_score <- na.omit(new.protein_score)
  ppi.pred.9606$top <- as.vector(na.omit.new.protein_score[1:top, 3])
  ppi.pred.9606$score <- as.numeric(na.omit.new.protein_score[1:top, 2])
  ppi.pred.9606
}





### inference for mouse
ppi.infer.mouse <- function (target, kernel, top = 10, classifier = net.infer,
                             input = "mgi_symbol", output = "mgi_symbol", ...) 
{
  ensembl <- useMart("ensembl")
  mouse.ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
  
  # input
  new.list <- getBM(attributes = c("ensembl_peptide_id", input), 
                    filters = input, values = target, mart = mouse.ensembl)[,1]
  new.list <- na.omit(new.list)
  
  # main
  ppi.pred.10090 <- classifier(new.list, kernel, ...)
  
  # output
  protein_score <- data.frame(ppi.pred.10090$top, ppi.pred.10090$score)
  format.protein <- getBM(attributes = c("ensembl_peptide_id", output),
                          filters = "ensembl_peptide_id", values = ppi.pred.10090$top, 
                          mart = mouse.ensembl)
  index <- match(protein_score[, 1], format.protein[, 1])
  new.protein_score <- cbind(protein_score, format.protein[index, 2])
  new.protein_score[new.protein_score == ""] <- NA
  na.omit.new.protein_score <- na.omit(new.protein_score)
  ppi.pred.10090$top <- as.vector(na.omit.new.protein_score[1:top, 3])
  ppi.pred.10090$score <- as.numeric(na.omit.new.protein_score[1:top, 2])
  ppi.pred.10090
}



