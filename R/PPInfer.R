




# kernel matrix
net.kernel=function(g, decay=0.5)
{
  L <- as.matrix(laplacian_matrix(g,normalized=TRUE))
  I <- diag(1,dim(L))
  K <- solve(I + decay*L)
  K
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





### inference for human
ppi.infer.human <- function(target, kernel, top=10, cross=0,
                            input='hgnc_symbol', output='hgnc_symbol',
                            C1 = 1, nu = 0.2, epsilon = 0.1,
                            cache1 = 40, tol1 = 0.001, shrinking1 = TRUE,
                            C2 = 1, cache2 = 40, tol2 = 0.001,
                            shrinking2 = TRUE)
{
  ensembl <- useMart("ensembl")
  human.ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
  
  # input
  new.list <-
    getBM(attributes=c('ensembl_peptide_id',input),
          filters=input, values = target,
          mart = human.ensembl)[,1]
  new.list <- na.omit(new.list)
  
  # main
  ppi.pred.9606 <- net.infer(new.list, kernel, cross=cross,
                             C1 = C1, nu = nu, epsilon = epsilon,
                             cache1 = cache1, tol1 = tol1,
                             shrinking1 = shrinking1, C2 = C2,
                             cache2 = cache2, tol2 = tol2,
                             shrinking2 = shrinking2)
  
  # output
  protein_score <- data.frame(ppi.pred.9606$top,ppi.pred.9606$score)
  format.protein <-
    getBM(attributes=c('ensembl_peptide_id',output),
          filters='ensembl_peptide_id',
          values = ppi.pred.9606$top,
          mart = human.ensembl)
  index <- match(protein_score[,1],format.protein[,1])
  new.protein_score <- cbind(protein_score,format.protein[index,2])
  new.protein_score[new.protein_score==""]<-NA
  na.omit.new.protein_score <- na.omit(new.protein_score)
  ppi.pred.9606$top <- as.vector(na.omit.new.protein_score[1:top,3])
  ppi.pred.9606$score <- as.numeric(na.omit.new.protein_score[1:top,2])
  
  ppi.pred.9606
}



### inference for mouse
ppi.infer.mouse <- function(target, kernel, top=10, cross=0,
                            input='mgi_symbol', output='mgi_symbol',
                            C1 = 1, nu = 0.2, epsilon = 0.1,
                            cache1 = 40, tol1 = 0.001, shrinking1 = TRUE,
                            C2 = 1, cache2 = 40, tol2 = 0.001,
                            shrinking2 = TRUE)
{
  ensembl <- useMart("ensembl")
  mouse.ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
  
  # input
  new.list <-
    getBM(attributes=c('ensembl_peptide_id',input),
          filters=input,values = target,
          mart = mouse.ensembl)[,1]
  new.list <- na.omit(new.list)
  
  # main
  ppi.pred.10090 <- net.infer(new.list,kernel,cross=cross,
                              C1 = C1, nu = nu, epsilon = epsilon,
                              cache1 = cache1, tol1 = tol1,
                              shrinking1 = shrinking1,
                              C2 = C2, cache2 = cache2, 
                              tol2 = tol2, shrinking2 = shrinking2)
  
  # output
  protein_score <- data.frame(ppi.pred.10090$top,ppi.pred.10090$score)
  format.protein <-
    getBM(attributes=c('ensembl_peptide_id',output),
          filters='ensembl_peptide_id',
          values = ppi.pred.10090$top,
          mart = mouse.ensembl)
  index <- match(protein_score[,1],format.protein[,1])
  new.protein_score <- cbind(protein_score,format.protein[index,2])
  new.protein_score[new.protein_score==""]<-NA
  na.omit.new.protein_score <- na.omit(new.protein_score)
  ppi.pred.10090$top <- as.vector(na.omit.new.protein_score[1:top,3])
  ppi.pred.10090$score <- as.numeric(na.omit.new.protein_score[1:top,2])
  
  ppi.pred.10090
}


