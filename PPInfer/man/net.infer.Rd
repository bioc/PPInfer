\name{net.infer}
\alias{net.infer}

\title{
Inferring functionally related proteins using networks
}
\description{
The OCSVM and classical SVM are sequentially used. First, we apply the OCSVM by training a one-class classifier using the data from the known class only and classify the background data. Let n be the number of proteins in the target class. This model is used to predict remaining N-n proteins in the background. Proteins with zero similarity with the target class or presence data are extracted. Then they are potentially defined as the other class by pseudo-absence selection methods from spatial statistics. For balanced data, it is assumed that two classes contain the same number of proteins. Next, by the classical SVM, these two classes can be used to infer remaining N-2n proteins and to find proteins closely related to the known class or interesting proteins. Note that the number of proteins is not greater than a half of total proteins in a kernel matrix due to N-2n>0. Also, the number of top proteins to be inferred is less than or equal to N-2n.
}
\usage{
net.infer(target, kernel, top = NULL, cross = 0,
          C1 = 1, nu = 0.2, epsilon = 0.1, cache1 = 40,
          tol1 = 0.001, shrinking1 = TRUE, C2 = 1,
          cache2 = 40, tol2 = 0.001, shrinking2 = TRUE)
}

\arguments{
  \item{target}{
list of interesting proteins or target class
}
  \item{kernel}{
the regularized Laplacian matrix for a graph
}
  \item{top}{
number of top proteins most closely related to target class  (default: all proteins except for target and pseudo-absence class)
}
  \item{cross}{
if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model
}
  \item{C1}{
cost of constraints violation (default: 1) for OCSVM
}
  \item{nu}{
The nu parameter for OCSVM (default: 0.2)
}
  \item{epsilon}{
epsilon in the insensitive-loss function for OCSVM (default: 0.1)
}
  \item{cache1}{
cache memory in MB for OCSVM (default: 40)
}
  \item{tol1}{
tolerance of termination criterion for OCSVM (default: 0.001)
}
  \item{shrinking1}{
option whether to use the shrinking-heuristics for OCSVM (default: TRUE)
}
  \item{C2}{
cost of constraints violation (default: 1) for SVM
}
  \item{cache2}{
cache memory in MB for SVM (default: 40)
}
  \item{tol2}{
tolerance of termination criterion for SVM (default: 0.001)
}
  \item{shrinking2}{
option whether to use the shrinking-heuristics for SVM (default: TRUE)
}
}


\value{
\item{list}{
list of a target class used in the model
}

\item{error}{
training error
}

\item{CVerror}{
cross validation error, (when cross > 0)
}

\item{top}{
top proteins
}

\item{score}{
decision values for top proteins
}
}

\references{
Senay, S. D. et al. (2013). Novel three-step pseudo-absence selection technique for improved species distribution modelling. PLOS ONE. 8(8), e71218.
}
\author{
Dongmin Jung, Xijin Ge
}



\seealso{
ksvm
}
\examples{
# example 1
\dontrun{
string.db.9606 <- STRINGdb$new(version='10',species=9606,
score_threshold = 999)
string.db.9606.graph <- string.db.9606$get_graph()
K.9606 <- net.kernel(string.db.9606.graph)
rownames(K.9606) <- substring(rownames(K.9606),6)
colnames(K.9606) <- substring(colnames(K.9606),6)
target <- colnames(K.9606)[1:100]
infer <- net.infer(target,K.9606,10)
}

# example 2
data(litG)
litG <- graph_from_graphnel(litG)
sg <- decompose(litG, min.vertices=50)
sg <- sg[[1]]
K <- net.kernel(sg)
litG.infer <- net.infer(names(V(sg))[1:10],K,top=20)
}


