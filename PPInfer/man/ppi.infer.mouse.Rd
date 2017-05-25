\name{ppi.infer.mouse}
\alias{ppi.infer.mouse}

\title{
Inferring functionally related proteins using protein networks for mouse
}
\description{
This function is designed for mouse protein-protein interaction from STRING database. Default format is 'mgi'. The number of proteins is 10 in default. Note that the number of proteins used as a target may be different from the number of proteins in the input since mapping between formats is not always one-to-one in getBM.
}
\usage{
ppi.infer.mouse(target, kernel, top = 10, cross = 0,
                input = "mgi_symbol", output = "mgi_symbol",
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
number of top proteins most closely related to target class (default: 10)
}
  \item{cross}{
if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model
}

  \item{input}{
  input format
}

  \item{output}{
  output format
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

\author{
Dongmin Jung, Xijin Ge
}




\seealso{
net.infer, getBM
}
\examples{
string.db.10090 <- STRINGdb$new(version='10',species=10090,
score_threshold = 999)
string.db.10090.graph <- string.db.10090$get_graph()
K.10090 <- net.kernel(string.db.10090.graph)
rownames(K.10090) <- substring(rownames(K.10090),7)
colnames(K.10090) <- substring(colnames(K.10090),7)
target <- colnames(K.10090)[1:100]
infer.mouse <- ppi.infer.mouse(target,K.10090,
input="ensembl_peptide_id")
}

