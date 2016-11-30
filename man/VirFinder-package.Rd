\name{VirFinder-package}
\alias{VirFinder-package}
\alias{VirFinder}
\docType{package}
\title{
  An R package for predicting viral sequences in a fasta file.
}
\description{
  The package provides functions to predict viral sequences in a fasta file, such as the assembled contigs from metagenomic data. The method has good prediction accuracy for short (~1kb) and noval viral sequences.
}
\details{
	The prediction is based on the sequence signatures (k-tuple word frequencies) that distinguish virus from host sequences. The model was trained using equal number of known viral and host sequences.
	For a query sequence, the number of occurrences of k-tuple words are first counted by a c++ program using a hash table.
	Then the sequence is predicted based on the k-tuple word frequencies using a logistic regression model trained with previously known sequences.

	The package contains,
	two functions: VF.pred, VF.qvalue, and

	two data: VF.mod8mer, LC.contigs.gt1000.top30.fa
}
\author{
Jie Ren, Nathan Ahlgren, Jed Fuhrman, Fengzhu Sun

Maintainer: Jie Ren <renj@usc.edu>
}
\references{
Ren J, Ahlgren N, Fuhrman J and Sun F (2017). VirFinder: identifying viral sequences from metagenomic data using sequence signatures
}
%\keyword{ package }
%\seealso{
%  Optional links to other man pages
%}
%\examples{
%  \dontrun{
%     ## Optional simple examples of the most important functions
%     ## These can be in \dontrun{} and \donttest{} blocks.
%  }
%}
