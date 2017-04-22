\name{VF.pred}
\alias{VF.pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify viral sequences in a fasta file
}
\description{
The function identifies viral sequences in the input fasta file using the trained model.
A score (higher value suggests higher possibility) and a p-value (lower value suggests higher possibility) will be assigned to each query sequence.
For a query sequence of length L, if L < 1 kb, the model trained using 0.5 kb sequence is used for prediction;
if 1 kb <= L < 3 kb, the model trained using 1000 bp sequence is used; if L>=3 kb, the model trained using 3000 bp sequences is used for prediction.
The prediction model uses the feature of k-tuple word frequency of a sequence to predict if it is a viral sequence. The k-tuple length used in the model is 8.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VF.pred(inFaFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inFaFile}{
The full file name (including path) of the fasta file
}
}

\value{
%%  ~Describe the value returned
The function returns a data frame. The rows correspond to sequences, and the columns are from the left to the right, sequence name (name), sequence length (length), prediction score (score), and prediction p-value (pvalue).
}

\examples{

## (1) set the input fasta file name. 
library(VirFinder)
inFaFile <- system.file("data", "contigs.fa", package="VirFinder")

## (2) prediction
predResult <- VF.pred(inFaFile)
predResult

## sort sequences by p-value in ascending order
predResult[order(predResult$pvalue),]

## (3) predict for crAssphage
inFaFile <- system.file("data", "crAssphage.fa", package="VirFinder")
VF.pred(inFaFile)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ default }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
