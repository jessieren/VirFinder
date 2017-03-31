\name{VF.pred.user}
\alias{VF.pred.user}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify viral sequences in a fasta file using user's trained prediction model
}
\description{
The function identifies viral sequences in the input fasta file using the model trained based on the user's database. 
A score (higher value suggests higher possibility) and a p-value (lower value suggests higher possibility) will be assigned to each query sequence.
For a query sequence of length L, if L < 1 kb, the model trained using 0.5 kb sequence is used for prediction;
if 1 kb <= L < 3 kb, the model trained using 1 kb sequence is used; if L>=3 kb, the model trained using 3000 bp sequences is used for prediction.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VF.pred.user(inFaFile, VF.trainModUser)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inFaFile}{
The full file name (including path) of the fasta file
}
  \item{VF.trainModUser}{
The trained model using user's database
}
}

\value{
%%  ~Describe the value returned
The function returns a data frame. The rows correspond to sequences, and the columns are from the left to the right, sequence name (name), sequence length (length), prediction score (score), and prediction p-value (pvalue).
}

\note{
To train models using user's database of viral and prokaryotic sequences, please use the function "VF.train.user".
}

\examples{

## (1) specifiy the fasta files of the training contigs
#### (1.1) one for virus and one for prokaryotic hosts
trainFaFileHost <- system.file("data", "tara_host.fa", package="VirFinder")
trainFaFileVirus <- system.file("data", "tara_virus.fa", package="VirFinder")

#### (1.2) specify the directory where the trained model will be saved to, and the name of the model
userModDir <- file.path(find.package("VirFinder"))
userModName <- "modTara"

## (2) train the model using user's database
w <- 6  # the length of the k-tuple word
VF.trainModUser <- VF.train.user(trainFaFileHost, trainFaFileVirus, userModDir, userModName, w)

## (3) predict the contigs using the customized model
#### (3.1) specify the fasta file containing contigs for prediction
inFaFile <- system.file("data", "contigs.fa", package="VirFinder")

#### (3.2) prediction
predResultUser <- VF.pred.user(inFaFile, VF.trainModUser)
predResultUser

#### (3.3) sort sequences by p-value in ascending order
predResultUser[order(predResultUser$pvalue),]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ customization }
