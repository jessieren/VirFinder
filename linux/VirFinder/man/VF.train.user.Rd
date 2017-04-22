\name{VF.train.user}
\alias{VF.train.user}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Train virus prediction model using user's database
}
\description{
The function trains the viral contig prediction model using user's customized database 
containing viral sequences and prokaryotic sequences. 

The sequences are first splitted into fragments of fixed lengths, such as 0.5 kb, 1 kb and 3kb.
, and the k-tuple word frequency are counted for each fragment.
The prediction model is trained based on the k-tuple frequencies of viral and prokaryotic fragments. 

}
\usage{
VF.train.user(trainFaFileHost, trainFaFileVirus, userModDir, userModName, w, equalSize)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainFaFileHost}{
The fasta file (including path) of the prokaryotic sequences. 
}
  \item{trainFaFileVirus}{
The fasta file (including path) of the viral sequences.
}
  \item{userModDir}{
The directory where the trained model will be saved.
}
  \item{userModName}{
The name of the model file.
}
  \item{w}{
The length of the k-tuple word.
}
  \item{equalSize}{
If use the same number of virus and host fragments to train? Default is FALSE.
}
}

\value{
The function returns the trained model. There models are trained based on three different fragment lengths, 
0.5 kb, 1 kb and 3 kb, respectively. The model for 0.5 kb will be used for predicting query sequences of length < 1 kb,
the model for 1 kb will be used for predicting sequences of length ranging from 1 kb to 3 kb, and the model for 3 kb will be used for predicting 
sequences of length > 3kb. The model is trained based on k-tuple frequencies, which can be specified using w.

The trained model will also be saved into the user's specified directory, 
for easily loading and calling it in the future.

}

\note{
  One must be careful that the sequences in the database can be trusted as truly viral and prokaryotic, respectively. 
  The assembled contigs in metagenome could possibly have viral sequences, and similarity the contigs in metavirome could also
  have host sequences. 
  
  The longer k-tuple can describe better the difference between virus and host sequences, 
but if the data is not enough, it can lead to an overfitting problem. 
Given the database, users are suggested to test different lengths of k-tuple in order to get the best model.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

## (1) train the model using user's database
#### (1.1) specifiy the fasta files of the training contigs, one for virus and one for prokaryotic hosts
trainFaFileHost <- system.file("data", "tara_host.fa", package="VirFinder")
trainFaFileVirus <- system.file("data", "tara_virus.fa", package="VirFinder")

#### (1.2) specify the directory where the trained model will be saved to, and the name of the model
userModDir <- file.path(find.package("VirFinder"))
userModName <- "modTara"

## (2) train the model using user's database
w <- 4  # the length of the k-tuple word
VF.trainModUser <- VF.train.user(trainFaFileHost, trainFaFileVirus, userModDir, 
userModName, w, equalSize=TRUE)

## (3) load the trained model based on user's database
modFile <- list.files(userModDir, userModName, full.names=TRUE)
load(modFile)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ customization }
