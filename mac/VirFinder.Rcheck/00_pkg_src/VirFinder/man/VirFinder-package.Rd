\name{VirFinder-package}
\alias{VirFinder-package}
\alias{VirFinder}
\docType{package}
\title{
  An R package for predicting viral sequences in a fasta file.
}
\description{
  The package provides functions to predict viral sequences in a fasta file, such as the assembled contigs from metagenomic data. 
  The method has good prediction accuracy for short (~1kb) and novel viral sequences.
}
\details{
	The prediction is based on the sequence signatures (k-tuple word frequencies) that distinguish virus from host sequences. 
  The model was trained using all known viral sequences obtained from complete virus refSeqs and the equal number of host sequences in NCBI.
	For a query sequence, the number of occurrences of k-tuple words are first counted by a c++ program using a hash table.
	Then the sequence is predicted based on the k-tuple word frequencies using a logistic regression model trained with previously known sequences.
  
  Users can also train the prediction model using their own database of choice. 
  One must be careful that the sequences in the database can be trusted as truly viral and prokaryotic, respectively. 

	The package contains the following functions: "VF.pred", "VF.qvalue", "VF.train.user", and "VF.pred.user".
  "VF.pred" predicts viral sequences using the defualt trained model. 
  "VF.train.user" and "VF.pred.user" trains and then predicts viral sequences using the user's customized database.

	The package contains five data files: "VF.trainMod8mer.rda", "contigs.fa", "crAssphage.fa", "tara_virus.fa", and "tara_host.fa". 
  "VF.trainMod8mer.rda" is the default model trained using all known viral sequences and equal number of host sequences in NCBI.
  "contigs.fa" is a testing example with 30 contigs asembled from a human gut metagenome. 
  "crAssphage.fa" is the genomic sequence of crAssphage, a newly discovered virus ubiquitous in healthy human gut.
  "tara_virus.fa" and "tara_host.fa" are small subsets of viral and prokaryotic contigs respectively, assembled using samples from Tara Ocean Expeditions. 
  The two fasta file are used for testing the function of training model with customized database.
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
