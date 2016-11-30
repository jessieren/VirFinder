# VirFinder: R package for identifying viral sequences from metagenomic data using sequence signatures
Version: 1.0

Author: Jie Ren, Nathan Ahlgren, Jed Fuhrman, Fengzhu Sun

Maintainer: Jie Ren <renj@usc.edu>

Description
----------------

The package provides functions to predict viral sequences in a fasta file, such as the assembled contigs from metagenomic data. The method has good prediction accuracy for short (~1kb) and noval viral sequences.

The prediction method is based on the sequence signatures (k-tuple word frequencies) that distinguish virus from host sequences. The model was trained using equal number of known viral and host sequences. For a query sequence, the number of occurrences of k-tuple words are first counted by a c++ program using a hash table. Then the sequence is predicted based on the k-tuple word frequencies using a logistic regression model trained with previously known sequences.

Please refer to VirFinder-manual.pdf for usage instruction.


Dependencies
---------------
R packages "glmnet", "qvalue" and "Rcpp" are needed to be installed before intall VirFinder.

To install "glmnet" and "Rcpp", start R and enter,
	
	install.packages("glmnet")
	install.packages("Rcpp")


To install qvalue, start R and enter,

	## try http:// if https:// URLs are not supported
	source("https://bioconductor.org/biocLite.R")
	biocLite("qvalue")




Installation
---------------
To install the R package VirFinder, follow the instuctions on http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages.

To quick start, first download the package file VirFinder_1.0.tar.gz/VirFinder_1.0.zip according to your operating system.

For Mac/Linux users, if you have a Graphic User Interfaces (GUI) of R, you fire up a R graphic window and type 

	install.packages("<path_to_the_file>/VirFinder_1.0.tar.gz", repos = NULL, type="source")

	library(VirFinder)


If you are not using GUI of R, you can install the package from the command line. Simply type the following to the command line,

	R CMD INSTALL <path_to_the_file>/VirFinder_1.0.tar.gz



For Windows users, if you have a Graphic User Interfaces (GUI) of R, you fire up a R graphic window and type 

	install.packages("<path_to_the_file>/VirFinder_1.0.zip", repos = NULL, type="source")

	library(VirFinder)


If you are not using GUI of R, you can install the package from the command line. Simply type the following to the command line,

	Rcmd INSTALL <path_to_the_file>\VirFinder_1.0.zip


