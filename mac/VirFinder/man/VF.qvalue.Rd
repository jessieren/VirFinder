\name{VF.qvalue}
\alias{VF.qvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the false discovery rates (q-values) using p-values
%%  ~~function to do ... ~~
}
\description{
The function estimates the false discovery rates (q-values) for a given set of p-values.
The q-value measures the proportion of false positives incurred when predicting viral sequences using the corresponding p-value as a threshold.
The function VF.qvalue uses the function "qvalue" in the R package "qvalue" by John D. Storey.
}
\usage{
VF.qvalue(pvalue, fdr.level = NULL, pfdr = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalue}{
%%     ~~Describe \code{pvalue} here~~
The predicted p-values for sequences obtained from the function "VF.pred"
}
  \item{fdr.level}{
A level at which to control the FDR. Must be in (0,1].
Optional; if this is selected, a vector of TRUE and FALSE is
returned that specifies whether each q-value is less than
fdr.level or not.
}
  \item{pfdr}{
An indicator of whether it is desired to make the estimate
more robust for small p-values and a direct finite sample
estimate of pFDR - optional.
}
  \item{\dots}{
Additional arguments passed to ‘pi0est’ and ‘lfdr’.
}
}
\details{
Please refer to the function "qvalue" in the package "qvalue" for details.

<URL:http://www.bioconductor.org/packages/release/bioc/html/qvalue.html>
}
\value{
the q-values for controling the porportion of host sequences in the predicted viral sequences.
}
\references{
Storey JD. (2002) A direct approach to false discovery rates.
Journal of the Royal Statistical Society, Series B, 64: 479-498.
<URL:
http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00346/abstract>
Storey JD and Tibshirani R. (2003) Statistical significance for
genome-wide experiments. Proceedings of the National Academy of
Sciences, 100: 9440-9445.
<URL: http://www.pnas.org/content/100/16/9440.full>
}

\examples{

## predict for contigs
inFaFile <- system.file("data", "contigs.fa", package="VirFinder")
predResult <- VF.pred(inFaFile)

## estimate q-values based on p-values
predResult$qvalue <- VF.qvalue(predResult$pvalue)

## sort sequences by q-value in ascending order
predResult[order(predResult$qvalue),]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
