reverse <-
function (seqA, ZI)
{
n <- length(seqA)
seqB <- ZI+1 - seqA
seqC <- rep(0, n)
for (i in 1:n) {
seqC[i] <- seqB[n - i + 1]
}
seqC
}
