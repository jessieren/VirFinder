convert1234ToACGT <-
function(word1234)
{
wordACGT <- word1234
wordACGT[word1234==1] <- "A"
wordACGT[word1234==2] <- "C"
wordACGT[word1234==3] <- "G"
wordACGT[word1234==4] <- "T"
wordACGT

}
