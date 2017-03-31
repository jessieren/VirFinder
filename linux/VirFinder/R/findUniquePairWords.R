findUniquePairWords <-
function(WWW, ZI)
{
#uniquePairIDs <- NULL
uniquePairWords <- NULL
for(ID in 1:ZI^WWW)
{
word1234 <- ten2fourCplusOrder(ID, WWW, ZI)
wordACGT <- convert1234ToACGT(word1234)
complimentWord <- reverse(word1234, ZI)
complimentWordID <- four2tenCplusOrder(complimentWord, ZI)
if(complimentWordID >= ID)
{
#uniquePairIDs <- c(uniquePairIDs, ID)
uniquePairWords <- c(uniquePairWords, paste(wordACGT, collapse=""))
}
}
names(uniquePairWords) <- seq(1,length(uniquePairWords))
uniquePairWords
}
