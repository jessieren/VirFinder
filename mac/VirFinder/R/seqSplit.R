seqSplit <-
function(seqTrain, subLength, w)
{
  featureOut_kmerCount <- NULL
  seqPos0 <- 1
  seqPos1 <- seqPos0 + subLength - 1
  while( seqPos1 <= length(seqTrain) )
  {
    #print(seqPos0)
    seqSub <- seqTrain[seqPos0:seqPos1]
    if( mean(seqSub == "N") < 0.3 )
    {
      featureOut <- countSeqFeatureCpp(seqSub, w)
      featureOut_kmerCount <- rbind(featureOut_kmerCount, featureOut$kmerCount)
    }
    
    seqPos0 <- seqPos1 + 1
    seqPos1 <- seqPos0 + subLength - 1
  }
  featureOut_kmerCount
}
