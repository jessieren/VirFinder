trainDataCollect <-
function(trainFaFile, subLength, w)
{
  print(paste("...splitting ", trainFaFile, " into ", subLength, " bp...", sep=""))
  seqTrainKmerCount <- NULL
  countSeq <- 0
  countSubSeq <- 0
  flag <- 0
  seqFa <- NULL
  
  con <- file(trainFaFile, open = "r")
  while ( length(oneLine <- readLines(con, n = 1)) > 0 ) 
  {
    # 0628 last line is empty
    if( is.na(oneLine) | oneLine == "" )
    {
      next
    }
    
    seqTrain <- strsplit(oneLine, "")[[1]]
    
    if( seqTrain[1] == ">" && flag == 0 )
    {
      #### start the first contig ####
      flag <- 1
      
    }else if ( seqTrain[1] == ">" && flag > 1 )
    {
      countSeq <- countSeq + 1
      if( countSeq %% 10 == 0)
      {
        print(paste(".....seq ", countSeq, ".....", sep=""))
      }
      seqSplitOut <- seqSplit(seqFa, subLength, w)

      if( !is.null(seqSplitOut) )
      {
        #print("in")
        seqSplitOutSparse <- Matrix(seqSplitOut, sparse=TRUE)
        if( is.null(seqTrainKmerCount) )
        {
          seqTrainKmerCount <- seqSplitOutSparse
        }else{
          seqTrainKmerCount <- rBind(seqTrainKmerCount, seqSplitOutSparse)
        }
        countSubSeq <- countSubSeq + nrow(seqSplitOut)
      }
      
      seqFa <- NULL
      flag <- 1
      
    }else{
    
      seqFa <- c(seqFa, seqTrain)
      #write(line, file=currentContigSeqFile, append=TRUE)
      #sink()
      flag <- flag + 1
    }
  } 
  close(con)
  
  countSeq <- countSeq + 1
  if( countSeq %% 10 == 0)
  {
    print(paste(".....seq ", countSeq, ".....", sep=""))
  }
  seqSplitOut <- seqSplit(seqFa, subLength, w)
  if( !is.null(seqSplitOut) )
  {
    #print("in")
    seqSplitOutSparse <- Matrix(seqSplitOut, sparse=TRUE)
    if( is.null(seqTrainKmerCount) )
    {
      seqTrainKmerCount <- seqSplitOutSparse
    }else{
      seqTrainKmerCount <- rBind(seqTrainKmerCount, seqSplitOutSparse)
    }
    countSubSeq <- countSubSeq + nrow(seqSplitOut)
  }

  #print(paste("trainSeqFile:", trainFaFile, sep=""))
  print(paste("#seq:", countSeq, ", #subseq:", countSubSeq, sep=""))
  #print(seqTrainKmerCount)
  
  seqTrainKmerCount
}
