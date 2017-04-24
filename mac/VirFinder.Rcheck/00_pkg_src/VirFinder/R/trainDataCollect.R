trainDataCollect <-
function(trainFaFile, subLength, w, userModDir)
{
  print(paste("...splitting ", trainFaFile, " into ", subLength, " bp...", sep=""))
  
  seqTrainKmerCountDir <- file.path(userModDir, "seqTrainKmerCount")
  dir.create(seqTrainKmerCountDir, recursive = TRUE, showWarnings = FALSE)

  seqTrainKmerCount <- NULL
  countSeq <- 0
  countSubSeq <- 0
  flag <- 0
  seqFa <- NULL
  fileNum <- 0
  
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
      if( countSeq %% 100 == 0)
      {
        print(paste(".....processing '>' seq ", countSeq, " in ", basename(trainFaFile), ".....", sep=""))
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
        #if( countSubSeq %% 100 == 0)
        #{
        #  print(paste(".....obstaining ", subLength/1000, "kb subSeq ", countSubSeq, ".....", sep=""))
        #}
      }
      
      ## split seqTrainKmerCount into several small documents to speed up
      if( !is.null(seqTrainKmerCount) && nrow(seqTrainKmerCount) >= 2000 )
      {
        fileNum <- fileNum + 1
        print(paste(".....collecting ", subLength/1000, "kb subSeq ", countSubSeq, ", saving kmerCount file", fileNum, "...", sep=""))
        save(seqTrainKmerCount, file=file.path(seqTrainKmerCountDir, paste("VF.trainKmer.", basename(trainFaFile), ".subLen", subLength, ".k", w, ".file", fileNum, ".RData", sep="")))
        seqTrainKmerCount <- NULL
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
  if( countSeq %% 100 == 0)
  {
    print(paste(".....processing '>' seq ", countSeq, ".....", sep=""))
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
    
    #if( countSubSeq %% 1000 == 0)
    #{
    #  print(paste(".....obtaining ", subLength/1000, "kb subSeq ", countSubSeq, " in ", basename(trainFaFile), ".....", sep=""))
    #}
  }
  
  if( !is.null(seqTrainKmerCount) )
  {
    fileNum <- fileNum + 1
    print(paste(".....collecting ", subLength/1000, "kb subSeq ", countSubSeq, ", saving kmerCount file", fileNum, "...", sep=""))
    save(seqTrainKmerCount, file=file.path(seqTrainKmerCountDir, paste("VF.trainKmer.", basename(trainFaFile), ".subLen", subLength, ".k", w, ".file", fileNum, ".RData", sep="")))
  }


  #print(paste("trainSeqFile:", trainFaFile, sep=""))
  print(paste("...done! #seq:", countSeq, ", #subseq:", countSubSeq, sep=""))
  #print(seqTrainKmerCount)
  
  ## reload seqTrainKmerCount
  seqTrainKmerCount <- NULL
  seqTrainKmerCountAll <- NULL
  for(num in 1:fileNum)
  {
    load(file.path(seqTrainKmerCountDir, paste("VF.trainKmer.", basename(trainFaFile), ".subLen", subLength, ".k", w, ".file", num, ".RData", sep="")))
    #print(dim(seqTrainKmerCount))
    if( is.null(seqTrainKmerCountAll) )
    {
      seqTrainKmerCountAll <- seqTrainKmerCount
    }else{
      seqTrainKmerCountAll <- rBind(seqTrainKmerCountAll, seqTrainKmerCount)
    }
  }
  
  seqTrainKmerCountAll
}
