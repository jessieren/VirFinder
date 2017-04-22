VF.pred <-
function(inFaFile)
{
  data(VF.trainMod8mer)
  w <- VF.trainMod8mer
  
  predResult <- NULL
  flag <- 0
  seqLength <- 0
  seqFa <- NULL
  contigCount <- 0
  lineNum <- 0
  line0 <- 0
  
  con <- file(inFaFile, open = "r")
  while ( length(line <- readLines(con, n = 1)) > 0 ) 
  {
    lineNum <- lineNum + 1
    
    # 0628 last line is empty
    if( is.na(line) | line == "" )
    {
      next
    }
    
    if( strsplit(line, "")[[1]][1] == ">" && flag == 0 )
    {
      #### start the first contig ####
      flag <- 1
      currentFileName <- strsplit(line, ">")[[1]][2]
      line0 <- lineNum
      
    }else if ( strsplit(line, "")[[1]][1] == ">" && flag > 1 )
    {
      ## count sequence feature
      featureOut <- countSeqFeatureCpp(seqFa, w)
      featureOut_kmerCount <- featureOut$kmerCount
      #names(featureOut_kmerCount) <- attr(VF.trainMod8mer, "uniquePairWords")
      
      ## predict
      ## seqlength ##
      seqLength <- length(seqFa)
      
      #print("predict using lasso")
      if( seqLength < 1*1000 )
      {
        lasso.mod <- attr(VF.trainMod8mer, "lasso.mod_0.5k")
        rmWordID <- attr(VF.trainMod8mer, "rmWordID_0.5k")
        nullDis <- attr(VF.trainMod8mer, "nullDis_0.5k")
      }else if( seqLength < 3*1000 ){
        lasso.mod <- attr(VF.trainMod8mer, "lasso.mod_1k")
        rmWordID <- attr(VF.trainMod8mer, "rmWordID_1k")
        nullDis <- attr(VF.trainMod8mer, "nullDis_1k")
      }else {
        lasso.mod <- attr(VF.trainMod8mer, "lasso.mod_3k")
        rmWordID <- attr(VF.trainMod8mer, "rmWordID_3k")
        nullDis <- attr(VF.trainMod8mer, "nullDis_3k")
      }
      
      lasso.pred <- predict(lasso.mod, t(as.matrix(featureOut_kmerCount[-rmWordID])), type="response")
      pvalue <- mean(nullDis > as.numeric(lasso.pred) )
      #write(paste(currentFileName, seqLength, lasso.pred, pvalue, sep=","), file=attr(objFa, "predFile"), append=TRUE)
      predResult <- rbind(predResult, c(currentFileName, seqLength, lasso.pred, pvalue))
      contigCount <- contigCount + 1
      print(paste(basename(inFaFile), paste("[", contigCount, "]", sep=""), paste("line", line0, "-", lineNum, sep=""), currentFileName, "len", seqLength, "score", round(lasso.pred, 4), "pvalue", round(pvalue, 4)))
      
      #### start the current contig ####
      seqFa <- NULL
      flag <- 1
      currentFileName <- strsplit(line, ">")[[1]][2]
      line0 <- lineNum
      
    }else{
      #### continue writing contig line ####
      seqFa <- c(seqFa, strsplit(line, "")[[1]])
      #write(line, file=currentContigSeqFile, append=TRUE)
      #sink()
      flag <- flag + 1
    }
  }
  close(con)
  
  ## count sequence feature
  featureOut <- countSeqFeatureCpp(seqFa, w)
  featureOut_kmerCount <- featureOut$kmerCount
  #names(featureOut_kmerCount) <- attr(VF.trainMod8mer, "uniquePairWords")
  
  ## predict
  ## seqlength ##
  seqLength <- length(seqFa)
  
  #print("predict using lasso")
  if( seqLength < 1*1000 )
  {
    lasso.mod <- attr(VF.trainMod8mer, "lasso.mod_0.5k")
    rmWordID <- attr(VF.trainMod8mer, "rmWordID_0.5k")
    nullDis <- attr(VF.trainMod8mer, "nullDis_0.5k")
  }else if( seqLength < 3*1000 ){
    lasso.mod <- attr(VF.trainMod8mer, "lasso.mod_1k")
    rmWordID <- attr(VF.trainMod8mer, "rmWordID_1k")
    nullDis <- attr(VF.trainMod8mer, "nullDis_1k")
  }else {
    lasso.mod <- attr(VF.trainMod8mer, "lasso.mod_3k")
    rmWordID <- attr(VF.trainMod8mer, "rmWordID_3k")
    nullDis <- attr(VF.trainMod8mer, "nullDis_3k")
  }
  
  lasso.pred <- predict(lasso.mod, t(as.matrix(featureOut_kmerCount[-rmWordID])), type="response")
  pvalue <- mean(nullDis > as.numeric(lasso.pred) )
  #write(paste(currentFileName, seqLength, lasso.pred, pvalue, sep=","), file=attr(objFa, "predFile"), append=TRUE)
  predResult <- rbind(predResult, c(currentFileName, seqLength, lasso.pred, pvalue))
  
  contigCount <- contigCount + 1
  print(paste(basename(inFaFile), paste("[", contigCount, "]", sep=""), paste("line", line0, "-", lineNum, sep=""), currentFileName, "len", seqLength, "score", round(lasso.pred, 4), "pvalue", round(pvalue, 4)))
  
  colnames(predResult) <- c("name", "length", "score", "pvalue")
  predResult_df <- as.data.frame(predResult)
  predResult_df$length <- as.numeric(as.character(predResult_df$length))
  predResult_df$score <- as.numeric(as.character(predResult_df$score))
  predResult_df$pvalue <- as.numeric(as.character(predResult_df$pvalue))
  predResult_df
  
}
