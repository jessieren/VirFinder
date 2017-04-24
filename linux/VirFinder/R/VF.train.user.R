VF.train.user <- function(trainFaFileHost, trainFaFileVirus, userModDir, userModName, w, equalSize=FALSE)
{
  #w <- 8
  pairWords <- findUniquePairWords(w, 4)

  VF.trainModUser <- w
  subLengthAll <- c(0.5, 1, 3) * 10^3
  for(subLength in subLengthAll)
  {
    trainDataHost <- trainDataCollect(trainFaFileHost, subLength, w, userModDir)
    colnames(trainDataHost) <- pairWords

    trainDataVirus <- trainDataCollect(trainFaFileVirus, subLength, w, userModDir)
    colnames(trainDataVirus) <- pairWords
    
    if(equalSize == TRUE)
    {
      size <- min(nrow(trainDataHost), nrow(trainDataVirus))
      hostID <- sample(1:nrow(trainDataHost), size, replace=FALSE)
      virusID <- sample(1:nrow(trainDataVirus), size, replace=FALSE)
      trainDataHost <- trainDataHost[hostID, ]
      trainDataVirus <- trainDataVirus[virusID, ] 
    }
    
    trainModOut <- trainMod(trainDataHost, trainDataVirus)
    nullDisOut <- nullDis(trainDataHost, trainDataVirus, trainModOut)
    
    if(subLength == 0.5*10^3)
    {
      attr(VF.trainModUser, "lasso.mod_0.5k") <- trainModOut
      attr(VF.trainModUser, "nullDis_0.5k") <- nullDisOut
    }else if (subLength == 1*10^3){
      attr(VF.trainModUser, "lasso.mod_1k") <- trainModOut
      attr(VF.trainModUser, "nullDis_1k") <- nullDisOut
    }else if (subLength == 3*10^3){
      attr(VF.trainModUser, "lasso.mod_3k") <- trainModOut
      attr(VF.trainModUser, "nullDis_3k") <- nullDisOut
    }
  }
  
  save(VF.trainModUser, file=file.path(userModDir, paste("VF.trainModUser.", userModName, ".rda", sep="")))
  
  VF.trainModUser
}


