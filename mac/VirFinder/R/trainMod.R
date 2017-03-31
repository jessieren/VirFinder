trainMod <-
function(trainDataHost, trainDataVirus)
{
  train <- rBind(trainDataHost, trainDataVirus)
  label <- c(rep(0, nrow(trainDataHost)), rep(1, nrow(trainDataVirus))) 
  
  set.seed(1)
  grid=10^seq(2,-10,length=100)
  print("...training the model...")
  cv.mod <- cv.glmnet(train,label,alpha=1,family="binomial",type.measure="auc", lambda=grid, nfolds=5)
  cv.mod
}
