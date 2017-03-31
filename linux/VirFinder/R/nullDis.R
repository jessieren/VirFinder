nullDis <- function(trainDataHost, trainDataVirus, trainModOut)
{
  set.seed(1)
  trainIDHost <- sample(nrow(trainDataHost), nrow(trainDataHost)*3/4, replace=FALSE)
  trainIDVirus <- sample(nrow(trainDataVirus), nrow(trainDataVirus)*3/4, replace=FALSE)
  train <- rBind(trainDataHost[trainIDHost, ], trainDataVirus[trainIDVirus, ])
  label <- c(rep(0, length(trainIDHost)), rep(1, length(trainIDVirus))) 
  
  ## no need for remove the least significant word. 
  ## it takes too long and for prediction it does not matter
  #print("...training the model...")
  lasso.mod <- glmnet(train,label,alpha=1, family="binomial", lambda=trainModOut$lambda.min)
  scores_host <- predict(lasso.mod, type="response", trainDataHost[-trainIDHost, ])
  #scores_virus <- predict(lasso.mod, type="response", trainDataVirus[-trainIDVirus, ])
  
  #pred <- prediction(c(scores_host, scores_virus), c(rep(0, length(scores_host)), rep(1, length(scores_virus))))
	#perf <- performance(pred,"tpr","fpr")
  scores_host
}

