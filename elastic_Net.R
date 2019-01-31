setwd('//ATKK/home/m/mukhan/Documents/assignment_drugTargets/')


library('glmnet')

randomIterations = c(842,1408,1203,2016)
for(r in 1:length(randomIterations))
{
  load("drugResponse.RData")
  Y = as.matrix(Y)
  numOfSamples = dim(Y)[1]
  numOfFeatures = dim(X)[2]
  #log transformation of the data so as to enhance the signal as most 
  # values close to zero can be spread out so that we can see the distribution
  #hist(X)
  X=log2(1+X)
  #hist(X)
  #cbind(rownames(X),rownames(Y))]]
  numOfFolds = 10 #cv performance show
  set.seed(randomIterations[r])
  cvFolds = sample(numOfFolds,numOfSamples,replace = T)
  
  #alpha grid is established to test the values
  alphas = c(seq(0,1,by=0.01))[2:100]
  alphas = c(seq(0,1,by=0.1))[2:10]
  mse=matrix(NA,length(alphas),1)
  spCor=matrix(NA,length(alphas),1)
  prCor=matrix(NA,length(alphas),1)
  
  for(a in 1:length(alphas))
  {
    
    YPred = matrix(NA,nrow=nrow(Y),ncol=ncol(Y))
    rownames(YPred)=rownames(Y)
    for(fold in 1:numOfFolds){
      
      trainingIndex = which(fold!=cvFolds)
      testIndex = which(fold==cvFolds)
      
      Xtrain = X[trainingIndex,]
      Xtest = X[testIndex,]
      
      Ytrain = as.matrix(Y[trainingIndex,])
      Ytest = as.matrix(Y[testIndex,])
      
      #preprocessing Z-transformation
      meanXtrain = apply(Xtrain,2,mean) #feature-wise mean
      sdXtrain = apply(Xtrain,2,sd) #feature-wise sd
      
      Xtrain = (Xtrain -  matrix(meanXtrain, nrow=dim(Xtrain)[1],ncol=dim(Xtrain)[2],byrow = TRUE)/matrix(sdXtrain, nrow=dim(Xtrain)[1],ncol=dim(Xtrain)[2],byrow = TRUE))
      Xtest = (Xtest - matrix(meanXtrain, nrow=dim(Xtest)[1],ncol=dim(Xtest)[2],byrow = TRUE)/matrix(sdXtrain, nrow=dim(Xtest)[1],ncol=dim(Xtest)[2],byrow = TRUE))
      
      meanYtrain = apply(Ytrain,2,mean) #feature-wise mean
      sdYtrain = apply(Ytrain,2,sd) #deature-wise sd
      
      Ytrain = (Ytrain -  matrix(meanYtrain, nrow=dim(Ytrain)[1],ncol=dim(Ytrain)[2],byrow = TRUE))/(matrix(sdYtrain, nrow=dim(Ytrain)[1],ncol=dim(Ytrain)[2],byrow = TRUE))
      Ytest = (Ytest - matrix(meanYtrain, nrow=dim(Ytest)[1],ncol=dim(Ytest)[2],byrow = TRUE))/(matrix(sdYtrain, nrow=dim(Ytest)[1],ncol=dim(Ytest)[2],byrow = TRUE))
      
      set.seed(randomIterations[r])
      innercvFolds = sample(numOfFolds,dim(Ytrain)[1],replace = T)
      cvModel = cv.glmnet(x = Xtrain, y = Ytrain,
                          alpha = alphas[a],
                          standardize = FALSE,
                          intercept = FALSE,
                          family = 'gaussian',
                          type.measure = 'mse',
                          nfolds = 10,
                          foldid=innercvFolds)
      
      lambda.1se=cvModel$lambda.1se
      
      model=glmnet(x = Xtrain, y = Ytrain,
                   alpha = alphas[a],
                   lambda = lambda.1se,
                   standardize = FALSE,
                   intercept = FALSE,
                   family = 'gaussian')
      
      YPred.tmp = predict(model,Xtest)
      #reverting to orignal space
      YPred[testIndex,1] = (YPred.tmp * matrix(sdYtrain, nrow=dim(Ytest)[1],ncol=dim(Ytest)[2],byrow = TRUE)) + (matrix(meanYtrain, nrow=dim(Ytest)[1],ncol=dim(Ytest)[2],byrow = TRUE))
      
    }
    
    mse[a,1] = mean((Y-YPred)^2)
    spCor[a,1] = cor(Y,YPred,method = 'spearman')
    prCor[a,1] = cor(Y,YPred,method = 'pearson')
    
  }
  
  
  #index = #c(which.min(alphas),which.max(spCor),which.max(prCor))
  alphaChosen = alphas[which.max(spCor)]
  
  
  #Full data and feature selection
  meanX = apply(X,2,mean) #feature-wise mean
  sdX = apply(X,2,sd) #deature-wise sd
  X = (X -  matrix(meanX, nrow=dim(X)[1],ncol=dim(X)[2],byrow = TRUE)/matrix(sdX, nrow=dim(X)[1],ncol=dim(X)[2],byrow = TRUE))
  
  meanY = apply(Y,2,mean) #feature-wise mean
  sdY = apply(Y,2,sd) #deature-wise sd
  Y = (Y -  matrix(meanY, nrow=dim(Y)[1],ncol=dim(Y)[2],byrow = TRUE))/(matrix(sdY, nrow=dim(Y)[1],ncol=dim(Y)[2],byrow = TRUE))
  
  
  cvModelFullData = cv.glmnet(x = X, y = Y,
                              alpha = alphaChosen,
                              standardize = FALSE,
                              intercept = FALSE,
                              family = 'gaussian',
                              type.measure = 'mse',
                              nfolds = 10)
  lambda.1se=cvModelFullData$lambda.1se
  
  modelFullData = glmnet(x = X, y = Y,
                         alpha = alphaChosen,
                         lambda = lambda.1se,
                         standardize = FALSE,
                         intercept = FALSE,
                         family = 'gaussian')
  importantGenes = modelFullData$beta[which(modelFullData$beta!=0),1]
  #index=which(coef(modelFullData,s="lambda.min")!=0)
  #coef(modelFullData,s="lambda.min")[index,1]
  #modelFullData$beta[coefficients,1]
  featureSelectionResults = list()
  featureSelectionResults$cvMSE=mse
  featureSelectionResults$cvSpCor=spCor
  featureSelectionResults$cvprCor=prCor
  featureSelectionResults$cvAlpha=alphaChosen
  featureSelectionResults$cvLambda.1se=lambda.1se
  featureSelectionResults$modelFullData=modelFullData
  featureSelectionResults$selectedGenes=importantGenes
  
  
  # if you want to see the plots and list of genes uncomment the following lines
  #windows()
  #plot(alphas,featureSelectionResults$cvSpCor, ylab = "spearman correlation", xlab = "alpha values", main = "Optimal alpha selection")
  #pnt <- identify(alphas, featureSelectionResults$cvSpCor, plot = T)
  #This colors those points red
  #points(alphas[pnt], featureSelectionResults$cvSpCor[pnt], col = "red",  pch=(19), cex=(1.5))
  
  #print(names(featureSelectionResults$selectedGenes))
  
  
  save(featureSelectionResults,file=paste0('FeatureSelection_Results_Iteration_',randomIterations[r],'.RData'))
  print(paste0(r,' Iteration ',randomIterations[r]))
}