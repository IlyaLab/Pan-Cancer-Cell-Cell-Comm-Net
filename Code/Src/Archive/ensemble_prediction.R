
  # http://xgboost.readthedocs.io/en/latest/R-package/xgboostPresentation.html
  require(xgboost)
  library(ggplot2)
  library(stringr)

# load the c2cWts
  load("dat_for_prediction.rda")
  scores <- dat

  c2cWtsT <- t(c2cWts)
  colnames(c2cWtsT) <- c2cNet$Path
  
  all(rownames(c2cWtsT) == resp$SampleBarcode)
  #[1] TRUE
  all(tcga$ShortBarcode == rownames(c2cWtsT))
  #[1] TRUE
  
  # correlate paths with clin features of interest.
  
  
  pairWise <- function(a,b) {
    
    if (length(unique(a)) <= 1) {
      return(1)
    }
    
    if (length(unique(b)) <= 1) {
      return(1)
    }
    
    if(is.numeric(a) & is.numeric(b)) {
      return(cor.test(a, b, method="spearman", use = 'pairwise.complete.obs')$p.value)
    }
    
    if(is.numeric(a) & is.character(b)) {
      bb <- as.factor(b)
      # anova
      res0 <- summary(aov(a~bb))
      return(as.numeric(res0[[1]][[5]])[1])
    }
    
  }
  
  # get the clinical data in
  #clinfeature <- tcga$mRNABaileyCluster
  clinfeature <- tcga$NumLymphNodes
  
  featCor <- sapply(1:ncol(c2cWtsT), function(i) pairWise(c2cWtsT[,i], clinfeature) )
  
  selFeat <- colnames(c2cWtsT)[featCor < 0.01]
  
  dat <- c2cWtsT[,selFeat]  
  
  superheat(X = t(dat), yt = clinfeature,
            dist.method = 'manhattan', scale = F, 
            pretty.order.rows = T, pretty.order.cols = T)
  
  save(dat, tcga, file='dat_for_prediction.rda')
  
  
# do some feature pre-selection #
  testPreds <- list()
  modelList <- list()
  featList <- list()
  fdx <- 0

for (idx in 1:nrow(scores)) {

    # leave one out
    train <- scores[-idx,]
    test <- scores[idx,]
    trainlabel <- ifelse(clinfeature[-idx]==fdx, yes=1, no=0) 
    testlabel  <- ifelse(clinfeature[idx] ==fdx, yes=1, no=0)  #clinfeature[idx]   # ifelse(clinfeature[idx]  == fdx, yes=0, no=1)

    trainlabel[is.na(trainlabel)] <- 0
    testlabel[is.na(testlabel)] <- 0
    
    # build an model for each
    #for (jdx in 1:nrow(train)) {
    #  trainJ <- train[-jdx,]
    #  trainlabelJ <- trainlabel[-jdx]
    modJ <- xgboost(data = as.matrix(train), label = trainlabel, max.depth = 3, eta = 0.3, nthread = 2, nround = 21, subsample=0.7, gamma=0.2)
      #modJ <- xgboost(data = as.matrix(trainJ), label = trainlabelJ, max.depth = 10, eta = 0.1, nthread = 20, nround = 35, subsample=0.7, gamma=0.1)
    modelList[[idx]] <- modJ
    featList[[idx]] <- xgb.importance(model=modJ)$Feature[1:10]
    #fdx <- fdx+1
    #}

    # then we predict on this one hold out.
    testPreds[[idx]] <- predict(modJ, matrix(test,nrow=1))
    #preds <- lapply(modelList, function(m) predict(m, matrix(test,nrow=1)))
    #testPreds[[idx]] <- unlist(preds)
}


#meanPreds <- unlist(lapply(testPreds, function(a) mean(a)))
res0 <- data.frame(Raw=unlist(testPreds), Rnd=round(unlist(testPreds)), Ys=ifelse(clinfeature == fdx, yes=0, no=1), Clus=clinfeature)

library(ROCR)
pred <- prediction(res0$Raw, res0$Ys)
roc.perf = performance(pred, measure = "prec", x.measure = "rec")
plot(roc.perf)
abline(a=0, b= 1)

table(res0$Rnd, Y=res0$Ys)


top_features <- sort(table( unlist(featList) ), decreasing=T)


