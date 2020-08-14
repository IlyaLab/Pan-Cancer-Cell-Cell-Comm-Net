
  # http://xgboost.readthedocs.io/en/latest/R-package/xgboostPresentation.html
  require(xgboost)
  library(ggplot2)
  library(stringr)
  library(superheat)
  library(ROCR)

# load the c2cWts
  load("dat_for_prediction.rda")
  load("assets/predicting_data.rda")
  
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
  clinfeature <- resp$PFI_1
  
  featCor <- sapply(1:ncol(c2cWtsT), function(i) pairWise(c2cWtsT[,i], clinfeature) )
  
  selFeat <- colnames(c2cWtsT)[featCor < 0.01]
  
  dat <- c2cWtsT[,selFeat]  
  
  superheat(X = dat, yr = clinfeature,
            dist.method = 'manhattan', scale = F, 
            pretty.order.rows = T, pretty.order.cols = T)
  
  save(dat, tcga, file='dat_for_prediction_bailey_clusters.rda')
  

pred <- list()
perf <- list()
allres <- list()
modelList <- list()
featList <- list()
testPreds <- list()
fdx <- 1

for (idx in 1:nrow(dat)) {

    # leave one out
    train <- dat[-idx,]
    test <- dat[idx,]
    trainlabel <- ifelse(clinfeature[-idx] == fdx, yes=1, no=0) 
    testlabel  <- ifelse(clinfeature[idx]  == fdx, yes=1, no=0) 

    modJ <- xgboost(data = as.matrix(train), label = trainlabel, max.depth = 3, eta = 0.3, nthread = 2, nround = 51, subsample=0.7, gamma=0.2)
    modelList[[idx]] <- modJ
    featList[[idx]] <- xgb.importance(model=modJ)$Feature[1:10]

    # then we predict on this one hold out.
    testPreds[[idx]] <- predict(modJ, matrix(test,nrow=1))
  }

res0 <- data.frame(Raw=unlist(testPreds), Rnd=round(unlist(testPreds)), Ys=ifelse(clinfeature == fdx, yes=1, no=0), Clus=clinfeature)
  
pred1 <- prediction(res0$Raw, res0$Ys)
predrand <- prediction(res0$Raw, sample(x = c(0,1), size = length(res0$Raw), replace = T))

plot(performance(pred1, measure = "prec", x.measure = "rec"), ylim=c(0,1), col=2)
plot(performance(predrand, measure = "prec", x.measure = "rec"), add=T, col=1)


plot(performance(pred1, measure = "tpr", x.measure = "fpr"), ylim=c(0,1), col=2)
plot(performance(predrand, measure = "tpr", x.measure = "fpr"), add=T, col=1)


top_features <- sort(table( unlist(featList) ), decreasing=T)

featuredf <- data.frame(Edge=names(top_features), Count=top_features)
write.table(featuredf, file="pfi_top_features_v2.csv", sep=',', quote=F, row.names = F)


