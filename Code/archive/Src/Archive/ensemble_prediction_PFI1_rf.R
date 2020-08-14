
  library(randomForest)
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
  

pred <- list()
perf <- list()
allres <- list()
modelList <- list()
accList <- list()
giniList <- list()
testPreds <- list()
fdx <- 1

for (idx in 1:nrow(dat)) {
    print(idx)
    # leave one out
    train <- dat[-idx,]
    test <- dat[idx,]
    trainlabel <- ifelse(clinfeature[-idx] == fdx, yes=1, no=0) 
    testlabel  <- ifelse(clinfeature[idx]  == fdx, yes=1, no=0) 

    modJ <- randomForest(x = as.matrix(train), y = as.factor(trainlabel), importance = T, ntree=2000)
    modelList[[idx]] <- modJ
    
    imptable <-  as.data.frame(modJ$importance)  # xgb.importance(model=modJ)$Feature[1:10]
    
    accList[[idx]] <- head(imptable[order(imptable$MeanDecreaseAccuracy, decreasing = T),])
    giniList[[idx]] <- head(imptable[order(imptable$MeanDecreaseGini, decreasing = T),])
    
    testPreds[[idx]] <- predict(modJ, matrix(test,nrow=1), type = "prob")[2]
  }

res0 <- data.frame(Raw=unlist(testPreds), Ys=ifelse(clinfeature == fdx, yes=1, no=0))
  
pred1    <- prediction(res0$Raw, res0$Ys)
predrand <- prediction(res0$Raw, sample(x = c(0,1), size = length(res0$Raw), replace = T))

plot(performance(pred1, measure = "prec", x.measure = "rec"), ylim=c(0,1), col=2)
plot(performance(predrand, measure = "prec", x.measure = "rec"), add=T, col=1)

plot(performance(pred1, measure = "tpr", x.measure = "fpr"), ylim=c(0,1), col=2)
abline(a=0, b=1)

acc_features  <- sort(table(unlist(lapply(accList, rownames))), decreasing=T)
gini_features <- sort(table(unlist(lapply(giniList, rownames))), decreasing=T)

featuredf <- data.frame(Acc=names(acc_features), Count=acc_features)
write.table(featuredf, file="pfi_top_features_acc.csv", sep=',', quote=F, row.names = F)

featuredf <- data.frame(Gini=names(gini_features), Count=gini_features)
write.table(featuredf, file="pfi_top_features_gini.csv", sep=',', quote=F, row.names = F)


x <- data.frame(PFI=resp$PFI_1,  EdgeWts=c2cWtsT[,'ductal_LAMC2_ITGA2_Plasma.cells'])
qplot(data=x, y=EdgeWts, x=as.factor(PFI), geom="boxplot") + theme_bw() + xlab("PFI")

