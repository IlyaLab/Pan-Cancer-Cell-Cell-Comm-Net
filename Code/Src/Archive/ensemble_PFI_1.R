

  # http://xgboost.readthedocs.io/en/latest/R-package/xgboostPresentation.html
  require(xgboost)
  library(ggplot2)
  library(stringr)

# load the c2cWts
  scores <- read.table("paad_gene_set_scores.tsv",sep='\t',stringsAsFactors=F)
  scores <- t(scores)
  rownames(scores) <- str_replace_all(rownames(scores), '\\.', '-')
  rownames(scores) <- str_sub(rownames(scores), 1, 12)

# get the clinical data in
  clin <- read.table("Survival_pancan_results-20171213-095944.csv", sep=',', stringsAsFactors=F, header=T)
  rownames(clin) <- clin$ParticipantBarcode
  clin <- clin[rownames(scores),]
  all(rownames(clin) == rownames(scores))
  #[1] TRUE

# choose the feature of interest (Bailey clusters)
  clinfeature <- clin$PFI_1


# some feature selection maybe    
  x <- sapply(1:ncol(scores), function(i) (t.test(scores[,i] ~ clinfeature)$statistic))
  #x <- cor(scores, clinfeature, method="spearman", use="pairwise.complete.obs")

  scores2 <- scores[,abs(x) > 1]

  clinfeature2 <- clinfeature

  # do some feature pre-selection #

  testPreds <- list()
  modelList <- list()
  featList <- list()
  fdx <- 1

for (idx in 1:nrow(scores2)) {

    # leave one out
    train <- scores2[-idx,]
    test <- scores2[idx,]
    trainlabel <- clinfeature2[-idx]
    testlabel  <- clinfeature2[idx]

    # build an model for each
    for (jdx in 1:nrow(train)) {
      trainJ <- train[-jdx,]
      trainlabelJ <- trainlabel[-jdx]
      #modJ <- xgboost(data = as.matrix(trainJ), label = trainlabelJ, max.depth = 3, eta = 0.3, nthread = 2, nround = 22, subsample=0.7, gamma=0.2)
      modJ <- xgboost(data = as.matrix(trainJ), label = trainlabelJ, max.depth = 10, eta = 0.1, nthread = 20, nround = 35, subsample=0.7, gamma=0.1)
      modelList[[jdx]] <- modJ
      featList[[fdx]] <- xgb.importance(model=modJ)$Feature[1:10]
      fdx <- fdx+1
    }

    # then we predict on this one hold out.
    preds <- lapply(modelList, function(m) predict(m, matrix(test,nrow=1)))
    testPreds[[idx]] <- unlist(preds)
}



meanPreds <- unlist(lapply(testPreds, function(a) mean(a)))


library(ROCR)
pred <- prediction(meanPreds, clinfeature)
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc.perf)
abline(a=0, b= 1)

top_features <- sort(table( unlist(featList) ), decreasing=T)


testPreds[[1]]
