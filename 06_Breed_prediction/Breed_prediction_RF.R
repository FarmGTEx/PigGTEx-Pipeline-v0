#!/usr/bin/env Rscript
# library module ----------------------------------------------------------

library("data.table")
library("stringr")
library("plyr")
library("dplyr")
library("randomForest")
library("e1071")
library("caret")
library("reshape")

# Model training-------------------------------------------------------------------------

# Machine learning with a training set of 2237 individuals of purebred non-DLY
folds_2237_pure_oth <- createFolds(y=traindata_2237_pure_oth$FID,k = 3,list = T, returnTrain = FALSE)
folds_2237_pure_DLY <- createFolds(y=traindata_2237_pure_DLY$FID,k = 3,list = T, returnTrain = FALSE)
folds_2237_cross <- createFolds(y=traindata_2237_cross$FID,k = 3,list = T, returnTrain = FALSE)
for (i in 1:10) {
  for (j in 1:3) {
    # Machine learning with a training set of 2237 individuals of purebred non-DLY
    assign(paste("TestIndex_2237_pure_oth",i,j,sep = "_"),melt(data = folds_2237_pure_oth[j])$value)
    assign(paste("CV_test_2237_pure_oth",i,j,sep = "_"),traindata_2237_pure_oth[get(paste("TestIndex_2237_pure_oth",i,j,sep = "_")),])
    assign(paste("CV_train_2237_pure_oth",i,j,sep = "_"),traindata_2237_pure_oth[-get(paste("TestIndex_2237_pure_oth",i,j,sep = "_")),])
    assign(paste("Model_RF_pure_oth",i,j,sep = "_"),
           randomForest(x = get(paste("CV_train_2237_pure_oth",i,j,sep = "_"))[-c(1:6)],y = get(paste("CV_train_2237_pure_oth",i,j,sep = "_"))$FID))
    
    # Machine learning with a training set of 2237 individuals of purebred DLY
    assign(paste("TestIndex_pure_DLY",i,j,sep = "_"),melt(data = folds_2237_pure_DLY[j])$value)
    assign(paste("CV_test_pure_DLY",i,j,sep = "_"),traindata_2237_pure_DLY[get(paste("TestIndex_pure_DLY",i,j,sep = "_")),])
    assign(paste("CV_train_pure_DLY",i,j,sep = "_"),traindata_2237_pure_DLY[-get(paste("TestIndex_pure_DLY",i,j,sep = "_")),])
    assign(paste("Model_RF_pure_DLY",i,j,sep = "_"),
           randomForest(x = get(paste("CV_train_pure_DLY",i,j,sep = "_"))[-c(1:6)],y = get(paste("CV_train_pure_DLY",i,j,sep = "_"))$FID))
    
    # Machine learning with a training set of 2237 individuals of purebred DLY
    assign(paste("TestIndex_pure_DLY",i,j,sep = "_"),melt(data = folds_2237_pure_DLY[j])$value)
    assign(paste("CV_test_pure_DLY",i,j,sep = "_"),traindata_2237_pure_DLY[get(paste("TestIndex_pure_DLY",i,j,sep = "_")),])
    assign(paste("CV_train_pure_DLY",i,j,sep = "_"),traindata_2237_pure_DLY[-get(paste("TestIndex_pure_DLY",i,j,sep = "_")),])
    assign(paste("Model_RF_pure_DLY",i,j,sep = "_"),
           randomForest(x = get(paste("CV_train_pure_DLY",i,j,sep = "_"))[-c(1:6)],y = get(paste("CV_train_pure_DLY",i,j,sep = "_"))$FID))
  }
}


# Prediction results and secondary statistical analysis with a training set of 2,237 individuals of purebred non-DLY ------------------------------------------------

for (i in 1:10) {
  for (j in 1:3) {
    # Prediction of results
    # Prediction results and secondary statistical analysis with a training set of 2,237 individuals of purebred non-DLY
    assign(paste("prediction_pure_oth_class_7008",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_oth",i,j,sep = "_")),newdata = testdata_7008[-c(1:6)]))
    assign(paste("prediction_pure_oth_prob_7008",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_oth",i,j,sep = "_")),newdata = testdata_7008[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_pure_oth_prob_max_7008",i,j,sep = "_"),apply(X = get(paste("prediction_pure_oth_prob_7008",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_pure_oth_class_2237",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_oth",i,j,sep = "_")),newdata = testdata_2237[-c(1:6)]))
    assign(paste("prediction_pure_oth_prob_2237",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_oth",i,j,sep = "_")),newdata = testdata_2237[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_pure_oth_prob_max_2237",i,j,sep = "_"),apply(X = get(paste("prediction_pure_oth_prob_2237",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_pure_oth_class_4617",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_oth",i,j,sep = "_")),newdata = testdata_4617[-c(1:6)]))
    assign(paste("prediction_pure_oth_prob_4617",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_oth",i,j,sep = "_")),newdata = testdata_4617[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_pure_oth_prob_max_4617",i,j,sep = "_"),apply(X = get(paste("prediction_pure_oth_prob_4617",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_pure_oth_class_154",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_oth",i,j,sep = "_")),newdata = testdata_154[-c(1:6)]))
    assign(paste("prediction_pure_oth_prob_154",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_oth",i,j,sep = "_")),newdata = testdata_154[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_pure_oth_prob_max_154",i,j,sep = "_"),apply(X = get(paste("prediction_pure_oth_prob_154",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_pure_oth_class_pure_oth",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_oth",i,j,sep = "_")),newdata = testdata_2237_pure_oth[-c(1:6)]))
    assign(paste("prediction_pure_oth_prob_pure_oth",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_oth",i,j,sep = "_")),newdata = testdata_2237_pure_oth[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_pure_oth_prob_max_pure_oth",i,j,sep = "_"),apply(X = get(paste("prediction_pure_oth_prob_pure_oth",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    
    # Prediction results and secondary statistical analysis with a training set of 2,237 individuals of purebred DLY
    assign(paste("prediction_pure_DLY_class_7008",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_DLY",i,j,sep = "_")),newdata = testdata_7008[-c(1:6)]))
    assign(paste("prediction_pure_DLY_prob_7008",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_DLY",i,j,sep = "_")),newdata = testdata_7008[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_pure_DLY_prob_max_7008",i,j,sep = "_"),apply(X = get(paste("prediction_pure_DLY_prob_7008",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_pure_DLY_class_2237",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_DLY",i,j,sep = "_")),newdata = testdata_2237[-c(1:6)]))
    assign(paste("prediction_pure_DLY_prob_2237",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_DLY",i,j,sep = "_")),newdata = testdata_2237[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_pure_DLY_prob_max_2237",i,j,sep = "_"),apply(X = get(paste("prediction_pure_DLY_prob_2237",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_pure_DLY_class_4617",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_DLY",i,j,sep = "_")),newdata = testdata_4617[-c(1:6)]))
    assign(paste("prediction_pure_DLY_prob_4617",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_DLY",i,j,sep = "_")),newdata = testdata_4617[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_pure_DLY_prob_max_4617",i,j,sep = "_"),apply(X = get(paste("prediction_pure_DLY_prob_4617",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_pure_DLY_class_154",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_DLY",i,j,sep = "_")),newdata = testdata_154[-c(1:6)]))
    assign(paste("prediction_pure_DLY_prob_154",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_DLY",i,j,sep = "_")),newdata = testdata_154[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_pure_DLY_prob_max_154",i,j,sep = "_"),apply(X = get(paste("prediction_pure_DLY_prob_154",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_pure_DLY_class_pure_DLY",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_DLY",i,j,sep = "_")),newdata = testdata_2237_pure_DLY[-c(1:6)]))
    assign(paste("prediction_pure_DLY_prob_pure_DLY",i,j,sep = "_"),predict(object = get(paste("Model_RF_pure_DLY",i,j,sep = "_")),newdata = testdata_2237_pure_DLY[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_pure_DLY_prob_max_pure_DLY",i,j,sep = "_"),apply(X = get(paste("prediction_pure_DLY_prob_pure_DLY",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    
    # Prediction results and secondary statistical analysis with a training set of 2,237 individuals of cross breed
    assign(paste("prediction_cross_class_7008",i,j,sep = "_"),predict(object = get(paste("Model_RF_cross",i,j,sep = "_")),newdata = testdata_7008[-c(1:6)]))
    assign(paste("prediction_cross_prob_7008",i,j,sep = "_"),predict(object = get(paste("Model_RF_cross",i,j,sep = "_")),newdata = testdata_7008[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_cross_prob_max_7008",i,j,sep = "_"),apply(X = get(paste("prediction_cross_prob_7008",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_cross_class_2237",i,j,sep = "_"),predict(object = get(paste("Model_RF_cross",i,j,sep = "_")),newdata = testdata_2237[-c(1:6)]))
    assign(paste("prediction_cross_prob_2237",i,j,sep = "_"),predict(object = get(paste("Model_RF_cross",i,j,sep = "_")),newdata = testdata_2237[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_cross_prob_max_2237",i,j,sep = "_"),apply(X = get(paste("prediction_cross_prob_2237",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_cross_class_4617",i,j,sep = "_"),predict(object = get(paste("Model_RF_cross",i,j,sep = "_")),newdata = testdata_4617[-c(1:6)]))
    assign(paste("prediction_cross_prob_4617",i,j,sep = "_"),predict(object = get(paste("Model_RF_cross",i,j,sep = "_")),newdata = testdata_4617[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_cross_prob_max_4617",i,j,sep = "_"),apply(X = get(paste("prediction_cross_prob_4617",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_cross_class_154",i,j,sep = "_"),predict(object = get(paste("Model_RF_cross",i,j,sep = "_")),newdata = testdata_154[-c(1:6)]))
    assign(paste("prediction_cross_prob_154",i,j,sep = "_"),predict(object = get(paste("Model_RF_cross",i,j,sep = "_")),newdata = testdata_154[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_cross_prob_max_154",i,j,sep = "_"),apply(X = get(paste("prediction_cross_prob_154",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
    assign(paste("prediction_cross_class_cross",i,j,sep = "_"),predict(object = get(paste("Model_RF_cross",i,j,sep = "_")),newdata = testdata_2237_cross[-c(1:6)]))
    assign(paste("prediction_cross_prob_cross",i,j,sep = "_"),predict(object = get(paste("Model_RF_cross",i,j,sep = "_")),newdata = testdata_2237_cross[-c(1:6)],type = "prob") %>% as.data.frame(.))
    assign(paste("prediction_cross_prob_max_cross",i,j,sep = "_"),apply(X = get(paste("prediction_cross_prob_cross",i,j,sep = "_")), MARGIN = 1, function(x){max(x)}))
  }
}


#  -------------------------------------------------------------------






















