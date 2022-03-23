library(randomForest)
library(dplyr)

# ----- set random seed
set.seed(808)

# ----- files
input_data_path <- "path/to/datadir/dataset_minus_confidential.txt"
results_output_path <- "path/to/outputdir/filename.txt"
workspace_img_path <- "path/to/outputdir/filename.RData"

# ----- Load and Randomly shuffle data
masterSet <- read.csv(input_data_path, sep = "\t")
shuffled_masterSet <- masterSet[sample(nrow(masterSet)),] # skip if holdout set used

# ----- set k folds
k <- 10

# ----- Create k equally size folds
folds <- cut(seq(1,nrow(shuffled_masterSet)),breaks=k,labels=FALSE)

# ----- pre-allocate vectors for storing results
rf_store <- vector("list",k)
pred_result_store <- vector("list",k)
pred_mae_store <- vector("list",k)
pred_rmse_store <- vector("list",k)
rebuild_set <- data.frame()

# ----- Start timer
start_time <- proc.time()

# ----- Perform k-fold cross validation
for(i in 1:k){
  print(i)
  
  # Split train/test
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testSet <- shuffled_masterSet[testIndexes, ]
  trainingSet <- shuffled_masterSet[-testIndexes, ]
  trainingResponse<- trainingSet$solubility_g_100g_log
  testResponse <- testSet$solubility_g_100g_log
  
  # Remove cols not for ML (metadata, response, etc)
  cols_to_exc <- c(1:13,15:20) # this leaves COSMO-RSD predictions in
  tr_exc_cols <- trainingSet[,1:20]
  te_exc_cols <- testSet[,1:20]
  trainingSet <- select(trainingSet, -cols_to_exc)
  testSet <- select(testSet, -cols_to_exc)
  
  # RF
  data.rf <- randomForest(trainingSet,
                          trainingResponse,
                          ntree=1000,
                          type=regression,
                          replace=TRUE,
                          keep.forest=TRUE,
                          importance=FALSE,
                          proximity=FALSE
                          )
  
  # Store raw results
  rf_store[[i]] <- data.rf
  data.pred <- predict(data.rf, testSet)
  results <- data.frame(te_exc_cols, data.pred)
  pred_result_store[[i]] <- results
  pred_mae_store[i] <- sum(abs((results$data.pred-testResponse)/length(data.pred)))
  pred_rmse_store[i] <- sqrt(mean((data.pred - testResponse)^2))
  
  # Append fold predictions to final result dataframe
  if (i == 1) {rebuild_set <- results} else {rebuild_set <- rbind(rebuild_set, results)}
}

# ----- End timer
elapsed_time <- proc.time() - start_time

# ----- Performance metrics

r2_general <-function(preds,actual){ 
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

r2_COSMO <- r2_general(rebuild_set$cosmo_solubility_g_100g_log, rebuild_set$solubility_g_100g_log)
rmse_COSMO <- sqrt(mean((rebuild_set$cosmo_solubility_g_100g_log - rebuild_set$solubility_g_100g_log)^2))
mae_COSMO <- mean(abs(rebuild_set$cosmo_solubility_g_100g_log - rebuild_set$solubility_g_100g_log))
r2_RF <- r2_general(rebuild_set$data.pred, rebuild_set$solubility_g_100g_log)
rmse_RF <- sqrt(mean((rebuild_set$data.pred - rebuild_set$solubility_g_100g_log)^2))
mae_RF <- mean(abs(rebuild_set$data.pred - rebuild_set$solubility_g_100g_log))
count_imp <- length(rebuild_set$data.pred[abs(rebuild_set$data.pred - rebuild_set$solubility_g_100g_log) < abs(rebuild_set$cosmo_solubility_g_100g_log - rebuild_set$solubility_g_100g_log)])
frac_imp <- count_imp/nrow(rebuild_set)

writeLines(
  c(
    paste("R2 (COSMO) =", round(r2_COSMO, 3)),
    paste("RMSE (COSMO) =", round(rmse_COSMO, 3)),
    paste("MAE (COSMO) =", round(mae_COSMO, 3)),
    paste("R2 (RF) =", round(r2_RF, 3)),
    paste("RMSE (RF) =", round(rmse_RF, 3)),
    paste("MAE (RF) =", round(mae_RF, 3)),
    paste("Fraction Improved =", round(frac_imp, 3))
  )
)

# ----- Write to file
write.table(rebuild_set, results_output_path, row.names = FALSE, quote = FALSE, sep = "\t" )
save.image(workspace_img_path)

