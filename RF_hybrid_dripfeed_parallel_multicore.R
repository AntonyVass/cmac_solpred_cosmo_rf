library(randomForest)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

# ----- set random seed
set.seed(1989)


# ----- I/O
input_data_path <- "path/to/datadir/dataset_minus_confidential.txt"
workspace_img_path <- "path/to/outputdir/results_rf_hybrid_dripfeed.RData"

# ----- 
masterSet <- read.csv(input_data_path, sep = "\t")

# ----- get number of solutes in dataset matching your condition (in this case, min. 5 instances)
count_Solute <- as.data.frame(table(masterSet$solute_name))
count_Solute <- count_Solute[which(count_Solute$Freq > 4),]
nsolu <- length(count_Solute[,1])

# ----- set fixed max number of solute data points to dripfeed into training set, even if more are possible
# ----- (should be at least 1 less than total so there is something left to test)
dripfeed_limit <- 4

### PARALLEL SETUP
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

myfunction <- function(iter_c_, iter_d_, all_combinations, solute_subset, init_trainingSet, pred_result_store, pred_error_store){
  # packages needed on workers
  library(randomForest)
  library(dplyr)
  
  # Split train/test
  print(paste("    iter_c_ = ", iter_c_))
  if (iter_d_ == 0) { 
    solutes_to_return <- all_combinations[,iter_c_]
    testSet <- solute_subset[solutes_to_return,]
    rest_solute <- solute_subset[-solutes_to_return,]
    trainingSet <- rbind(init_trainingSet, rest_solute)
  }  else {
    solutes_to_return <- all_combinations[,iter_c_]
    testSet <- solute_subset[-solutes_to_return,]
    rest_solute <- solute_subset[solutes_to_return,]
    trainingSet <- rbind(init_trainingSet, rest_solute)
  }
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
  data.pred <- predict(data.rf, testSet)
  te_exc_cols$PREDICTION <- data.pred
  
  if (iter_d_ == 0) {
    te_exc_cols$Training_Solvents <- rep("none", length(te_exc_cols$solute_name))
  } else {
    te_exc_cols$Training_Solvents <- rep(paste(rest_solute$solvent_name, sep = '', collapse = ' '), length(te_exc_cols$solute_name))
  }
  
  return(te_exc_cols)
}

# ----- Results storage
top_result_store <- vector("list",nsolu)
top_exectime_store <- vector("list",nsolu)

# ----- loop 1: iterate over number of CANDIDATE SOLUTES, "iter_s_"

overall_start_time <- Sys.time()
for (iter_s_ in 1:nsolu) {
  cat("----------\n","solute number", iter_s_, "of", nsolu, "\n----------\n")
  
  # ----- specify test solute and determine how many data points can be dripfed into TRAINING SET, given constraints
  testSolute <- count_Solute$Var1[iter_s_]
  fsolu <- count_Solute$Freq[iter_s_]
  if (dripfeed_limit > count_Solute$Freq[iter_s_]) {
    dripfeed_max <- fsolu-1
  } else {
    dripfeed_max <- dripfeed_limit
  }
  
  # ----- initial prep
  solute_subset <- masterSet[masterSet$solute_name==testSolute,]
  init_trainingSet <- anti_join(masterSet, solute_subset, by=colnames(masterSet))
  per_dripfeed_size_result_store <- vector("list",dripfeed_max+1)
  per_dripfeed_size_exectime_store <- vector("list",dripfeed_max+1)
  plotData <- data.frame()
  results <- data.frame ()

  # ----- loop 2: iterate over number of SOLUTE DATA POINTS to place in TRAINING SET,  "iter_d_"
  # ----- NOTE: all remaining data points are placed in test set

  for(iter_d_ in 0:dripfeed_max) {
    start_time <- Sys.time()
    cat("     dripfeeding", iter_d_, "solute data point(s) into training\n")
    
    # ----- calculate number of combinations possible
    if (iter_d_ == 0) {
      all_combinations <- cbind(1:fsolu)
    } else {
      all_combinations <- combn(fsolu,iter_d_)
    }
    ncombos <- length(all_combinations[1,])
 
    # ----- loop 3A: iterate over number of COMBINATIONS FOR THIS NUMBER OF DATA POINTS  "iter_c_"
    # ----- NOTE: contains the ML
    
    cat("     no. of combinations to calculate: ",ncombos,"\n")
    
    plotData <- foreach(iter_c_ = 1:ncombos, .combine='rbind') %dopar% myfunction(iter_c_, iter_d_, all_combinations, solute_subset, init_trainingSet, pred_result_store, pred_error_store)
    
    end_time <- Sys.time()
    exec_time <- end_time - start_time
    print(exec_time)
    
    per_dripfeed_size_result_store[[iter_d_+1]] <- plotData
    per_dripfeed_size_exectime_store[[iter_d_+1]] <- list(exec_time, exec_time[[1]])

    } # loop 2 close
  
  top_result_store[[iter_s_]] <- per_dripfeed_size_result_store
  top_exectime_store[[iter_s_]] <- per_dripfeed_size_exectime_store
  
  } # loop 1 close

overall_end_time <- Sys.time()
overall_exec_time <- overall_end_time - overall_start_time
print(overall_exec_time)

# ----- Write to file (crucial: nothing is saved until you do this)
save.image(workspace_img_path)

# ----- Release workers
stopCluster(cl)
