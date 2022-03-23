library(plotly)
library(viridis)
library(dplyr)


r2_general <-function(preds,actual){ 
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

RData_file <- "path/to/outputdir/results_rf_hybrid_dripfeed.RData"
load(RData_file)
masterSet <- read.csv("path/to/datadir/dataset_minus_confidential.txt", sep = "\t")

# --- flatten nested list of data
flat_list <- setNames(data.frame(matrix(ncol = 22, nrow = 0)), colnames(uber_result_store[[1]][[1]]))
for (h in 1:57) {
  print(h)
  for (i in 1:5) {
    flat_list <- rbind(flat_list, uber_result_store[[h]][[i]])
  }
}

# --- calculate direct prediction error
flat_list$DirectPred_Error <- flat_list$solubility_g_100g_log - flat_list$PREDICTION

# --- get frequency tables of each solute and each solvent in original dataset (masterSet)
real_solu_freq <- sort(table(masterSet$Solute),decreasing=T)
real_solv_freq <- sort(table(masterSet$Solvent),decreasing=T)

# --- get frequency tables of each solvent in flat_list - each instance is effectively multiplied by the number of combinations it was in
count_test <- sort(table(flat_list$solvent_name),decreasing=T)

# --- misc prep
solvent_list <- unique(flat_list$solvent_name)
solute_list <- unique(flat_list$solute_name)
solvent_combos_all <- names(count_test2)

# --- WARNING!
# --- Some solvents have spaces in their name... this needs to be addressed for the count to work
# --- Workaround: identify then search/replace each solvent with a space in name

bad_solv_names <- as.character(solvent_list[grepl(" ",  solvent_list)])
new_solv_names <- gsub(" ", "", bad_solv_names)

for (x in 1:length(bad_solv_names)) {
  print(x)
  flat_list$Training_Solvents <- gsub(bad_solv_names[x], new_solv_names[x], flat_list$Training_Solvents)
}

flat_list$Training_Solvents_split <- strsplit(flat_list$Training_Solvents, " +")
flat_list$Training_Solvents_count <- sapply(flat_list$Training_Solvents_split,function(x) {
  x <- ifelse(x[1] == "none", 0, length(unlist(x)))
  x
})

# --- 1 solvent
# by MAE
flat_list_dripfeed1 <- flat_list[flat_list$Training_Solvents_count == 1,]
flat_list_dripfeed1$DirectPred_Error_abs <- abs(flat_list_dripfeed1$DirectPred_Error)
grouped_means_abs_dripfeed1 <- aggregate(flat_list_dripfeed1$DirectPred_Error_abs, list(Training_Solvents = flat_list_dripfeed1$Training_Solvents), mean)

dripfeed1_freq <- sort(table(flat_list_dripfeed1$Training_Solvents),decreasing=T)
temp <- as.data.frame(dripfeed1_freq)
names(temp)[names(temp) == 'Var1'] <- 'Training_Solvents'
grouped_means_abs_dripfeed1 <- merge(grouped_means_abs_dripfeed1, temp, by = "Training_Solvents")
# by RMSE
list_of_rmse <- lapply(split(flat_list_dripfeed1, flat_list_dripfeed1$Training_Solvents), function (z) { sqrt(mean((z$DirectPred_Error)^2)) } )
rmse_vals <- do.call("rbind", list_of_rmse)
grouped_rmse_dripfeed1 <- data.frame(Training_Solvents = row.names(rmse_vals), data.frame(rmse_vals), row.names = NULL)
grouped_rmse_dripfeed1 <- merge(grouped_rmse_dripfeed1, temp, by = "Training_Solvents")

# --- 2 solvents
# by MAE
flat_list_dripfeed2 <- flat_list[flat_list$Training_Solvents_count == 2,]
flat_list_dripfeed2$DirectPred_Error_abs <- abs(flat_list_dripfeed2$DirectPred_Error)
grouped_means_abs_dripfeed2 <- aggregate(flat_list_dripfeed2$DirectPred_Error_abs, list(Training_Solvents = flat_list_dripfeed2$Training_Solvents), mean)

dripfeed2_freq <- sort(table(flat_list_dripfeed2$Training_Solvents),decreasing=T)
temp <- as.data.frame(dripfeed2_freq)
names(temp)[names(temp) == 'Var1'] <- 'Training_Solvents'
grouped_means_abs_dripfeed2 <- merge(grouped_means_abs_dripfeed2, temp, by = "Training_Solvents")
# by RMSE
list_of_rmse <- lapply(split(flat_list_dripfeed2, flat_list_dripfeed2$Training_Solvents), function (z) { sqrt(mean((z$DirectPred_Error)^2)) } )
rmse_vals <- do.call("rbind", list_of_rmse)
grouped_rmse_dripfeed2 <- data.frame(Training_Solvents = row.names(rmse_vals), data.frame(rmse_vals), row.names = NULL)
grouped_rmse_dripfeed2 <- merge(grouped_rmse_dripfeed2, temp, by = "Training_Solvents")

# --- 3 solvents
# by MAE
flat_list_dripfeed3 <- flat_list[flat_list$Training_Solvents_count == 3,]
flat_list_dripfeed3$DirectPred_Error_abs <- abs(flat_list_dripfeed3$DirectPred_Error)
grouped_means_abs_dripfeed3 <- aggregate(flat_list_dripfeed3$DirectPred_Error_abs, list(Training_Solvents = flat_list_dripfeed3$Training_Solvents), mean)

dripfeed3_freq <- sort(table(flat_list_dripfeed3$Training_Solvents),decreasing=T)
temp <- as.data.frame(dripfeed3_freq)
names(temp)[names(temp) == 'Var1'] <- 'Training_Solvents'
grouped_means_abs_dripfeed3 <- merge(grouped_means_abs_dripfeed3, temp, by = "Training_Solvents")
# by RMSE
list_of_rmse <- lapply(split(flat_list_dripfeed3, flat_list_dripfeed3$Training_Solvents), function (z) { sqrt(mean((z$DirectPred_Error)^2)) } )
rmse_vals <- do.call("rbind", list_of_rmse)
grouped_rmse_dripfeed3 <- data.frame(Training_Solvents = row.names(rmse_vals), data.frame(rmse_vals), row.names = NULL)
grouped_rmse_dripfeed3 <- merge(grouped_rmse_dripfeed3, temp, by = "Training_Solvents")

# --- 4 solvents
# by MAE
flat_list_dripfeed4 <- flat_list[flat_list$Training_Solvents_count == 4,]
flat_list_dripfeed4$DirectPred_Error_abs <- abs(flat_list_dripfeed4$DirectPred_Error)
grouped_means_abs_dripfeed4 <- aggregate(flat_list_dripfeed4$DirectPred_Error_abs, list(Training_Solvents = flat_list_dripfeed4$Training_Solvents), mean)

dripfeed4_freq <- sort(table(flat_list_dripfeed4$Training_Solvents),decreasing=T)
temp <- as.data.frame(dripfeed4_freq)
names(temp)[names(temp) == 'Var1'] <- 'Training_Solvents'
grouped_means_abs_dripfeed4 <- merge(grouped_means_abs_dripfeed4, temp, by = "Training_Solvents")
# by RMSE
list_of_rmse <- lapply(split(flat_list_dripfeed4, flat_list_dripfeed4$Training_Solvents), function (z) { sqrt(mean((z$DirectPred_Error)^2)) } )
rmse_vals <- do.call("rbind", list_of_rmse)
grouped_rmse_dripfeed4 <- data.frame(Training_Solvents = row.names(rmse_vals), data.frame(rmse_vals), row.names = NULL)
grouped_rmse_dripfeed4 <- merge(grouped_rmse_dripfeed4, temp, by = "Training_Solvents")


# --- get a mean prediction for each case, per dripfeed size
flat_list_dripfeed0 <- flat_list[flat_list$Training_Solvents_count == 0,]
flat_list_dripfeed0$DirectPred_Error_abs <- abs(flat_list_dripfeed0$DirectPred_Error)
flat_list_dripfeed1 <- flat_list[flat_list$Training_Solvents_count == 1,]
flat_list_dripfeed2 <- flat_list[flat_list$Training_Solvents_count == 2,]
flat_list_dripfeed3 <- flat_list[flat_list$Training_Solvents_count == 3,]
flat_list_dripfeed4 <- flat_list[flat_list$Training_Solvents_count == 4,]

# this generates a results file for a single dripfeed size, e.g. dripfeed = 0.
##### edit these lines accordingly to write a file for each of the above
temp_outfile_path <- "path/to/outputdir/dripfeed__df0_mean_pred_per_case.txt"
temp <- flat_list_dripfeed0 %>%
###### nothing to edit below here
  group_by(solvent_name, solute_name) %>%
  summarise_at(c("PREDICTION"), mean)
temp2 <- left_join(temp, masterSet[,1:20], by = c("solvent_name", "solute_name"))

r2_COSMO <- r2_general(temp2$cosmo_solubility_g_100g_log, temp2$solubility_g_100g_log)
r2_RF <- r2_general(temp2$PREDICTION, temp2$solubility_g_100g_log)
rmse_COSMO <- sqrt(mean((temp2$cosmo_solubility_g_100g_log - temp2$solubility_g_100g_log)^2))
mae_COSMO <- mean(abs(temp2$cosmo_solubility_g_100g_log - temp2$solubility_g_100g_log))
rmse_RF <- sqrt(mean((temp2$PREDICTION - temp2$solubility_g_100g_log)^2))
mae_RF <- mean(abs(temp2$PREDICTION - temp2$solubility_g_100g_log))
count_imp <- length(temp2$PREDICTION[abs(temp2$PREDICTION - temp2$solubility_g_100g_log) < abs(temp2$cosmo_solubility_g_100g_log - temp2$solubility_g_100g_log)])
frac_imp <- count_imp/nrow(temp2)

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

write.table(temp2, temp_outfile_path, row.names = FALSE, quote = FALSE, sep = "\t" )


# --- Get progression for one specific solvent sequence # ---

flat_list_dripfeed0 <- flat_list[flat_list$Training_Solvents_count == 0,]
flat_list_dripfeed0$DirectPred_Error_abs <- abs(flat_list_dripfeed0$DirectPred_Error)
t0 <- flat_list_dripfeed0[flat_list_dripfeed0$solute_name == "ibuprofen",]
t1 <- flat_list_dripfeed1[flat_list_dripfeed1$solute_name == "ibuprofen" & flat_list_dripfeed1$Training_Solvents == "Ethanol",]
t2 <- flat_list_dripfeed2[flat_list_dripfeed2$solute_name == "ibuprofen" & flat_list_dripfeed2$Training_Solvents == "Acetone Ethanol",]
t3 <- flat_list_dripfeed3[flat_list_dripfeed3$solute_name == "ibuprofen" & flat_list_dripfeed3$Training_Solvents == "Acetone Ethanol 2-Propanol ",]
t4 <- flat_list_dripfeed4[flat_list_dripfeed4$solute_name == "ibuprofen" & flat_list_dripfeed4$Training_Solvents == "Acetone Ethanol 2-Propanol Heptane",]

t <- full_join(t0,t1)
t <- full_join(t,t2)
t <- full_join(t,t3)
t <- full_join(t,t4)
t$Training_Solvents_split <- NULL

t_path <- "path/to/outputdir/ibu_0_eth_ace_2prop_hex.txt"
write.table(t, t_path, row.names = FALSE, quote = FALSE, sep = "\t" )



