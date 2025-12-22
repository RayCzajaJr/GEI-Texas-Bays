library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(gridExtra)
library(grid)
library(dsem)
library(ggpubr)
library(ggraph)
library(phylopath)
library(ggdag)
library(readxl)
library(qgraph)
library(DHARMa)



########## Get DHARMA plots and df 
########## Assess prediction skill (via a custom loo function) for each lowest AIC model
########## Get predictions for each year, variable, minor bay, trophic system and major bay, with subsequent plots



# calculcate df for each of four final models
logLik(fit_semAB_Sciaenid_notrophics)
logLik(fit_semGB_Sciaenid_fullbottomup)
logLik(fit_semGB_Pred_fullbottomup)
logLik(fit_semAB_Pred_fulltopdown)

# calculcate df for each of all  models
logLik(fit_semAB_Sciaenid_notrophics)
logLik(fit_semAB_Sciaenid_abiotic)
logLik(fit_semAB_Sciaenid_bottomup_noDD)
logLik(fit_semAB_Sciaenid_fullbottomup)
logLik(fit_semAB_Sciaenid_topdown_noDD)
logLik(fit_semAB_Sciaenid_fulltopdown)

logLik(fit_semAB_Pred_notrophics)
logLik(fit_semAB_Pred_abiotic)
logLik(fit_semAB_Pred_bottomupnoDD)
logLik(fit_semAB_Pred_fullbottomup)
logLik(fit_semAB_Pred_topdown_noDD)
logLik(fit_semAB_Pred_fulltopdown)

logLik(fit_semGB_Sciaenid_notrophics)
logLik(fit_semGB_Sciaenid_abiotic)
logLik(fit_semGB_Sciaenid_bottomup_noDD)
logLik(fit_semGB_Sciaenid_fullbottomup)
logLik(fit_semGB_Sciaenid_topdown_noDD)
logLik(fit_semGB_Sciaenid_fulltopdown)

logLik(fit_semGB_Pred_notrophics)
logLik(fit_semGB_Pred_abiotic)
logLik(fit_semGB_Pred_bottomupnoDD)
logLik(fit_semGB_Pred_fullbottomup)
logLik(fit_semGB_Pred_topdown_noDD)
logLik(fit_semGB_Pred_fulltopdown)



# inspect residuals/dharma plots for each of four final models
samples <- loo_residuals(fit_semAB_Sciaenid_notrophics,  what="samples", track_progress=FALSE)
which_use = which(!is.na(AB_Sciaenid_TS))
fitResp = loo_residuals(fit_semAB_Sciaenid_notrophics, what="loo", track_progress=FALSE)[,'est']
simResp = apply(samples, MARGIN=3, FUN=as.vector)[which_use,]
res = DHARMa::createDHARMa(
  simulatedResponse = simResp,
  observedResponse = unlist(AB_Sciaenid_TS)[which_use],
  fittedPredictedResponse = fitResp )
plot(res)

samples <- loo_residuals(fit_semGB_Sciaenid_fulltopdown,  what="samples", track_progress=FALSE)
which_use = which(!is.na(GB_Sciaenid_TS))
fitResp = loo_residuals(fit_semGB_Sciaenid_fulltopdown, what="loo", track_progress=FALSE)[,'est']
simResp = apply(samples, MARGIN=3, FUN=as.vector)[which_use,]
res = DHARMa::createDHARMa(
  simulatedResponse = simResp,
  observedResponse = unlist(GB_Sciaenid_TS)[which_use],
  fittedPredictedResponse = fitResp )
plot(res)

samples <- loo_residuals(fit_semGB_Pred_fullbottomup,  what="samples", track_progress=FALSE)
which_use = which(!is.na(GB_Pred_TS))
fitResp = loo_residuals(fit_semGB_Pred_fullbottomup, what="loo", track_progress=FALSE)[,'est']
simResp = apply(samples, MARGIN=3, FUN=as.vector)[which_use,]
res = DHARMa::createDHARMa(
  simulatedResponse = simResp,
  observedResponse = unlist(GB_Pred_TS)[which_use],
  fittedPredictedResponse = fitResp )
plot(res)

samples <- loo_residuals(fit_semAB_Pred_fulltopdown,  what="samples", track_progress=FALSE)
which_use = which(!is.na(AB_Pred_TS))
fitResp = loo_residuals(fit_semAB_Pred_fulltopdown, what="loo", track_progress=FALSE)[,'est']
simResp = apply(samples, MARGIN=3, FUN=as.vector)[which_use,]
res = DHARMa::createDHARMa(
  simulatedResponse = simResp,
  observedResponse = unlist(AB_Pred_TS)[which_use],
  fittedPredictedResponse = fitResp )
plot(res)


# custome loo/cv function for assessing prediction skills
cross_validation_ts <- function(tsdata, sem, fit_function, loo_function, chunk_size = 4, seed = 42, output_df_name = "loodf_cv", results_df_name = "results_semAB_Pred") {
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Extract time index from the ts object
  time_index <- as.numeric(time(tsdata))  # Converts ts time attribute into a numeric vector
  
  # Convert tsdata into a dataframe (keeping time index as first column)
  tsdata_df <- as.data.frame(tsdata) %>%
    mutate(Time = time_index)  # Add extracted time column
  
  # Identify the start and end year
  min_year <- min(tsdata_df$Time)
  max_year <- max(tsdata_df$Time)
  total_obs <- nrow(tsdata_df)
  
  # Define chunk start years (adjusting the last chunk if needed)
  chunk_starts <- seq(min_year, max_year - 4, by = chunk_size)  # Skip last 4 years for chunking
  
  # Initialize dataframe to store final test-set predictions
  final_loodf_cv <- data.frame()
  
  # Loop over each chunk (1983 to 2018 or as per your data)
  for (start_year in chunk_starts) {
    
    # Determine the chunk size dynamically
    current_chunk_size <- ifelse(start_year == max_year - 2, 3, chunk_size)
    
    # Define the chunk to remove
    test_years <- seq(start_year, start_year + current_chunk_size - 1)
    
    # Ensure we don’t go beyond the dataset range
    test_data <- tsdata_df %>% filter(Time %in% test_years)
    train_data <- tsdata_df %>% filter(!Time %in% test_years)
    
    # Convert train data back to time series (ts object)
    train_data_ts <- ts(train_data %>% select(-Time), start = start(tsdata), frequency = frequency(tsdata))
    
    # Fit the model using the training data (now as ts object)
    fit_model <- fit_function(sem = sem,
                              tsdata = train_data_ts,
                              control = dsem_control(quiet = TRUE))
    
    # Get model predictions using loo_function
    loodf_cv <- loo_function(fit_model, what = "loo", track_progress = FALSE)
    
    # Keep only test set predictions (where those years were missing in training)
    test_set_predictions <- loodf_cv %>%
      filter(Var1 %in% test_years)  # Only keep years that were omitted
    
    # Store these test set predictions
    final_loodf_cv <- bind_rows(final_loodf_cv, test_set_predictions)
  }
  
  # For the last 4 years (2019 to 2022), modify the approach: Predict for each year by excluding that year alone
  last_years <- 2019:2022
  for (year in last_years) {
    # Exclude the current year and create the training set
    train_data_last4 <- tsdata_df %>% filter(Time != year)
    
    # Convert to time series
    train_data_last4_ts <- ts(train_data_last4 %>% select(-Time), start = start(tsdata), frequency = frequency(tsdata))
    
    # Fit the model
    fit_model_last4 <- fit_function(sem = sem,
                                    tsdata = train_data_last4_ts,
                                    control = dsem_control(quiet = TRUE))
    
    # Get model predictions for that year (exclude the year we are predicting)
    loodf_last4 <- loo_function(fit_model_last4, what = "loo", track_progress = FALSE)
    
    # Filter predictions for the current year
    test_set_last4 <- loodf_last4 %>%
      filter(Var1 == year)
    
    # Add to the final predictions
    final_loodf_cv <- bind_rows(final_loodf_cv, test_set_last4)
  }
  
  # Remove NA values before calculating metrics
  final_loodf_cv <- final_loodf_cv %>%
    filter(!is.na(obs) & !is.na(est))
  
  # Calculate RMSE, correlation, and R-squared per variable
  final_metrics <- final_loodf_cv %>%
    group_by(Var2) %>%
    summarise(
      correlation = cor(obs, est, use = "complete.obs"),
      spearman_rho = cor(obs, est, method = "spearman", use = "complete.obs"),
      RMSE = sqrt(mean((obs - est)^2, na.rm = TRUE)),
      R_squared = 1 - (sum((obs - est)^2, na.rm = TRUE) / sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # Assign the predictions and metrics to the global environment with the user-defined names
  assign(output_df_name, final_loodf_cv, envir = .GlobalEnv)  # Assign predictions df
  assign(results_df_name, final_metrics, envir = .GlobalEnv)  # Assign metrics df
  
  # Return both the predictions and the metrics as a list
  return(list(predictions = final_loodf_cv, metrics = final_metrics))
}


# Get predictions for each trophic system and major bay
results <- cross_validation_ts(
  tsdata = AB_Pred_TS,
  sem = semAB_Pred_fulltopdown,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_AB_Pred",  
  results_df_name = "results_semAB_Pred"  
)
print(results_semAB_Pred) 

results <- cross_validation_ts(
  tsdata = GB_Pred_TS,
  sem = semGB_Pred_fullbottomup,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_GB_Pred",  
  results_df_name = "results_semGB_Pred"  
)
print(results_semGB_Pred)

results <- cross_validation_ts(
  tsdata = GB_Sciaenid_TS,
  sem = semGB_Sciaenid_fulltopdown,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_GB_Sciaenid",
  results_df_name = "results_semGB_Sciaenid"  
)
print(results_semGB_Sciaenid)

results <- cross_validation_ts(
  tsdata = AB_Sciaenid_TS,
  sem = semAB_Sciaenid_notrophics,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_AB_Sciaenid",
  results_df_name = "results_semAB_Sciaenid"  
)
print(results_semAB_Sciaenid)


# calculcate spearmen rho and adjusted R2 values
results_semAB_Pred <- as.data.frame(results_semAB_Pred)
results_semGB_Pred <- as.data.frame(results_semGB_Pred)
results_semGB_Sciaenid <- as.data.frame(results_semGB_Sciaenid)
results_semAB_Sciaenid <- as.data.frame(results_semAB_Sciaenid)

add_spearman_rho <- function(pred_df, results_df) {
  results_df$spearman_rho <- NA
  
  for (i in seq_len(nrow(results_df))) {
    var2_value <- results_df$Var2[i]
    
    df_subset <- pred_df %>% 
      filter(Var2 == var2_value)
    
    if (nrow(df_subset) >= 2) {
      rho <- cor(df_subset$obs, df_subset$est, method = "spearman", use = "complete.obs")
      results_df$spearman_rho[i] <- rho
    }
  }
  
  return(results_df)
}

results_semAB_Pred <- add_spearman_rho(loodf_AB_Pred, results_semAB_Pred)
results_semGB_Pred <- add_spearman_rho(loodf_GB_Pred, results_semGB_Pred)
results_semAB_Sciaenid <- add_spearman_rho(loodf_AB_Sciaenid, results_semAB_Sciaenid)
results_semGB_Sciaenid <- add_spearman_rho(loodf_GB_Sciaenid, results_semGB_Sciaenid)

results_semAB_Pred <- results_semAB_Pred %>% select(Var2, correlation, RMSE, R_squared, spearman_rho)
results_semGB_Pred <- results_semGB_Pred %>% select(Var2, correlation, RMSE, R_squared, spearman_rho)
results_semGB_Sciaenid <- results_semGB_Sciaenid %>% select(Var2, correlation, RMSE, R_squared, spearman_rho)
results_semAB_Sciaenid <- results_semAB_Sciaenid %>% select(Var2, correlation, RMSE, R_squared, spearman_rho)

results_semAB_Pred$adjusted_R2 <- NA
results_semGB_Pred$adjusted_R2 <- NA
results_semGB_Sciaenid$adjusted_R2 <- NA
results_semAB_Sciaenid$adjusted_R2 <- NA


calculate_adjusted_r2 <- function(df, predictors_count, observations_count, specific_var2) {
  df_filtered <- df %>%
    filter(Var2 == specific_var2) %>%
    mutate(
      adjusted_R2 = 1 - ((1 - R_squared) * (observations_count - 1)) / (observations_count - predictors_count - 1)
    )
  df$adjusted_R2[df$Var2 == specific_var2] <- df_filtered$adjusted_R2
  return(df)
}

results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_AransasBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_CopanoBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_MesquiteBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_AransasBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_CopanoBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_MesquiteBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMenhaden_AransasBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMenhaden_CopanoBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMenhaden_MesquiteBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMullet_AransasBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMullet_CopanoBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMullet_MesquiteBay")

results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_GalvestonBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_TrinityBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_WestBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_EastBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_GalvestonBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_TrinityBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_WestBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_EastBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Menhaden_GalvestonBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Menhaden_TrinityBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Menhaden_WestBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Menhaden_EastBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Mullet_GalvestonBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Mullet_TrinityBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Mullet_WestBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Mullet_EastBay")

results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "SpottedSeatrout_GalvestonBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "SpottedSeatrout_TrinityBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "SpottedSeatrout_WestBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "SpottedSeatrout_EastBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_GalvestonBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_TrinityBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_WestBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_EastBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "RedDrum_GalvestonBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "RedDrum_TrinityBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "RedDrum_WestBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "RedDrum_EastBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_GalvestonBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_TrinityBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_WestBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_EastBay")

results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "SpottedSeatrout_AransasBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "SpottedSeatrout_CopanoBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "SpottedSeatrout_MesquiteBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_AransasBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_CopanoBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_MesquiteBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "RedDrum_AransasBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "RedDrum_CopanoBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "RedDrum_MesquiteBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_AransasBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_CopanoBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_MesquiteBay")

# exporting all 8 of the prediction DFs for easy future access because the prediction models take a while to run.
library(openxlsx)
desktop_path <- "~/Desktop/2025 DSEM/May/LOO Data Files"
write.xlsx(loodf_GB_Pred, file.path(desktop_path, "loodf_GB_Pred.xlsx"))
write.xlsx(results_semGB_Pred, file.path(desktop_path, "results_semGB_Pred.xlsx"))
write.xlsx(loodf_GB_Sciaenid, file.path(desktop_path, "loodf_GB_Sciaenid.xlsx"))
write.xlsx(results_semGB_Sciaenid, file.path(desktop_path, "results_semGB_Sciaenid.xlsx"))

write.xlsx(loodf_AB_Pred, file.path(desktop_path, "loodf_AB_Pred.xlsx"))
write.xlsx(results_semAB_Pred, file.path(desktop_path, "results_semAB_Pred.xlsx"))
write.xlsx(loodf_AB_Sciaenid, file.path(desktop_path, "loodf_AB_Sciaenid.xlsx"))
write.xlsx(results_semAB_Sciaenid, file.path(desktop_path, "results_semAB_Sciaenid.xlsx"))

# bringing in exported files 
folder_path <- "~/Desktop/2025 DSEM/May/LOO Data Files"
excel_files <- list.files(path = folder_path, pattern = "\\.xlsx?$", full.names = TRUE)
data_list <- lapply(excel_files, read_excel)
names(data_list) <- tools::file_path_sans_ext(basename(excel_files))
list2env(data_list, envir = .GlobalEnv)


########## Time Series Plots


plot_observed_vs_predicted <- function(df, response_var, yaxis_title, plot_title, predictions_color, results_df) {
  # Filter the data for the specified response variable
  plot_data <- df %>% 
    filter(Var2 == response_var) %>%
    select(Var1, obs, est, se)
  # Ensure Var1 (Year) is numeric
  plot_data <- plot_data %>% mutate(Var1 = as.numeric(Var1))
  # Calculate confidence intervals (95% CI) for the predictions
  plot_data <- plot_data %>%
    mutate(
      lower_ci = est - 1.96 * se,
      upper_ci = est + 1.96 * se
    )
  # Extract spearman rho value using base R indexing (avoiding `pull()`)
  rho_value <- results_df$spearman_rho[results_df$Var2 == response_var]
  
  min_ci_value <- min(plot_data$lower_ci, na.rm = TRUE)
  
  # Create the base plot
  p <- ggplot(plot_data, aes(x = Var1)) +
    # Observed values (dots + line)
    geom_point(aes(y = obs), color = "black", size =0.3) +
    geom_line(aes(y = obs), color = "black", linewidth = 0.1) +
    # Predicted values (dots + line)
    geom_point(aes(y = est), color = predictions_color, size=0.5) +
    geom_line(aes(y = est), color = predictions_color, linewidth = 0.3) +
    # Confidence interval shading
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = predictions_color, alpha = 0.25) +
    # Formatting
    scale_x_continuous(breaks = seq(min(plot_data$Var1), max(plot_data$Var1), by = 2)) +  # Show every other year
    labs(
      title = plot_title,
      y = ifelse(is.na(yaxis_title), "", yaxis_title),  # Set y-axis title to "" if NA
      x = "Year"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),  # Rotate x-axis labels
      axis.text.y = element_text(size = 4),
      axis.title.y = element_text(size = 6.7),
      axis.title.x = element_text(size = 5),
      plot.title = element_text(size = 7, hjust = 0.5)
    ) +
    # Add spearman rho annotation in the top right corner
    annotate("text", 
             x = 1999,  # Adjust for positioning near the right edge
             y = min_ci_value * .95, 
             label = paste0("Sp. Rho = ", round(rho_value, 2)), 
             size = 1.9, 
             fontface = "bold", 
             color = "black")
  # If yaxis_title is NA, remove the y-axis title completely
  if (is.na(yaxis_title)) {
    p <- p + theme(axis.title.y = element_blank())
  }
  # If plot title is NA, remove it completely
  if (is.na(plot_title)) {
    p <- p + theme(plot.title = element_blank())
  }
  # Return the plot
  return(p)
}
# Function to remove x-axis text while keeping the optional y-axis title removal
plot_observed_vs_predicted_no_xtext <- function(df, response_var, yaxis_title, plot_title, predictions_color, results_df) {
  p <- plot_observed_vs_predicted(df, response_var, yaxis_title, plot_title, predictions_color, results_df) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),  
    )
  return(p)
}


# Make plots


#AB PREDATORS
AlligatorGar_CopanoBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AlligatorGar_CopanoBay",
  yaxis_title = NA, plot_title = NA,
  predictions_color = "#2a9d8f", results_df = results_semAB_Pred)

AlligatorGar_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AlligatorGar_AransasBay",
  predictions_color = "#2a9d8f",yaxis_title = "Alligator Gar CPUE",
  plot_title = NA, results_df = results_semAB_Pred)

AlligatorGar_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AlligatorGar_MesquiteBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Pred)

BullShark_CopanoBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "BullShark_CopanoBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "Copano Bay",results_df = results_semAB_Pred)

BullShark_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "BullShark_AransasBay",
  predictions_color = "#2a9d8f", yaxis_title = "Bull Shark CPUE",
  plot_title = "Aransas Bay",results_df = results_semAB_Pred)

BullShark_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "BullShark_MesquiteBay",
  predictions_color = "#2a9d8f",yaxis_title = NA,
  plot_title = "Mesquite Bay",results_df = results_semAB_Pred)

Mullet_CopanoBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AllMullet_CopanoBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Pred)

Mullet_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AllMullet_AransasBay",
  predictions_color = "#2a9d8f",yaxis_title = "Mullet CPUE",
  plot_title = NA, results_df = results_semAB_Pred)

Mullet_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AllMullet_MesquiteBay",
  predictions_color = "#2a9d8f",yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Pred)

Menhaden_CopanoBay <- plot_observed_vs_predicted(
  df = loodf_AB_Pred, response_var = "AllMenhaden_CopanoBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Pred)

Menhaden_AransasBay <- plot_observed_vs_predicted(
  df = loodf_AB_Pred, response_var = "AllMenhaden_AransasBay",
  predictions_color = "#2a9d8f", yaxis_title = "Menhaden CPUE",
  plot_title = NA, results_df = results_semAB_Pred)

Menhaden_MesquiteBay <- plot_observed_vs_predicted(
  df = loodf_AB_Pred, response_var = "AllMenhaden_MesquiteBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Pred)


#GB PREDATORS
AlligatorGar_GalvestonBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "AlligatorGar_GalvestonBay",
  yaxis_title = NA, plot_title = NA,
  predictions_color = "#f28482", results_df = results_semGB_Pred)

AlligatorGar_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "AlligatorGar_TrinityBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

AlligatorGar_WestBay <-plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "AlligatorGar_WestBay",
  predictions_color = "#f28482", yaxis_title = "Alligator Gar CPUE",
  plot_title = NA, results_df = results_semGB_Pred)

AlligatorGar_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "AlligatorGar_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

BullShark_GalvestonBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "BullShark_GalvestonBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = "Galveston Bay", results_df = results_semGB_Pred)

BullShark_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "BullShark_TrinityBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = "Trinity Bay", results_df = results_semGB_Pred)

BullShark_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "BullShark_WestBay",
  predictions_color = "#f28482", yaxis_title = "Bull Shark CPUE",
  plot_title = "West Bay", results_df = results_semGB_Pred)

BullShark_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "BullShark_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = "East Bay", results_df = results_semGB_Pred)

Mullet_GalvestonBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Mullet_GalvestonBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

Mullet_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Mullet_TrinityBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

Mullet_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Mullet_WestBay",
  predictions_color = "#f28482", yaxis_title = "Mullet CPUE",
  plot_title = NA, results_df = results_semGB_Pred)

Mullet_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Mullet_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

Menhaden_GalvestonBay <- plot_observed_vs_predicted(
  df = loodf_GB_Pred, response_var = "Menhaden_GalvestonBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

Menhaden_TrinityBay <- plot_observed_vs_predicted(
  df = loodf_GB_Pred, response_var = "Menhaden_TrinityBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

Menhaden_WestBay <- plot_observed_vs_predicted(
  df = loodf_GB_Pred, response_var = "Menhaden_WestBay",
  predictions_color = "#f28482", yaxis_title = "Menhaden CPUE",
  plot_title = NA, results_df = results_semGB_Pred)

Menhaden_EastBay <- plot_observed_vs_predicted(
  df = loodf_GB_Pred, response_var = "Menhaden_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)


##GB SCIAENID
RedDrum_GalvestonBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "RedDrum_GalvestonBay",
  yaxis_title = NA, plot_title = "Galveston Bay",
  predictions_color = "#f28482", results_df = results_semGB_Sciaenid)

RedDrum_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "RedDrum_TrinityBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = "Trinity Bay", results_df = results_semGB_Sciaenid)

RedDrum_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "RedDrum_WestBay",
  predictions_color = "#f28482", yaxis_title = "Red Drum CPUE",
  plot_title = "West Bay", results_df = results_semGB_Sciaenid)

RedDrum_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "RedDrum_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = "East Bay", results_df = results_semGB_Sciaenid)

SpottedSeatrout_GalvestonBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "SpottedSeatrout_GalvestonBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

SpottedSeatrout_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "SpottedSeatrout_TrinityBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

SpottedSeatrout_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "SpottedSeatrout_WestBay",
  predictions_color = "#f28482", yaxis_title = "Spotted Seatrout CPUE",
  plot_title = NA, results_df = results_semGB_Sciaenid)

SpottedSeatrout_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "SpottedSeatrout_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

BlueCrabSmall_GalvestonBay <- plot_observed_vs_predicted(
  df = loodf_GB_Sciaenid, response_var = "BlueCrabSmall_GalvestonBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

BlueCrabSmall_TrinityBay <- plot_observed_vs_predicted(
  df = loodf_GB_Sciaenid, response_var = "BlueCrabSmall_TrinityBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

BlueCrabSmall_WestBay <- plot_observed_vs_predicted(
  df = loodf_GB_Sciaenid, response_var = "BlueCrabSmall_WestBay",
  predictions_color = "#f28482", yaxis_title = "Blue Crab CPUE",
  plot_title = NA, results_df = results_semGB_Sciaenid)

BlueCrabSmall_EastBay <- plot_observed_vs_predicted(
  df = loodf_GB_Sciaenid, response_var = "BlueCrabSmall_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

Atlanticcroaker_GalvestonBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "Atlanticcroaker_GalvestonBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

Atlanticcroaker_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "Atlanticcroaker_TrinityBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

Atlanticcroaker_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "Atlanticcroaker_WestBay",
  predictions_color = "#f28482", yaxis_title = "Atlantic Croaker CPUE",
  plot_title = NA, results_df = results_semGB_Sciaenid)

Atlanticcroaker_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "Atlanticcroaker_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

##AB SCIAENID
RedDrum_CopanoBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "RedDrum_CopanoBay",
  yaxis_title = NA, plot_title = "Copano Bay",
  predictions_color = "#2a9d8f", results_df = results_semAB_Sciaenid)

RedDrum_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "RedDrum_AransasBay",
  predictions_color = "#2a9d8f", yaxis_title = "Red Drum CPUE",
  plot_title = "Aransas Bay", results_df = results_semAB_Sciaenid)

RedDrum_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "RedDrum_MesquiteBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "Mesquite Bay", results_df = results_semAB_Sciaenid)

SpottedSeatrout_CopanoBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "SpottedSeatrout_CopanoBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

SpottedSeatrout_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "SpottedSeatrout_AransasBay",
  predictions_color = "#2a9d8f", yaxis_title = "Spotted Seatrout CPUE",
  plot_title = NA, results_df = results_semAB_Sciaenid)

SpottedSeatrout_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "SpottedSeatrout_MesquiteBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

BlueCrabSmall_CopanoBay <- plot_observed_vs_predicted(
  df = loodf_AB_Sciaenid, response_var = "BlueCrabSmall_CopanoBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

BlueCrabSmall_AransasBay <- plot_observed_vs_predicted(
  df = loodf_AB_Sciaenid, response_var = "BlueCrabSmall_AransasBay",
  predictions_color = "#2a9d8f", yaxis_title = "Blue Crab CPUE",
  plot_title = NA, results_df = results_semAB_Sciaenid)

BlueCrabSmall_MesquiteBay <- plot_observed_vs_predicted(
  df = loodf_AB_Sciaenid, response_var = "BlueCrabSmall_MesquiteBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

Atlanticcroaker_CopanoBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "Atlanticcroaker_CopanoBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

Atlanticcroaker_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "Atlanticcroaker_AransasBay",
  predictions_color = "#2a9d8f", yaxis_title = "Atlantic Croaker CPUE",
  plot_title = NA, results_df = results_semAB_Sciaenid)

Atlanticcroaker_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "Atlanticcroaker_MesquiteBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

# predator plots
GB_keystone <- grid.arrange(
  BullShark_WestBay, BullShark_GalvestonBay,  BullShark_EastBay, BullShark_TrinityBay,
  AlligatorGar_WestBay, AlligatorGar_GalvestonBay, AlligatorGar_EastBay,AlligatorGar_TrinityBay,
  Mullet_WestBay, Mullet_GalvestonBay,  Mullet_EastBay, Mullet_TrinityBay,
  Menhaden_WestBay, Menhaden_GalvestonBay, Menhaden_EastBay, Menhaden_TrinityBay,
                                        ncol = 4, nrow = 4)

ggsave("GB_keystone.tiff", GB_keystone, dpi = 300, bg = "white",
       width = 170,
       height = 100,
       units = "mm")

# prey plots 
GB_sciaenid <- grid.arrange(
  RedDrum_WestBay, RedDrum_GalvestonBay, RedDrum_EastBay, RedDrum_TrinityBay,
  SpottedSeatrout_WestBay, SpottedSeatrout_GalvestonBay, SpottedSeatrout_EastBay, SpottedSeatrout_TrinityBay,
  Atlanticcroaker_WestBay, Atlanticcroaker_GalvestonBay, Atlanticcroaker_EastBay, Atlanticcroaker_TrinityBay, 
BlueCrabSmall_WestBay, BlueCrabSmall_GalvestonBay, BlueCrabSmall_EastBay, BlueCrabSmall_TrinityBay,
  ncol = 4, nrow = 4)

ggsave("GB_sciaenid.tiff", GB_sciaenid,  dpi = 300, bg = "white",
       width = 170,
       height = 100,
       units = "mm")

AB_keystone <- grid.arrange(
                          BullShark_AransasBay, BullShark_MesquiteBay,BullShark_CopanoBay,
                          AlligatorGar_AransasBay, AlligatorGar_MesquiteBay, AlligatorGar_CopanoBay,
                          Mullet_AransasBay, Mullet_MesquiteBay, Mullet_CopanoBay,
                          Menhaden_AransasBay, Menhaden_MesquiteBay, Menhaden_CopanoBay, 
                          ncol = 3, nrow = 4)

ggsave("AB_keystone.tiff", AB_keystone,  dpi = 300, bg = "white",
       width = 170,
       height = 100,
       units = "mm")


AB_sciaenid  <- grid.arrange(
  RedDrum_AransasBay, RedDrum_MesquiteBay, RedDrum_CopanoBay,
  SpottedSeatrout_AransasBay,SpottedSeatrout_MesquiteBay, SpottedSeatrout_CopanoBay, 
  Atlanticcroaker_AransasBay, Atlanticcroaker_MesquiteBay,Atlanticcroaker_CopanoBay,
  BlueCrabSmall_AransasBay, BlueCrabSmall_MesquiteBay, BlueCrabSmall_CopanoBay,
  ncol = 3, nrow = 4)

ggsave("AB_sciaenid.tiff", AB_sciaenid,  dpi = 300, bg = "white",
       width = 170,
       height = 100,
       units = "mm")



#summary plot for adjusted R2 values 

#merge the data frames into one
squared_df <- bind_rows(results_semAB_Pred, results_semGB_Pred, results_semAB_Sciaenid, results_semGB_Sciaenid)
squared_df <- squared_df %>%
  separate(Var2, into = c("species", "bay"), sep = "_", remove = FALSE) %>%
  mutate(bay = gsub("([a-z])([A-Z])", "\\1 \\2", bay)) %>%  
  mutate(species = case_when(
    species == "AllMullet" ~ "Mullet",
    species == "AllMenhaden" ~ "Menhaden",
    species == "BlueCrabSmall" ~ "BlueCrab",
    TRUE ~ species  
  )) %>%
  filter(!species %in% c("PDSI", "Salinity")) %>%  
  mutate(major_bay = case_when(
    bay %in% c("Aransas Bay", "Copano Bay", "Mesquite Bay") ~ "Aransas Bay",
    bay %in% c("Galveston Bay", "West Bay", "East Bay", "Trinity Bay") ~ "Galveston Bay",
    TRUE ~ NA_character_  
  )) %>%
  mutate(trophic_system = case_when(
    species %in% c("BullShark", "AlligatorGar", "Mullet", "Menhaden") ~ "Keystone Predator",
    species %in% c("RedDrum", "SpottedSeatrout", "BlueCrab", "Atlanticcroaker") ~ "Sciaenid",
    TRUE ~ NA_character_  
  ))

bay_order <- c("Aransas Bay", "Copano Bay", "Mesquite Bay", "Galveston Bay", "Trinity Bay", "East Bay", "West Bay")
squared_df$bay <- factor(squared_df$bay, levels = bay_order)
squared_df$major_bay <- factor(squared_df$major_bay, levels = c("Aransas Bay", "Galveston Bay"))

species_colors <- c(
  "BullShark" = "#f28482",
  "AlligatorGar" = "#e9c46a",
  "RedDrum" = "#f28482",
  "SpottedSeatrout" = "#e9c46a",
  "BlueCrab" = "#2a9d8f",
  "AtlanticCroaker" = "#90a955",
  "Mullet" = "#2a9d8f",
  "Menhaden" = "#90a955"
)

species_shapes <- c(
  "BullShark" = 21,
  "AlligatorGar" = 21,
  "RedDrum" = 21,
  "SpottedSeatrout" = 21,
  "BlueCrab" = 24,
  "AtlanticCroaker" = 24,
  "Mullet" = 24,
  "Menhaden" = 24
)

squared_df <- squared_df %>%
  mutate(species = ifelse(species == "Atlanticcroaker", "AtlanticCroaker", species))

# Add a new size variable
squared_df <- squared_df %>%
  mutate(point_size = ifelse(species %in% c("BlueCrab", "AtlanticCroaker", "Mullet", "Menhaden"), 6, 9))
# 6 = smaller triangles, 9 = default circles

bay_order <- c(
  "Aransas Bay",
  "Mesquite Bay",
  "Copano Bay",
  "West Bay",
  "Galveston Bay",
  "East Bay",
  "Trinity Bay"
)

# Keystone Predator plot
plot_keystone <- ggplot(
  squared_df %>% filter(trophic_system == "Keystone Predator"), 
  aes(x = bay, y = spearman_rho, fill = species, shape = species, size = point_size)
) +
  geom_point(alpha = 0.7, color = "black", stroke = 1.2) +   # size is mapped from aes
  scale_size_identity() +   # uses the size values literally
  scale_fill_manual(values = species_colors) +
  scale_shape_manual(values = species_shapes) +
  scale_y_continuous(limits = c(0.10, 0.90), breaks = seq(0.10, 0.90, by = 0.1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 12)) +
  guides(
    fill = guide_legend(override.aes = list(size = 6)),    # adjust legend point sizes
    shape = guide_legend(override.aes = list(size = 6)))+
   scale_x_discrete(limits = bay_order)+
  labs(title = "Keystone Predator Trophic System", x = "Bay", y = "Spearman Rho")

# Sciaenid plot (same approach)
plot_sciaenid <- ggplot(
  squared_df %>% filter(trophic_system == "Sciaenid"), 
  aes(x = bay, y = spearman_rho, fill = species, shape = species, size = point_size)
) +
  geom_point(alpha = 0.7, color = "black", stroke = 1.2) +
  scale_size_identity() +
  scale_fill_manual(values = species_colors) +
  scale_shape_manual(values = species_shapes) +
  scale_y_continuous(limits = c(0.10, 0.90), breaks = seq(0.10, 0.90, by = 0.1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 12)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 6)),    # adjust legend point sizes
    shape = guide_legend(override.aes = list(size = 6)))+
  scale_x_discrete(limits = bay_order)+
  labs(title = "Sciaenid Trophic System", x = "Bay", y = "Spearman Rho")

# Arrange
rhoplot <- grid.arrange(plot_keystone, plot_sciaenid, ncol = 2)
ggsave("rhoplot.png", rhoplot, dpi = 150, bg = "white",
       width = 2000,
       height = 1000,
       units = "px")

