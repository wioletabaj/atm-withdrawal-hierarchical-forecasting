# Load required libraries
library(feasts)
library(lubridate)
library(dplyr)
library(hts)
library(fable)
library(fabletools)
library(forecast)
library(tidyr)
library(reshape2)
library(zoo)

# Load data from a CSV file
data <- read.csv('data_to_r.csv')

#--------------------------------------------------------------------------------------------------------
# HARMONIC REGRESSION MODEL SETUP

# Ensure 'date' column is in Date format
data$date <- as.Date(data$date)

# Define explanatory variables for the model (paydays and public holidays)

# Define paydays (dates when paychecks are issued)
paydays <- c("2020-01-24", "2020-02-25", "2020-03-25", "2020-04-24", "2020-05-25", 
             "2020-06-25", "2020-07-24", "2020-08-25", "2020-09-25", "2020-10-23", 
             "2020-11-25", "2020-12-24", "2021-01-25", "2021-02-25", "2021-03-25", 
             "2021-04-23", "2021-05-25", "2021-06-25", "2021-07-23", "2021-08-25", 
             "2021-09-24", "2021-10-25", "2021-11-25", "2021-12-24", "2022-01-25", 
             "2022-02-25", "2022-03-25", "2022-04-25", "2022-05-25", "2022-06-24", 
             "2022-07-25", "2022-08-25", "2022-09-23", "2022-10-25", "2022-11-25", 
             "2022-12-23")

# Define public holidays
public_holidays <- c("2020-01-01", "2020-04-10", "2020-04-13", "2020-05-21", 
                     "2020-06-01", "2020-08-01", "2020-12-25", "2020-12-26", 
                     "2021-01-01", "2021-04-02", "2021-04-05", "2021-05-13", 
                     "2021-05-24", "2021-08-01", "2020-12-25", "2020-12-26",
                     "2022-01-01", "2022-04-15", "2022-04-18", "2022-05-26",
                     "2022-06-06", "2022-08-01", "2020-12-25", "2020-12-26")

# Convert data to a time series table (tsibble) format for easier forecasting operations
data <- as_tsibble(data, key = c(canton, region), index = date)

# Aggregate data by region and canton, then create binary variables for paydays and public holidays
data <- data %>% 
  aggregate_key(region / canton, y = sum(y)) %>% # Aggregate y by region and canton
  mutate(
    # Create binary variables for paydays and public holidays
    payday = as.factor(ifelse(as.Date(date) %in% as.Date(paydays), 1, 0)),
    public_holiday = as.factor(ifelse(as.Date(date) %in% as.Date(public_holidays), 1, 0))
  ) %>%
  as_tsibble() # Ensure the data remains in tsibble format for forecasting

# Split data into training and test sets
df_train <- data %>% filter(date < as.Date("2022-10-01"))
df_test <- data %>% filter(date >= as.Date("2022-10-01"))

#--------------------------------------------------------------------------------------------------------
# MODEL FITTING FOR THE TRAINING SET

#' Generate AIC Values for Harmonic Regression Models for a given time series
#'
#' This function generates a table of AIC (Akaike Information Criterion) values for different (specified by the user)
#' combinations of harmonic components for weekly, monthly, and yearly periods in a harmonic regression model.
#' Function enables the user to pass the number of harmonic components for each period (weekly, monthly, and yearly)
#' to be tested and the explanatory variables to be used. It then returns a data frame with the tested combinations 
#' and their corresponding AIC scores. This table can be used to assess the performance of different combinations 
#' but does not return the best combination automatically.
#' 
#' @param train_ts A tsibble containing the training time series data. The data must include 
#'   a dependent variable `y` and explanatory variables specified by user in `explanatory_vars` parameter.
#' @param weekly_harmonics A vector of integers specifying the number of weekly harmonics to test for 
#'   (default is 1 to 3). This corresponds to the number of Fourier terms for the weekly period (7 days).
#' @param monthly_harmonics A vector of integers specifying the number of monthly harmonics to test for 
#'   (default is 1 to 8). This corresponds to the number of Fourier terms for the monthly period (30.5 days).
#' @param yearly_harmonics A vector of integers specifying the number of yearly harmonics to test for 
#'   (default is 1 to 10). This corresponds to the number of Fourier terms for the yearly period (365.25 days).
#' @param explanatory_vars A character vector containing the names of the explanatory variables to include in 
#'   the model formula (default is `c("payday", "public_holiday")`).
#'
#' @return A data frame with the following columns:
#'   - `k7`: The number of weekly harmonics tested (Fourier terms for period = 7 days).
#'   - `k30`: The number of monthly harmonics tested (Fourier terms for period = 30.5 days).
#'   - `k365`: The number of yearly harmonics tested (Fourier terms for period = 365.25 days).
#'   - `aic`: The AIC value of the model corresponding to the combination of harmonic components.
optimize_harmonics <- function(train_ts, weekly_harmonics = 1:3, monthly_harmonics = 1:8, yearly_harmonics = 1:10, explanatory_vars = c("payday", "public_holiday")) {
  
  # Initialize an empty data frame to store results
  results <- data.frame()
  
  # Loop through combinations of harmonic components for weekly, monthly, and yearly periods
  for (num_weekly_harmonics in weekly_harmonics) {        # Weekly harmonics (period = 7 days)
    for (num_monthly_harmonics in monthly_harmonics) {    # Monthly harmonics (period = 30.5 days)
      for (num_yearly_harmonics in yearly_harmonics) {    # Yearly harmonics (period = 365.25 days)
        
        # Construct the formula dynamically based on explanatory variables
        formula <- as.formula(paste("y ~", paste(explanatory_vars, collapse = " + "), 
                                    "+ PDQ(0, 0, 0) + pdq(d = 0) +",
                                    "fourier(period = 7, K =", num_weekly_harmonics, ") +",
                                    "fourier(period = 30.5, K =", num_monthly_harmonics, ") +",
                                    "fourier(period = 365.25, K =", num_yearly_harmonics, ")"))
        
        # Fit the model with ARIMA and Fourier terms for each frequency combination
        fit <- train_ts %>%
          model(base = ARIMA(formula))
        
        # Store the current combination and the model's AIC score
        result <- data.frame(k7 = num_weekly_harmonics, 
                             k30 = num_monthly_harmonics, 
                             k365 = num_yearly_harmonics, 
                             aic = glance(fit)$AIC)
        results <- dplyr::bind_rows(results, result)
        
        # Print iteration details and AIC score for monitoring
        message(sprintf("k7 = %d, k30 = %d, k365 = %d, AIC = %f", 
                        num_weekly_harmonics, num_monthly_harmonics, num_yearly_harmonics, glance(fit)$AIC))
      }
    }
  }
  # Return the results data frame containing all tested combinations and their AIC scores
  return(results)
}

#-------------------------------------------------------------------------------
#' Find Optimal Harmonic Parameters for Each Time Series
#' 
#' This function finds the optimal harmonic regression parameters (for weekly, monthly, and yearly periods)
#' for each time series in the training dataset based on the lowest AIC score.
#' 
#' @param df_train A tsibble containing the training data with time series for different regions and cantons.
#'  The data must include a dependent variable `y` and explanatory variables.
#' @param unique_keys A data frame containing unique combinations of region and canton.
#' @param weekly_harmonics A vector specifying the range of weekly harmonics to test (default: 1:3).
#' @param monthly_harmonics A vector specifying the range of monthly harmonics to test (default: 1:8).
#' @param yearly_harmonics A vector specifying the range of yearly harmonics to test (default: 1:10).
#' @param explanatory_vars A character vector containing the names of explanatory variables (e.g., "payday", "public_holiday").
#' 
#' @return A data frame with the optimal harmonic parameters (k7, k30, k365) for each time series.
find_optimal_params <- function(df_train, weekly_harmonics = 1:3, monthly_harmonics = 1:8, yearly_harmonics = 1:10, explanatory_vars = c("payday", "public_holiday")) {
  
  # Generate unique region and canton combinations
  unique_keys <- df_train %>% distinct(region, canton)
  
  # Prepare a data frame to store optimal parameters for each time series
  optimal_params_df <- data.frame(k7 = rep(NA, nrow(unique_keys)), 
                                  k30 = rep(NA, nrow(unique_keys)), 
                                  k365 = rep(NA, nrow(unique_keys)),
                                  region = unique_keys$region, 
                                  canton = unique_keys$canton)
  
  # Loop over each unique region and canton combination
  for (i in seq_len(nrow(unique_keys))) {
    
    # Filter the training data for the specific region and canton
    ts <- df_train %>% filter(region == unique_keys$region[i] & canton == unique_keys$canton[i])
    
    # Run optimize_harmonics function to get AIC results for different harmonic components
    aic_results <- optimize_harmonics(ts, 
                                      weekly_harmonics = weekly_harmonics, 
                                      monthly_harmonics = monthly_harmonics, 
                                      yearly_harmonics = yearly_harmonics, 
                                      explanatory_vars = explanatory_vars)
    
    # Extract the combination of k parameters with the lowest AIC score
    optimal_harmonics_params <- aic_results[which.min(aic_results$aic), ]
    
    # Store the best parameters
    optimal_params_df[i, c("k7", "k30", "k365")] <- optimal_harmonics_params[, c("k7", "k30", "k365")]
  }
  
  # Return a data frame with the best k values for each time series
  return(optimal_params_df)
}

# Find the best k parameters for each region and canton in the training set
best_params <- find_optimal_params(df_train, 
                                   weekly_harmonics = 1:3, 
                                   monthly_harmonics = 1:8,
                                   yearly_harmonics = 1:10,
                                   explanatory_vars = c("payday", "public_holiday"))


#-------------------------------------------------------------------------------
# HARMONIC REGRESSION FORECASTING

#' Perform Harmonic Regression Fitting and Forecasting
#'
#' This function fits a harmonic regression model for each unique region-canton combination 
#' using Fourier terms and ARIMA model specifications (using passed harmonic parameters). It then 
#' generates both training and test set forecasts for the fitted models.
#'
#' @param df_train A data frame containing the training data. It must have columns for 
#' region, canton, date, and the response variable `y`.
#' @param df_test A data frame containing the test data with the same structure as `df_train`.
#' @param unique_keys A data frame containing unique combinations of region and canton, 
#' with columns `region` and `canton`.
#' @param optimal_params A data frame containing the optimal harmonic parameters (`k7`, `k30`, `k365`) 
#' for each region-canton combination (Should be the output of the `find_optimal_params` function).
#'
#' @return A list with the following components:
#'   \item{train_forecasts_df}{A data frame containing forecasted values on the training set for each region-canton.}
#'   \item{test_forecasts_df}{A data frame containing forecasted values on the test set for each region-canton.}
#'   \item{list_of_coeffs}{A list of coefficients from the harmonic regression models.}
#'   \item{list_of_AICs}{A list of AIC values for the fitted models for each region-canton.}
perform_harmonic_regression <- function(df_train, df_test, unique_keys, optimal_params) {
  
  # Initialize lists to store coefficients and AIC values
  list_of_coeffs <- list()
  list_of_AICs <- list()
  
  # Dataframes to store the forecasts for training and test sets
  train_forecasts_df <- data.frame()  
  test_forecasts_df <- data.frame()
  
  # Generate unique region and canton combinations
  unique_keys <- df_train %>% distinct(region, canton)
  
  # Loop through each combination of region and canton
  for (i in 1:nrow(unique_keys)) {
    
    # Filter training and test dataset  for the current region and canton
    train_ts <- df_train %>% filter(region == unique_keys$region[i] & canton == unique_keys$canton[i])
    test_ts <- df_test %>% filter(region == unique_keys$region[i] & canton == unique_keys$canton[i])
    
    # Fit the harmonic regression model using ARIMA and Fourier terms
    fit <- train_ts %>%
      model(
        base = ARIMA(y ~ payday + public_holiday +
                       PDQ(0, 0, 0) + pdq(d = 0) +
                       fourier(period = 7, K = best_params$k7[i]) + 
                       fourier(period = 30.5, K = best_params$k30[i]) + 
                       fourier(period = 365.25, K = best_params$k365[i]))
      )
    
    # Store the AIC and coefficients for the fitted model
    list_of_AICs[[i]] <- glance(fit)
    list_of_coeffs[[i]] <- coef(fit)
    
    # Forecast on the train forecasts
    train_forecasts <- fit %>% fitted() %>% as.data.frame()
    train_forecasts_df <- bind_rows(train_forecasts_df, train_forecasts)
    
    # Forecast on the test set
    test_forecasts <- forecast(fit, new_data = test_ts) 
    test_forecasts_df <- bind_rows(test_forecasts_df, as.data.frame(test_forecasts))
  }
  
  # Return the forecasted values and model parameters
  return(list(train_forecasts_df, test_forecasts_df, list_of_coeffs, list_of_AICs))
}

# Call the function to perform harmonic regression and get results
harmonic_regression_result <- perform_harmonic_regression(df_train, df_test, best_params)

train_forecasts_df <- harmonic_regression_result[[1]]
test_forecasts_df <- harmonic_regression_result[[2]]
coeffs <- harmonic_regression_result[[3]]
aics <- harmonic_regression_result[[4]]

# Function to generate unique identifier for each region-canton combination
generate_unique_id <- function(region, canton) {
  case_when(
    is_aggregated(region) & is_aggregated(canton) ~ "total",              # For total aggregation
    is_aggregated(canton) ~ paste("region", region, sep = ""),            # For region-level aggregation
    TRUE ~ paste("region", region, "-canton", canton, sep = "")           # For specific region-canton pair
  )
}

# Merge actual and fitted values for the training set
train_subset <- train_forecasts_df %>%
  left_join(df_train %>% select(date, region, canton, y), 
            by = c("date", "region", "canton")) %>%
  mutate(unique_id = generate_unique_id(region, canton)) %>%
  select(unique_id, ds = date, harmonic_reg = .fitted, y)

# Prepare forecasted values for the test set
test_subset <- test_forecasts_df %>%
  mutate(unique_id = generate_unique_id(region, canton)) %>%
  select(unique_id, ds = date, harmonic_reg = .mean)

# Filter data based on predefined unique IDs (regions and cantons of interest)
train_subset_final <- train_subset %>% filter(unique_id %in% unique_ids)
test_subset_final <- test_subset %>% filter(unique_id %in% unique_ids)

# Save the results as CSV files
write.csv(train_subset_final, 'harmonic_reg_fcst_train.csv', row.names = FALSE)
write.csv(test_subset_final, 'harmonic_reg_fcst_test.csv', row.names = FALSE)