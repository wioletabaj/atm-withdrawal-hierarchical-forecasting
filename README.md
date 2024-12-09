# Hierarchical Forecasting of ATM Withdrawals

## Comparison of local and global models for hierarchical time series forecasting of ATM withdrawals.

This project explores hierarchical time series forecasting of daily ATM withdrawals across Switzerland. The data is structured hierarchically (ATM → region → country), enabling analysis at multiple levels. Two forecasting approaches were applied:

- **Local Models**: Separate models were built for each time series (individual ATMs, regions, and the entire country) using harmonic regression with Fourier components and ARMA errors to capture complex seasonality (weekly and monthly patterns).
- **Global Model**: A single model for all time series combined, leveraging the LightGBM algorithm to uncover network-wide patterns.


### Directory and File Structure

- harmonic_regression.R - R script used to build, evaluate, and generate forecasts with the harmonic regression model. The forecasts and fitted values were saved as CSV files, which are later read into the Jupyter Notebook for further analysis.

- forecasting_analysis.ipynb - Jupyter Notebook containing the full workflow, including exploratory data analysis, LightGBM model training&forecasting, and performance evaluation of the local and global models.

- forecasting_analysis.html - HTML export of the Jupyter Notebook.

- data/ - This folder contains all the datasets necessary for the analysis:
    - atm_withdrawals_data.csv: Daily ATM withdrawal data in Switzerland.
    - harmonic_reg_fcst_test.csv: Test set forecasts generated using the harmonic regression model (local approach).
    - harmonic_reg_fcst_train.csv: Fitted values for the training set from the harmonic regression model.
      
    *Precomputed forecasts are included to save time and avoid the high computational cost of running the R code repeatedly.
