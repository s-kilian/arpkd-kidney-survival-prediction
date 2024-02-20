# arpkd-kidney-survival-prediction
R Code for developing and evaluating prediction models for kidney survival of ARPKD patients.

## Files
`prediction-functions.R` defines functions needed for the development and validation of prediction models.  
`data-split.R` contains code to split a total data set into development and validation data set, stratifying for some variables.  
`development.R` contains code to develop prediction models on the development data set.  
`validation.R` contains code to evaluate prediction models on the validation data set.  
`imputation.R` executes multiple imputation of missing data.  
`calculate-RSF-predictions.R` contains code to calculate predicted scores from pooled random survival forests.

## Data
Patient data is excluded due to privacy issues. However, performance data of different predictor sets is available. Furthermore, the final Cox model with cutoffs for risk groups is available. The final RSF model is excluded due to its large size. 

