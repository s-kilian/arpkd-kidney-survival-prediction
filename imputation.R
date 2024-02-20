# Generate imputed data sets
# Load imputation settings
load("prediction/renal-pred-val-imp-1.Rdata")

# Load packages
sapply(
  Imputation.export$packages,
  library,
  character.only = TRUE
)

# Generate imputations
do.call(
  what = Imputation.export$fun,
  args = Imputation.export$args
) ->
  mi.sets

# Save imputations
save(
  mi.sets,
  file = "prediction/renal-pred-val-imp-1-result.Rdata"
)
