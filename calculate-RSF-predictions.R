# This programs calculates predictions from pooled RSF model

# import functions
source(here::here("prediction-functions.R"))

# import model
load(here::here("final-model-rsf.Rdata"))
final.model <- pooled.best.rsf.model

# import data
load(here::here("mi-sets-validation.Rdata"))

# apply RSF model to calculate predictions
lapply(
  mi.sets,
  function(x){
    x %>%
      mutate(
        prediction = apply.model(final.model, x),
        case.id = row_number()
      )
  } 
) ->
  mi.sets.predictions

# export predictions
save(
  mi.sets.predictions,
  file = here::here("mi-sets-rsf-predictions-validation.Rdata")
)


