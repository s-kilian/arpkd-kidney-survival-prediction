# This programs estimates performances of different models on the development data set and defines final risk scores for the methods "Cox regression" and "Random survival forest".
# Also, the performances of the final risk scores on the development data set are assessed.

# Load packages
library(tidyverse)
library(survival)
library(mice)
library(Hmisc)
library(here)

# Source functions for the development
source(here("prediction-functions.R"))

# Load development data set (excluded for privacy issues)
load(file = "dev-set.Rdata")
df.dev %>%
  select(all_of(c(var.set$endpoint.vars, var.set$predictors))) ->
  df.dev

# Load set of variables
load("var.set.Rdata")

# Extract reference levels of variables
summarise(
  .data = df.dev,
  across(
    everything(),
    ~levels(.)[1]
  )
) %>%
  as.vector() %>%
  unlist() ->
  reference.levels

# Development data set
df <- df.dev

# Method to aggregate performances
aggr.perf.method <- "mean"

# Performance metric
performance.metric <- "harrells.c"

# k for cross validation
k <- 10

# Export imputation
Imputation.export <- list(
  fun = Impute,
  args = list(
    df = df %>%
      select(-starts_with("RRT")),
    maxit = 20,
    m = 50,
    seed = 20220511
  ),
  packages = c("mice")
)
save(
  Imputation.export,
  file = here("renal-pred-imp-1.Rdata")
)

# load imputed data sets named mi.sets
load(here("renal-pred-imp-1-result.Rdata"))

# Add endpoint
lapply(
  mi.sets,
  function(x) mutate(x, RRT.free.time = df$RRT.free.time, RRT.event = df$RRT.event)
) ->
  mi.sets

# Cross validation split
fold.vec <- Split.data.k.fold(nrow(df), k = k, seed = 11052022)

# Count number of events before age 1 in each fold
df %>%
  mutate(
    fold.set = fold.vec
  ) %>%
  group_by(fold.set) %>%
  summarise(
    n.ev.bef.1 = sum(RRT.event & RRT.free.time <= 1)
  ) %>%
  print()

# Define all predictor subsets of maximum length 7
pred.sets <- list()
n.max.pred <- 7
for (pred.set.len in 1:n.max.pred) {
  pred.sets <- c(
    pred.sets,
    combn(
      var.set$predictors,
      pred.set.len,
      simplify = FALSE
    )
  )
}

### Cox regression ####
#### Choice of predictor set

# Performance evaluation of every predictor set with Harrell's C
performances <- c()
for (pred.set in pred.sets) {
  # Estimate performance
  Eval.cv.model.imp(
    mi.sets = mi.sets,
    fold.vec = fold.vec,
    model.function = list(
      fun = Fit.cox.regr,
      args = list(
        time.var = "RRT.free.time",
        event.var = "RRT.event",
        pred.set = pred.set
      )
    ),
    performance.object = list(
      fun = Calc.performance,
      args = list(
        time.var = "RRT.free.time",
        event.var = "RRT.event",
        method = performance.metric
      )
    )
  ) ->
    performance
  
  performances <- c(
    performances,
    performance
  )
}

# Save results
tuning.results <- list(
  pred.sets = pred.sets,
  performances = performances
)
save(
  tuning.results,
  file = here("renal-risk-pred-rs2-performances.Rdata")
)

# Load performances
load(here("renal-risk-pred-rs2-performances.Rdata"))

# Plot performances with respect to number of predictors
data.frame(
  n.predictors = lapply(
    tuning.results$pred.sets,
    length
  ) %>%
    unlist,
  performance = tuning.results$performances
) %>%
  ggplot() +
  geom_boxplot(
    aes(
      x = n.predictors,
      y = performance,
      group = n.predictors
    )
  ) +
  theme_bw()

# Following figure shows the predictor sets of the best models (Harrell's c > 0.697) with a maximum of 5 predictors in descending performance. 
performance.threshold <- 0.697
best.models.selection <- tuning.results$performances > performance.threshold & 
  (sapply(tuning.results$pred.sets, length) <= 5)
best.models.predictor.sets <- tuning.results$pred.sets[best.models.selection]
best.models.performance <- tuning.results$performances[best.models.selection]

data.frame() -> df.h
for(i in 1:length(best.models.predictor.sets)){
  df.h %>% bind_rows(
    data.frame(
      id.old = i,
      predictors = best.models.predictor.sets[[i]],
      performance = best.models.performance[i]
    )
  ) ->
    df.h
}
df.h %>%
  arrange(desc(performance), id.old) %>%
  group_by(id.old) %>%
  mutate(
    id = cur_group_id()
  ) %>%
  ungroup() %>%
  mutate(
    id = factor(as.numeric(factor(id, levels = unique(id))))
  ) ->
  df.hh

df.hh %>%
  mutate(
    predictors = factor(
      predictors,
      levels = df.hh %>%
        group_by(predictors) %>%
        summarise(n = n()) %>%
        arrange(desc(n)) %>%
        pull(predictors) %>%
        unique()
    )
  ) %>%
  ggplot() +
  geom_point(
    aes(
      x = predictors,
      y = fct_rev(id)
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  ylab("Model") +
  expand_limits(x = 16) +
  annotate(
    geom = "text",
    y = df.hh %>% select(-predictors) %>% unique() %>% pull(id),
    label = df.hh %>% select(-predictors) %>% unique() %>% pull(performance) %>% round(3),
    x = 15
  )

# Following list shows the coefficients for the best models.
for (id. in unique(df.hh$id)) {
  cat("\n\nModel: ", id., "\n")
  lapply(
    mi.sets,
    function(x) Fit.cox.regr(
      df = x,
      time.var = "RRT.free.time",
      event.var = "RRT.event",
      pred.set = df.hh %>% filter(id == id.) %>% pull(predictors)
    )
  ) ->
    model.imp.list
  sapply(
    model.imp.list,
    function(x) x$args$cox.coefs
  ) %>%
    rowMeans() ->
    model.pooled.coefs
  
  print(model.pooled.coefs)
}

# Based on performances and clinical expertise, define final predictor set
final.pred.set <- c(
  "Variant.2",
  "prenatal_oligo_anhydramnios_gest_age_cat",
  "perinatal_gestational_age_cat.3",
  "art.hypertension.2m",
  "lab_blood_crea.2m"
)

##### Calculate final model
# Fit final model on every imputation set and pool coefficients
lapply(
  mi.sets,
  function(x) Fit.cox.regr(
    df = x,
    time.var = "RRT.free.time",
    event.var = "RRT.event",
    pred.set = final.pred.set
  )
) ->
  final.model.imp.list
sapply(
  final.model.imp.list,
  function(x) x$args$cox.coefs
) %>%
  rowMeans() ->
  final.model.pooled.coefs

# The score can be linearly transformed for better application in practice:
# Select reference levels of predictors
reference.level.selection <- reference.levels[names(reference.levels) %in% final.pred.set]

# Set coefficients of reference levels to 0
reference.coefs <- rep(0, length(reference.level.selection))
names(reference.coefs) <- paste0(names(reference.level.selection), reference.level.selection)

# Complete coefficients of final models with reference levels
final.model.pooled.coefs.completed <- 
  c(
    final.model.pooled.coefs,
    reference.coefs
  )
final.model.pooled.coefs.completed <- final.model.pooled.coefs.completed[order(names(final.model.pooled.coefs.completed))]
names(final.model.pooled.coefs.completed) <- names(final.model.pooled.coefs.completed)[order(names(final.model.pooled.coefs.completed))]

# consider a transformation of the model coefficients
tf.list <- list(
  function(x) x,
  function(x) round((x + c(0, 0, rep(0.575, 3), rep(0, 4), rep(1.45, 7)))*5.517)
)

# calculate performances of model without and with transformed coefficients
performances <- c()
for (tf in tf.list) {
  sapply(
    mi.sets,
    function(x){
      Calc.performance(
        prediction.model = list(
          fun = Calc.cox.risk.score,
          args = list(
            pred.set = final.pred.set,
            cox.coefs = tf(final.model.pooled.coefs.completed)
          )
        ),             # prediction object
        new.data = x,                     # data to predict
        time.var = "RRT.free.time",           # name of time variable
        event.var = "RRT.event",          # name of event variable
        method = performance.metric    # performance metric
      )
    }
  ) %>%
    Aggregate.performances(
      method = "mean"
    ) ->
    res
  
  performances <- c(performances, res)
}

# Show distribution of calculated scores for original and transformed model:
# Calculate predictions
final.model.predictions.list <- list()
for(tf in tf.list){
  sapply(
    mi.sets,
    function(x){
      do.call(
        what = Calc.cox.risk.score,
        args = list(
          new.data = x,
          pred.set = final.pred.set,
          cox.coefs = tf(final.model.pooled.coefs.completed)
        )
      )
    }
  ) %>%
    as.data.frame() %>%
    pivot_longer(
      cols = everything(),
      names_to = "mi.set",
      values_to = "prediction"
    ) ->
    final.model.predictions
  
  final.model.predictions.list <- c(
    final.model.predictions.list,
    list(final.model.predictions)
  )
}

df.h <- data.frame()
for(i in 1:length(final.model.predictions.list)){
  df.h %>%
    bind_rows(
      final.model.predictions.list[[i]] %>%
        mutate(
          method = i
        )
      ) ->
    df.h
}

# Plot distribution 
df.h %>%
  ggplot(
    aes(
      x = prediction
    )
  ) +
  geom_histogram(
    bins = 100
  ) +
  geom_vline(
    data = ~group_by(., method) %>% summarise(qu = quantile(prediction, c(0.25, 0.5, 0.75))),
    aes(xintercept = qu),
    linetype = "dashed",
    color = "red"
  ) +
  geom_text(
    data = ~group_by(., method) %>% 
      summarise(
        label = round(quantile(prediction, c(0.25, 0.5, 0.75)), 2),
        x = quantile(prediction, c(0.25, 0.5, 0.75)),
        y = rep(1500, 3)
      ),
    aes(
      label = label,
      x = x,
      y = y
    )
  ) +
  theme_bw() +
  ylab("Frequency") +
  facet_wrap(
    .~method,
    scales = "free",
    ncol = 1
  ) +
  scale_x_continuous(
    breaks = 0:20
  )

# Define integer split points for risk groups by distribution of score and clinical expertise
split.points.2 <- c(8, 10, 13)

# Calculate proportions of risk groups
proportions.2 <- sapply(
  split.points.2,
  function(x) mean(df.h$prediction[df.h$method == 2] < x)
)

# Plot distribution of transformed score with new split points
df.h %>%
  filter(method == 2) %>%
  ggplot(
    aes(
      x = prediction
    )
  ) +
  geom_histogram(
    bins = 100
  ) +
  geom_vline(
    xintercept = split.points.2,
    linetype = "dashed",
    color = "red"
  ) +
  annotate(
    geom = "text",
    label = split.points.2,
    x = split.points.2,
    y = rep(1500, 3)
  ) +
  theme_bw() +
  ylab("Frequency")

# Pool Kaplan-Meier curves in risk groups:
# Extract predicted scores for both scores
mi.sets.predictions.list <- lapply(
  final.model.predictions.list,
  function(x){
    lapply(
      1:length(mi.sets),
      function(i){
        df %>%
          mutate(
            prediction = x %>%
              filter(mi.set == i) %>%
              pull(prediction)
          )
      }
    )
  }
)
 
# Pool Kaplan Meier curves for risk groups with transformed scores
Pool.KM(
  df.list = mi.sets.predictions.list[[2]],
  time.var = "RRT.free.time",
  event.var = "RRT.event",
  split.points = split.points.2
) ->
  pooled.KM.2

Plot.KM(pooled.KM.2)

# Export final model
final.model.cox <- list(
  fun = Calc.cox.risk.score,
  args = list(
    pred.set = final.pred.set,
    cox.coefs = tf.list[[2]](final.model.pooled.coefs.completed)
  )
)

save(
  final.model.cox,
  file = here("final-model-cox.Rdata")
)
final.model.cox.cutoffs <- list(
  split.points.2,
  proportions.2
)
save(
  final.model.cox.cutoffs,
  file = here("final-model-cox-cutoffs.Rdata")
)


### Random survival forest ####
#### Hyperparameter tuning
# Following hyperparameters are tuned:
# 
# * number of trees
# * splitting rule
# * number of predictors considered for splitting (as a function of the total number of predictors)
# * number of splitting points considered ("0" means that all are considered)
# * minimum terminal node size

# Tuning of hyperparameters is done only with the full predictor set for computational reasons. 

# number of (randomly chosen) variables to try for splitting each node as a function
# of the total number of variables
mtry.fun.list <- list(
  "all" = function(x) x,
  "sqrt" = function(x) ceiling(x^(1/2)),
  "0.1" = function(x) ceiling(0.1*x),
  "0.3" = function(x) ceiling(0.3*x)
)
other.hyperparameters <- expand.grid(
  ntree = c(100, 500, 1000),
  splitrule = c("logrank", "logrankscore"),
  mtry.fun.name = names(mtry.fun.list),
  nodesize = c(5, 10, 15, 20),
  nsplit = c(0, 5, 10, 15)
)

# Set seeds for parallel computations
set.seed(09062022)
other.hyperparameters %>%
  mutate(
    seed = round(runif(length(ntree), 0, 10000))
  ) ->
  other.hyperparameters.seed

# Full hyperparameter grid performance evaluation
pred.set <- var.set$predictors

# packages needed
library(parallel)
library(foreach)
library(doParallel)

# specify number of cores you want to use
core.number <- 3

# execute some commands for parallel computation
cl <- makeCluster(core.number)
registerDoParallel(cl)

# execute parallel loop
foreach(
  i = 1:nrow(other.hyperparameters.seed),
  .packages = c("randomForestSRC", "dplyr", "Hmisc", "survival"),
  .combine = bind_rows
) %dopar% {
  set.seed(other.hyperparameters.seed[i, "seed"])
  # Estimate performance
  Eval.RSF.model.imp(
    mi.sets = mi.sets,
    model.function = list(
      fun = Fit.RSF,
      args = list(
        time.var = "RRT.free.time",
        event.var = "RRT.event",
        pred.set = pred.set,
        ntree = other.hyperparameters.seed[i, "ntree"],
        mtry = mtry.fun.list[[other.hyperparameters.seed[i, "mtry.fun.name"]]](length(pred.set)),
        nodesize = other.hyperparameters.seed[i, "nodesize"],
        nodedepth = NULL,
        nsplit = other.hyperparameters.seed[i, "nsplit"],
        splitrule = as.character(other.hyperparameters.seed[i, "splitrule"]),
        fun = rfsrc.fast,
        samptype = "swr"
      )
    ),
    performance.object = list(
      fun = Calc.performance,
      args = list(
        time.var = "RRT.free.time",
        event.var = "RRT.event",
        method = performance.metric
      )
    )
  ) ->
    performance

  data.frame(
    i = i,
    performance = performance
  )
} ->
  tuning.results

save(
  tuning.results,
  file = here("renal-risk-pred-rs2-RSF-tuning-performances.Rdata")
)

# Load tuning results
load(here("renal-risk-pred-rs2-RSF-tuning-performances.Rdata"))

# Plot distribution of performances
tuning.results %>%
  ggplot() +
  geom_histogram(
    aes(
      x = performance
    )
  ) +
  theme_bw()

# Plot computation time with respect to performance and number of trees
other.hyperparameters.seed %>%
  mutate(
    i = row_number()
  ) %>%
  left_join(
    tuning.results,
    by = "i"
  ) ->
  tuning.results.parameters

tuning.results.parameters %>%
  ggplot() +
  geom_point(
    aes(
      x = performance,
      y = as.numeric(comp.time),
      color = factor(ntree)
    )
  ) +
  theme_bw() +
  ylab("Computation time (sec)") +
  xlab("Performance") +
  scale_color_discrete(
    name = "Number of trees"
  )

# Plot performances with respect to different hyperparameters
for (var in names(other.hyperparameters)) {
  tuning.results.parameters %>%
    ggplot() +
    geom_boxplot(
      aes(
        x = factor(!!sym(var)),
        y = performance
      )
    ) +
    theme_bw() +
    xlab(var) ->
    plot
  
  print(plot)
}

# Definition of final hyperparameter values
final.hyperparameter.constellation <- list(
  ntree = 500,
  splitrule = "logrank",
  mtry.fun.name = "sqrt",
  nodesize = 5,
  nsplit = 0
)

# Set seeds for RSF
set.seed(09062022)
lapply(
  pred.sets,
  function(x){
    list(
      pred.set = x,
      seed = round(runif(1, 0, 10000))
    )
  }
) ->
  pred.sets.seed

# Calculate performances for each predictor set
pred.set.results <- data.frame()
for(i in 1:length(pred.sets.seed)){
  set.seed(pred.sets.seed[[i]]$seed)
  pred.set <- pred.sets.seed[[i]]$pred.set
  # Estimate performance
  Eval.RSF.model.imp(
    mi.sets = mi.sets,
    model.function = list(
      fun = Fit.RSF,
      args = list(
        time.var = "RRT.free.time",
        event.var = "RRT.event",
        pred.set = pred.set,
        ntree = final.hyperparameter.constellation[["ntree"]],
        mtry = mtry.fun.list[[final.hyperparameter.constellation[["mtry.fun.name"]]]](length(pred.set)),
        nodesize = final.hyperparameter.constellation[["nodesize"]],
        nodedepth = NULL,
        nsplit = final.hyperparameter.constellation[["nsplit"]],
        splitrule = final.hyperparameter.constellation[["splitrule"]],
        fun = rfsrc.fast,
        samptype = "swr"
      )
    ),
    performance.object = list(
      fun = Calc.performance,
      args = list(
        time.var = "RRT.free.time",
        event.var = "RRT.event",
        method = performance.metric
      )
    )
  ) ->
    performance
  
  pred.set.results %>%
    bind_rows(
      data.frame(
        i = i,
        performance = performance
      )
    ) ->
    pred.set.results
} 

save(
  pred.set.results,
  file = here("renal-risk-pred-rs2-RSF-pred-set-performances.Rdata")
)

# Load performances
load(here("renal-risk-pred-rs2-RSF-pred-set-performances.Rdata"))

# Plot performance with respect to number of predictors
data.frame(
  n.predictors = lapply(
    pred.sets.seed,
    function(x) length(x$pred.set)
  ) %>%
    unlist,
  performance = pred.set.results$performance
) %>%
  ggplot() +
  geom_boxplot(
    aes(
      x = n.predictors,
      y = performance,
      group = n.predictors
    )
  ) +
  theme_bw()

# Following figure shows the predictor sets of the best models (performance > 0.695) with a maximum of 5 predictors.
performance.threshold <- 0.695
best.models.selection <- pred.set.results$performance > performance.threshold & 
  (sapply(
    pred.sets.seed,
    function(x) length(x$pred.set)
  ) <= 5)
best.models.predictor.sets <- pred.sets[best.models.selection]
best.models.performance <- pred.set.results$performance[best.models.selection]

data.frame() -> df.h
for(i in 1:length(best.models.predictor.sets)){
  df.h %>% bind_rows(
    data.frame(
      id.old = i,
      predictors = best.models.predictor.sets[[i]],
      performance = best.models.performance[i]
    )
  ) ->
    df.h
}
df.h %>%
  arrange(desc(performance), id.old) %>%
  group_by(id.old) %>%
  mutate(
    id = cur_group_id()
  ) %>%
  ungroup() %>%
  mutate(
    id = factor(as.numeric(factor(id, levels = unique(id))))
  ) ->
  df.hh

df.hh %>%
  mutate(
    predictors = factor(
      predictors,
      levels = df.hh %>%
        group_by(predictors) %>%
        summarise(n = n()) %>%
        arrange(desc(n)) %>%
        pull(predictors) %>%
        unique()
    )
  ) %>%
  ggplot() +
  geom_point(
    aes(
      x = predictors,
      y = fct_rev(id)
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  ylab("Model") +
  expand_limits(x = 16) +
  annotate(
    geom = "text",
    y = df.hh %>% select(-predictors) %>% unique() %>% pull(id),
    label = df.hh %>% select(-predictors) %>% unique() %>% pull(performance) %>% round(3),
    x = 15
  )

# Define best predictor set for final RSF model
best.pred.set <- unlist(pred.sets[best.models.selection & pred.set.results$performance == max(pred.set.results$performance[best.models.selection])])

# Fit RSFs on multiple imputation sets
set.seed(9112022)
lapply(
  mi.sets,
  function(x) Fit.RSF(
    df = x,
    time.var = "RRT.free.time",
    event.var = "RRT.event",
    pred.set = best.pred.set,
    ntree = final.hyperparameter.constellation[["ntree"]],
    mtry = mtry.fun.list[[final.hyperparameter.constellation[["mtry.fun.name"]]]](length(best.pred.set)),
    nodesize = final.hyperparameter.constellation[["nodesize"]],
    nodedepth = NULL,
    nsplit = final.hyperparameter.constellation[["nsplit"]],
    splitrule = final.hyperparameter.constellation[["splitrule"]],
    fun = rfsrc,
    samptype = "swr"
  )
) ->
  best.model.imp.list

# Calculate performance of best predictor set by cross validation for better comparison with Cox model
set.seed(9112022)
Eval.cv.model.imp(
  mi.sets = mi.sets,
  fold.vec = fold.vec,
  model.function = list(
    fun = Fit.RSF,
    args = list(
      time.var = "RRT.free.time",
      event.var = "RRT.event",
      pred.set = best.pred.set,
      ntree = final.hyperparameter.constellation[["ntree"]],
      mtry = mtry.fun.list[[final.hyperparameter.constellation[["mtry.fun.name"]]]](length(best.pred.set)),
      nodesize = final.hyperparameter.constellation[["nodesize"]],
      nodedepth = NULL,
      nsplit = final.hyperparameter.constellation[["nsplit"]],
      splitrule = final.hyperparameter.constellation[["splitrule"]],
      fun = rfsrc,
      samptype = "swr"
    )
  ),
  performance.object = list(
    fun = Calc.performance,
    args = list(
      time.var = "RRT.free.time",
      event.var = "RRT.event",
      method = performance.metric
    )
  )
) ->
  best.model.cv.performance

save(
  best.model.cv.performance,
  file = here("renal-risk-pred-rs2-RSF-best-model-CV-performance.Rdata")
)

pooled.best.rsf.model <- Pool.RSF.models(
  best.model.imp.list
)

save(
  pooled.best.rsf.model,
  file = here("final-model-rsf.Rdata")
)

# Calculate predictions of final RSF model
sapply(
  mi.sets,
  function(x){
    do.call(
      what = pooled.best.rsf.model$fun,
      args = c(
        list(new.data = x),
        pooled.best.rsf.model$args
      )
    )
  }
) %>%
  as.data.frame() %>%
  pivot_longer(
    cols = everything(),
    names_to = "mi.set",
    values_to = "prediction"
  ) ->
  best.model.predictions

save(
  best.model.predictions,
  file = here("renal-risk-pred-rs2-RSF-best-model-predictions.Rdata")
)

# Load performance and predictions of final RSF model
load(here("renal-risk-pred-rs2-RSF-best-model-CV-performance.Rdata"))
load(here("renal-risk-pred-rs2-RSF-best-model-predictions.Rdata"))

# Following figure shows the distribution of the predicted risk score in the development data set
split.points <- quantile(best.model.predictions$prediction, c(0.25, 0.5, 0.75))
final.model.rsf.cutoffs <- list(
  unname(split.points),
  c(0.25, 0.5, 0.75)
)
save(final.model.rsf.cutoffs, file = here("final-rsf-cutoffs.Rdata"))

best.model.predictions %>%
  ggplot(
    aes(
      x = prediction
    )
  ) +
  geom_histogram(
    bins = 100
  ) +
  geom_vline(
    data = ~summarise(., qu = split.points),
    aes(xintercept = qu),
    linetype = "dashed",
    color = "red"
  ) +
  geom_text(
    data = ~summarise(.,
        label = round(split.points, 0),
        x = split.points,
        y = rep(1500, 3)
      ),
    aes(
      label = label,
      x = x,
      y = y
    )
  ) +
  theme_bw() +
  ylab("Frequency")

# Plot survival in quartile risk classes
lapply(
  1:length(mi.sets),
  function(i){
    df %>%
      mutate(
        prediction = best.model.predictions %>%
          filter(mi.set == i) %>%
          pull(prediction)
      )
  }
) ->
  mi.sets.predictions

Pool.KM(
  df.list = mi.sets.predictions,
  time.var = "RRT.free.time",
  event.var = "RRT.event",
  split.points = split.points
) ->
  pooled.KM.2

Plot.KM(pooled.KM.2)
