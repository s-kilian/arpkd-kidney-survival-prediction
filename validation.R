# This program evaluates the final prediction models on the validation data set

# Load packages
library(tidyverse)
library(survival)
library(mice)
library(Hmisc)
library(here)

# Load functions
source(here("prediction-functions.R"))

# Load validation set (excluded for privacy issues)
load(file = "val-set.Rdata")
df.val %>%
  select(all_of(c(var.set$endpoint.vars, var.set$predictors))) ->
  df.val

# Load variable set
load("var.set.Rdata")

# Extract reference levels of variables
reference.levels <- lapply(
  df.val.list,
  function(x){
    summarise(
      .data = x,
      across(
        everything(),
        ~levels(.)[1]
      )
    ) %>%
      as.vector() %>%
      unlist()
  }
)

# Validation data set
df <- df.val

# Export imputation
Imputation.export <- list(
  fun = Impute,
  args = list(
    df = df %>%
      select(-starts_with("RRT")),
    maxit = 20,
    m = 50,
    seed = 20221117
  ),
  packages = c("mice")
)
save(
  Imputation.export,
  file = here("renal-pred-val-imp-1.Rdata")
)

# Load imputed data set
load("renal-pred-val-imp-1-result.Rdata")

# Add endpoint
lapply(
  mi.sets,
  function(x) mutate(x, RRT.free.time = df$RRT.free.time, RRT.event = df$RRT.event)
) ->
  mi.sets

# Validation of Cox risk score ####
# Load final Cox model
load(here("final-model-cox.Rdata"))
final.model <- final.model.cox

# Distribution of predictions
# Histogram
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

mi.sets.predictions %>%
  bind_rows() %>%
  ggplot(
    aes(
      x = prediction
    )
  ) +
  geom_histogram(
    bins = 100
  ) +
  theme_bw() +
  ylab("Frequency") +
  scale_x_continuous(
    breaks = 0:20
  )

# Quantiles
proportions <- seq(0.01, 0.99, by = 0.01)
sapply(
  proportions,
  function(x){
    mi.sets.predictions %>%
      bind_rows() %>%
      pull(prediction) %>%
      quantile(., x)
  }
) ->
  percentiles

data.frame(
  proportion = proportions,
  percentile = percentiles
) %>%
  ggplot() +
  geom_point(
    aes(
      x = percentile,
      y = proportion
    )
  ) +
  xlab("Score") +
  ylab("Proportion") +
  theme_bw()

# Goodness of predictions
# Harrell's c
mi.sets.predictions %>%
  lapply(
    function(mi.set){
      Calc.performance(
        prediction.model = final.model,
        prediction.supplied = T,
        new.data = mi.set,
        time.var = "RRT.free.time",
        event.var = "RRT.event",
        method = "harrells.c",
        stand.err = T
      )
    }
  ) %>%
  Aggregate.performances.incl.se() %>%
  Est.se.to.CI() %>%
  print()

# Uno's C
# Define different values of tau
uno.tau.vec <- c(1, 3, 10, 18)

# Calculate Uno's C for each value of tau
for(uno.tau in uno.tau.vec){
  cat("\ntau = ", uno.tau, "\n")
  mi.sets.predictions %>%
    lapply(
      function(mi.set){
        Calc.performance(
          prediction.model = final.model,
          prediction.supplied = T,
          new.data = mi.set,
          time.var = "RRT.free.time",
          event.var = "RRT.event",
          method = "unos.c",
          uno.tau = uno.tau,
          stand.err = T
        )
      }
    ) %>%
    Aggregate.performances.incl.se() %>%
    Est.se.to.CI() %>%
    print()
}

# Royston's D
mi.sets.predictions %>%
  lapply(
    function(mi.set){
      Calc.performance(
        prediction.model = final.model,
        prediction.supplied = T,
        new.data = mi.set,
        time.var = "RRT.free.time",
        event.var = "RRT.event",
        method = "roystons.d",
        stand.err = T
      )
    }
  ) %>%
  Aggregate.performances.incl.se() %>%
  Est.se.to.CI() %>%
  print()

# Graphical assessment
# This plot shows for every event time (x axis) the rank of the predicted score of the patient at that event time. The rank is calculated by all patients at risk at that event time. A rank value of 1 means that the patient with event at that time had the highest score among those at risk (which indicates a good model). A rank value > 0.5 means that the patient with event had a higher score than most of other patients at risk. Decreasing ranks over time are expected and mean that prediction perform worse for later events.

# Calculate ranks
lapply(
  mi.sets.predictions,
  function(mi.set){
    mi.set %>%
      mutate(prediction = -prediction) %>%
      concordance(
        object = Surv(RRT.free.time, RRT.event) ~ prediction,
        data = .,
        ranks = T
      ) %>%
      `[[`("ranks") %>%
      mutate(
        rank = (rank+1)/2,
        case.id = row_number()
      )
  }
) %>%
  bind_rows(.id = "mi.set") ->
  df.ranks

# Plot ranks over time
df.ranks %>%
  group_by(case.id, time) %>%
  summarise(
    mean.rank = mean(rank)
  ) %>%
  ggplot() +
  geom_point(
    aes(
      x = time,
      y = mean.rank
    )
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    limits = c(0,1)
  ) +
  scale_x_log10() +
  theme_bw()

# Cox score risk groups ####
# Load cutoffs for risk groups
load(here("final-model-cox-cutoffs.Rdata"))
final.model.cutoffs <- final.model.cox.cutoffs

# Plot actual proportions vs. planned proportions in risk groups
mi.sets.predictions %>%
  bind_rows() %>%
  mutate(
    risk.group = cut(
      prediction,
      breaks = c(-Inf, final.model.cutoffs[[1]], Inf),
      labels = 1:4,
      include.lowest = TRUE
    )
  ) %>%
  group_by(risk.group) %>%
  summarise(n = n()) %>%
  mutate(actual.prop = cumsum(n)/sum(n)) %>%
  arrange(risk.group) %>%
  mutate(desired.prop = c(final.model.cutoffs[[2]], 1)) %>%
  pivot_longer(
    cols = ends_with("prop")
  ) %>%
  ggplot() +
  geom_bar(
    aes(
      x = risk.group,
      y = value,
      fill = name
    ),
    stat = "identity",
    position = "dodge"
  ) +
  theme_bw()

# Plot survival within risk groups
Pool.KM(
  df.list = mi.sets.predictions,
  time.var = "RRT.free.time",
  event.var = "RRT.event",
  split.points = final.model.cutoffs[[1]]
) ->
  pooled.KM

Plot.KM(pooled.KM)

# Plot survival times vs. risk groups until age 1
lapply(
  mi.sets.predictions,
  function(mi.set){
    mi.set %>%
      mutate(
        risk.group = as.numeric(cut(
          prediction,
          breaks = c(-Inf, final.model.cutoffs[[1]], Inf),
          labels = 1:4,
          include.lowest = TRUE
        )),
        case.id = row_number()
      )
  }
) %>%
  bind_rows(.id = "mi.set") ->
  df.risk.groups

df.risk.groups %>%
  group_by(case.id, RRT.free.time, RRT.event) %>%
  summarise(
    mean.rg = mean(risk.group)
  ) %>%
  ggplot() +
  geom_point(
    aes(
      x = RRT.free.time,
      y = mean.rg,
      alpha = RRT.event,
      shape = RRT.event
    )
  ) +
  scale_x_continuous(
    limits = c(2/12,1),
    name = "Event/censoring time"
  ) +
  ylab("Mean risk group") +
  theme_bw() +
  theme(
    legend.position = "none"
  ) ->
  plot.event.times

df.risk.groups %>%
  group_by(mi.set, risk.group) %>%
  summarise(
    prop.gr.1 = mean(RRT.free.time > 1)
  ) %>%
  group_by(risk.group) %>%
  summarise(
    prop.gr.1 = mean(prop.gr.1)
  ) %>%
  ggplot() +
  geom_bar(
    aes(
      y = prop.gr.1,
      x = risk.group
    ),
    stat = "identity"
  ) +
  ylab("Proportion of event-free age > 1 year") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  ) ->
  plot.prop.gr.1

ggpubr::ggarrange(
  plot.event.times,
  plot.prop.gr.1
) %>%
  print()

# Validation of RSF ####
# Save imputed data sets to calculate RSF predictions with R script calculate-RSF-predictions.R
save(
  mi.sets,
  file = here("prediction/mi-sets-validation.Rdata")
)

# Load imputed data sets with RSF predictions
load(here("mi-sets-rsf-predictions-validation.Rdata"))

# Distribution of predictions
# Histogram
mi.sets.predictions %>%
  bind_rows() %>%
  ggplot(
    aes(
      x = prediction
    )
  ) +
  geom_histogram(
    bins = 100
  ) +
  theme_bw() +
  ylab("Frequency")

# Quantiles
proportions <- seq(0.01, 0.99, by = 0.01)
sapply(
  proportions,
  function(x){
    mi.sets.predictions %>%
      bind_rows() %>%
      pull(prediction) %>%
      quantile(., x)
  }
) ->
  percentiles

data.frame(
  proportion = proportions,
  percentile = percentiles
) %>%
  ggplot() +
  geom_point(
    aes(
      x = percentile,
      y = proportion
    )
  ) +
  xlab("Score") +
  ylab("Proportion") +
  theme_bw()

# Goodness of predictions
# Harrell's c
mi.sets.predictions %>%
  lapply(
    function(mi.set){
      Calc.performance(
        prediction.model = final.model,
        prediction.supplied = T,
        new.data = mi.set,
        time.var = "RRT.free.time",
        event.var = "RRT.event",
        method = "harrells.c",
        stand.err = T
      )
    }
  ) %>%
  Aggregate.performances.incl.se() %>%
  Est.se.to.CI() %>%
  print()

# Uno's C
# Different values for tau
uno.tau.vec <- c(1, 3, 10, 18)

# Calculate Uno's C for each value of tau
for(uno.tau in uno.tau.vec){
  cat("\ntau = ", uno.tau, "\n")
  mi.sets.predictions %>%
    lapply(
      function(mi.set){
        Calc.performance(
          prediction.model = final.model,
          prediction.supplied = T,
          new.data = mi.set,
          time.var = "RRT.free.time",
          event.var = "RRT.event",
          method = "unos.c",
          uno.tau = uno.tau,
          stand.err = T
        )
      }
    ) %>%
    Aggregate.performances.incl.se() %>%
    Est.se.to.CI() %>%
    print()
}

# Royston's D
mi.sets.predictions %>%
  lapply(
    function(mi.set){
      Calc.performance(
        prediction.model = final.model,
        prediction.supplied = T,
        new.data = mi.set,
        time.var = "RRT.free.time",
        event.var = "RRT.event",
        method = "roystons.d",
        stand.err = T
      )
    }
  ) %>%
  Aggregate.performances.incl.se() %>%
  Est.se.to.CI() %>%
  print()

# Graphical assessment
# This plot shows for every event time (x axis) the rank of the predicted score of the patient at that event time. The rank is calculated by all patients at risk at that event time. A rank value of 1 means that the patient with event at that time had the highest score among those at risk (which indicates a good model). A rank value > 0.5 means that the patient with event had a higher score than most of other patients at risk. Decreasing ranks over time are expected and mean that prediction perform worse for later events.

# Calculate ranks
lapply(
  mi.sets.predictions,
  function(mi.set){
    mi.set %>%
      mutate(prediction = -prediction) %>%
      concordance(
        object = Surv(RRT.free.time, RRT.event) ~ prediction,
        data = .,
        ranks = T
      ) %>%
      `[[`("ranks") %>%
      mutate(
        rank = (rank+1)/2,
        case.id = row_number()
      )
  }
) %>%
  bind_rows(.id = "mi.set") ->
  df.ranks

# Plot ranks over time
df.ranks %>%
  group_by(case.id, time) %>%
  summarise(
    mean.rank = mean(rank)
  ) %>%
  ggplot() +
  geom_point(
    aes(
      x = time,
      y = mean.rank
    )
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    limits = c(0,1)
  ) +
  scale_x_log10() +
  theme_bw()

# RSF score risk groups ####
# Load cutoffs for risk groups
load(here("prediction/final-rsf-cutoffs.Rdata"))
final.model.cutoffs <- final.model.rsf.cutoffs

# Plot actual proportions vs. planned proportions of risk groups
mi.sets.predictions %>%
  bind_rows() %>%
  mutate(
    risk.group = cut(
      prediction,
      breaks = c(-Inf, final.model.cutoffs[[1]], Inf),
      labels = 1:4,
      include.lowest = TRUE
    )
  ) %>%
  group_by(risk.group) %>%
  summarise(n = n()) %>%
  mutate(actual.prop = cumsum(n)/sum(n)) %>%
  arrange(risk.group) %>%
  mutate(desired.prop = c(final.model.cutoffs[[2]], 1)) %>%
  pivot_longer(
    cols = ends_with("prop")
  ) %>%
  ggplot() +
  geom_bar(
    aes(
      x = risk.group,
      y = value,
      fill = name
    ),
    stat = "identity",
    position = "dodge"
  ) +
  theme_bw()

# Plot survival within risk groups
Pool.KM(
  df.list = mi.sets.predictions,
  time.var = "RRT.free.time",
  event.var = "RRT.event",
  split.points = final.model.cutoffs[[1]]
) ->
  pooled.KM

Plot.KM(pooled.KM)

# Plot survival times vs. risk groups until age 1
lapply(
  mi.sets.predictions,
  function(mi.set){
    mi.set %>%
      mutate(
        risk.group = as.numeric(cut(
          prediction,
          breaks = c(-Inf, final.model.cutoffs[[1]], Inf),
          labels = 1:4,
          include.lowest = TRUE
        )),
        case.id = row_number()
      )
  }
) %>%
  bind_rows(.id = "mi.set") ->
  df.risk.groups

df.risk.groups %>%
  group_by(case.id, RRT.free.time, RRT.event) %>%
  summarise(
    mean.rg = mean(risk.group)
  ) %>%
  ggplot() +
  geom_point(
    aes(
      x = RRT.free.time,
      y = mean.rg,
      alpha = RRT.event,
      shape = RRT.event
    )
  ) +
  scale_x_continuous(
    limits = c(2/12,1),
    name = "Event/censoring time"
  ) +
  ylab("Mean risk group") +
  theme_bw() +
  theme(
    legend.position = "none"
  ) ->
  plot.event.times

df.risk.groups %>%
  group_by(mi.set, risk.group) %>%
  summarise(
    prop.gr.1 = mean(RRT.free.time > 1)
  ) %>%
  group_by(risk.group) %>%
  summarise(
    prop.gr.1 = mean(prop.gr.1)
  ) %>%
  ggplot() +
  geom_bar(
    aes(
      y = prop.gr.1,
      x = risk.group
    ),
    stat = "identity"
  ) +
  ylab("Proportion of event-free age > 1 year") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  ) ->
  plot.prop.gr.1

ggpubr::ggarrange(
  plot.event.times,
  plot.prop.gr.1
) %>%
  print()
