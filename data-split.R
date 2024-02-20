# This program splits the total data set df.total into development data set and validation data set in a ratio fo 3:1.
# The distribution of key stratification variables is assessed to evaluate the split.

# Load packages
library(tidyverse)
library(survival)
library(here)

# load data set df.total (excluded for privacy issues)
load(here("df.total.Rdata"))

# define stratification variables for split
split.strat.vars <- c(
  "variant.1",
  "gest.age.sm.37",
  "ventilation.bin",
  "pren.oligo.anhy",
  "RRT.free.time_quart",
  "RRT.event"
)

## Perform split ####
# Function to iteratively randomize the strata in two groups
Split.group <- function(
  df,               # data frame of one stratum
  n.threshold,      # stratum size threshold to apply within stratum split
  df.rand,          # already randomized data frame
  prop.1 = 0.25,    # proportion of individuals in group 1
  strat.vars       # vector of stratification variables
){
  # size of data frame to be randomized
  n <- nrow(df)
  
  # calculate number of already randomized cases and the proportion in group 1
  if(is.null(df.rand)){
    metrics.rand <- c(
      n.rand = 0,
      prop.rand = prop.1
    )
  } else {
    df.rand %>%
      ungroup() %>%
      summarise(
        n.rand = n(),
        prop.rand = sum(split.group == 1)/n.rand
      ) %>%
      as.vector() ->
      metrics.rand
  }
  
  # split stratum randomly if the stratum is large enough.
  # use minimization for small strata
  if(n >= n.threshold){
    # random split
    # number of cases to be randomized to group 1 to yield desired proportion
    n.1 <- round(prop.1*(metrics.rand[["n.rand"]]+n)-metrics.rand[["prop.rand"]]*metrics.rand[["n.rand"]])
    
    # randomize n.1 cases to group 1
    df.rand %>%
      bind_rows(
        df %>%
          ungroup() %>%
          mutate(
            split.group = sample(
              x = c(
                rep(1, n.1),
                rep(0, n-n.1)
              ),
              size = n,
              replace = FALSE
            )
          )
      ) ->
      df.result
  } else {
    # split by randomization
    df.rand -> df.result
    
    # iteratively allocate cases to groups by minimization principle
    for (i in 1:n) {
      # calculate number of matches for the new case if added to each of the groups
      # allocate case with probability 0.8 to group with less matches
      df.result %>%
        bind_rows(
          df[i,] %>%
            full_join(
              data.frame(
                split.group = c(1, 0)
              ),
              by = character()
            )
        ) %>%
        group_by(split.group) %>%
        summarise(
          across(
            all_of(strat.vars),
            list(n.match = ~sum(!is.na(.[seq_len(n()-1)]) & !is.na(.[n()]) & (.[seq_len(n()-1)] == .[n()])))
          )
        ) %>%
        pivot_longer(
          -split.group
        ) %>%
        group_by(split.group) %>%
        summarise(
          sum.match = sum(value)
        ) %>%
        arrange(split.group) %>%
        summarise(
          rand.prob.1 = 0.5 - 0.3*sign(sum.match[2]/prop.1-sum.match[1]/(1-prop.1))
        ) %>%
        pull(rand.prob.1) ->
        rand.prob.1
      
      df.result %>%
        bind_rows(
          df[i, ] %>%
            mutate(
              split.group = rbinom(
                n = 1,
                size = 1,
                prob = rand.prob.1
              )
            )
        ) ->
        df.result
    }
  }
  return(df.result)
}

# set seed
set.seed(12052022)

# create data set for random split
df.total %>%
  mutate( # divide numeric stratification variables in quartiles
    across(
      where(is.numeric) & sub(x = split.strat.vars, replacement = "", pattern = "_quart"),
      list(
        quart = ~cut(
          .,
          breaks = c(-Inf, quantile(., c(0.25, 0.5, 0.75), na.rm = TRUE), Inf),
          labels = 1:4
        )
      )
    )
  ) %>%
  group_by(
    across(all_of(split.strat.vars))
  ) %>%
  mutate( # create group index and number of cases in strata
    ind = cur_group_id(),
    n.stratum = n()
  ) %>%
  ungroup() %>%
  arrange(desc(n.stratum)) ->
  df.to.rand

# extract strata indices
df.to.rand %>%
  pull(ind) %>%
  unique() ->
  strata

# random split of data
df.rand <- NULL

# iteratively split strata, beginning with the largest
for (i in strata) {
  df.to.rand %>%
    filter(ind == i) %>%
    Split.group(
      df = .,
      n.threshold = 4,
      df.rand = df.rand,
      prop.1 = 0.25,
      strat.vars = split.strat.vars
    ) ->
    df.rand
}

# extract and save development data set and validation data set
df.rand %>%
  filter(split.group == 0) ->
  df.dev

df.rand %>%
  filter(split.group == 1) ->
  df.val

save(df.dev, file = here("dev-set.Rdata"))
save(df.val, file = here("val-set.Rdata"))

## Evaluate split ####

load(here("dev-set.Rdata"))
load(here("val-set.Rdata"))

split.eval.var.sets <-  c(
  "Variant.2",
  "perinatal_gestational_age_cat.3",
  "postnatal_ventilation_2",
  "prenatal_oligo_anhydramnios_gest_age_cat"
)

df.dev %>%
  bind_rows(
    df.val
  ) ->
  df.rand
  
cat("\nProportion in validation set: ", mean(df.rand$split.group), "\n")

# Survival endpoint
km.plot.data <- survfit(Surv(RRT.free.time, RRT.event) ~ split.group, df.rand)
survminer::ggsurvplot(
  fit = km.plot.data,
  xlab = "Age (years)",
  ylab = "Survival probability",
  conf.int = TRUE,
  risk.table = TRUE,
  legend = "right",
  risk.table.height = 0.4
) ->
  plot

print(plot)

# Other endpoints
df.rand %>%
  select(all_of(c(split.eval.var.sets, "split.group"))) %>%
  pivot_longer(
    cols = split.eval.var.sets,
    values_transform = list(value = as.character)
  ) %>%
  group_by(split.group, name, value) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(split.group, name) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot() +
  geom_bar(
    stat = "identity",
    aes(
      x = value,
      y = prop,
      fill = factor(split.group)
    ),
    position = "dodge"
  ) +
  facet_wrap(.~name, scales = "free_x", ) ->
  plot

print(plot)
