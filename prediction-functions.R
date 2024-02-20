# This program defines functions for developing and evaluating prediction models

# Load packages
library(tidyverse)
library(randomForestSRC)
library(survival)

# Imputation ####
Impute <- function(
    df,                    # data set
    m = 5,                     # number of multiple imputed data sets
    maxit = 5,                # number of iterations
    seed = 01042022
){
  # Calculate imputations
  mice(
    df,
    m = m,
    maxit = maxit,
    seed = seed,
    printFlag = FALSE
  ) ->
    imp
  
  complete(
    imp,
    "all"
  ) -> result
  
  return(result)
}

# Functions to pool and plot Kaplan-Meier estimates ####
# Complementary log-log transformation
cloglog <- function(x){
  log(-log(1-x))
}

# Adjusted mean function that treats Inf as NA
mean.adj <- function(x){
  ifelse(
    Inf %in% x,
    NA,
    mean(x, na.rm = T)
  ) ->
    mean.value
  
  return(
    list(
      mean = mean.value,
      n = sum(!is.na(x))
    )
  )
}

# Adjusted sd function that treats Inf as NA
sd.adj <- function(x){
  ifelse(
    Inf %in% x,
    NA,
    sd(x, na.rm = T)
  ) ->
    sd.value
  
  return(
    list(
      sd = sd.value,
      n = sum(!is.na(x))
    )
  )
}

# Inverse function to cloglog
cloglog.inv <- function(x){
  1-exp(-exp(x))
}

# Calculate variance of transformed Kaplan-Meier estimate by the Delta-method
Compute.transformed.KM.var <- function(
  surv,
  se
){
  ifelse(
    surv == 1,
    0,
    (1/(log(1-surv)*(1-surv))*se)^2
  )
}

# Pool Kaplan-Meier estimates of risk groups defined by split.points
Pool.KM <- function(
    df.list,               # list of multiple imputed data sets
    time.var,              # name of time variable
    event.var,             # name of event variable
    split.points           # split points for group assignment
){
  # Calculate risk groups defined by split.points
  lapply(
    df.list,
    function(df){
      df %>%
        mutate(
          group = cut(
            prediction,
            breaks = c(-Inf, split.points, Inf),
            labels = seq_len(length(split.points)+1)
          )
        )
    }
  ) ->
    df.list.2
  
  # number of risk groups
  n.groups <- length(split.points)+1
  
  # number of imputed data sets
  n.imp <- length(df.list)
  
  # for each group, pool Kaplan-Meier estimates
  pooled.km.list <- list()
  for(i in 1:n.groups){
    # Kaplan-Meier estimate in each imputed data set
    lapply(
      df.list.2,
      function(df){
        survfit(
          Surv(eval(as.symbol(time.var)), eval(as.symbol(event.var))) ~ 1,
          filter(df, group == i)
        )
      }
    ) ->
      km.res.list
    
    # extract event times including 0 and maximum follow up
    lapply(
      km.res.list,
      function(res){
        unique(c(
          0,
          summary(res)$time,
          res$time[length(res$time)]
        ))
      }
    ) ->
      time.list
    
    time.list %>%
      unlist() %>%
      unique() %>%
      sort() ->
      time.vec
    
    # extract survival probabilities
    lapply(
      km.res.list,
      function(res){
        c(
          1,
          summary(res)$surv
        ) %>%
          `[`(., . != 0)
      } 
    ) ->
      surv.list
    
    # extract standard errors for survival probabilities
    lapply(
      km.res.list,
      function(res){
        c(
          0,
          summary(res)$std.err
        ) %>%
          `[`(., !is.nan(.))
      } 
    ) ->
      se.list
    
    # extract survival probabilities at each event time
    lapply(
      time.vec,
      function(time){
        lapply(
          1:n.imp,
          function(i) surv.list[[i]][
            time.list[[i]][1:(length(time.list[[i]])-1)] <= time &
              time < time.list[[i]][2:length(time.list[[i]])]
          ]
        ) 
      }
    ) ->
      surv.at.time.list
    
    # Pool survival probabilities after transformation
    lapply(
      surv.at.time.list,
      function(x){
        lapply(
          x,
          function(xx) ifelse(is.null(xx), NA, xx)
        ) %>%
          unlist() %>%
          cloglog() %>%
          mean.adj() %>%
          `[[`("mean")
      }
    ) %>%
      unlist() ->
      surv.transformed.pooled
    
    # Calculate number of imputed data sets used for pooling
    lapply(
      surv.at.time.list,
      function(x){
        lapply(
          x,
          function(xx) ifelse(is.null(xx), NA, xx)
        ) %>%
          unlist() %>%
          cloglog() %>%
          mean.adj() %>%
          `[[`("n")
      }
    ) %>%
      unlist() ->
      n.used.for.pooling

    # Calculated between imputation standard error of transformed survival probabilities
    lapply(
      surv.at.time.list,
      function(x){
        lapply(
          x,
          function(xx) ifelse(is.null(xx), NA, xx)
        ) %>%
          unlist() %>%
          cloglog() %>%
          sd.adj() %>%
          `[[`("sd")
      }
    ) %>%
      unlist() ->
      se.between.imp
    
    # At each event time, extract standard errors of survival probabilities
    lapply(
      time.vec,
      function(time){
        lapply(
          1:n.imp,
          function(i) se.list[[i]][
            time.list[[i]][1:(length(time.list[[i]])-1)] <= time &
              time < time.list[[i]][2:length(time.list[[i]])]
          ]
        ) 
      }
    ) ->
      se.at.time.list
    
    # Pool standard errors of transformed survival probabilities to estimate within imputation standard error
    lapply(
      1:length(time.vec),
      function(i){
        lapply(
          1:n.imp,
          function(j){
            Compute.transformed.KM.var(
              surv = surv.at.time.list[[i]][[j]],
              se = se.at.time.list[[i]][[j]]
            )
          }
        ) %>%
          unlist() %>%
          mean() %>%
          sqrt()
      }
    ) %>%
    unlist() ->
      se.within.imp
    
    # Estimate total standard error
    se.transformed.pooled <- sqrt(se.within.imp^2 + se.between.imp^2*(1+1/n.used.for.pooling))
    
    # Note: Pooling after maximum follow up of some of the imputation sets can introduce bias and can lead to non-monotonous survival curves
    pooled.km.list <- c(
      pooled.km.list,
      list(list(
        time = time.vec[1:(sum(n.used.for.pooling == n.imp)+1)],
        surv.transf = c(surv.transformed.pooled[1:(sum(n.used.for.pooling == n.imp))], NA),
        se.transf = c(se.transformed.pooled[1:(sum(n.used.for.pooling == n.imp))], NA)
      ))
    )
  }
  
  return(pooled.km.list)
}

# Validate function Pool.KM
if(F){
  # Add predictions to imputed data sets
  lapply(
    1:n.imp,
    function(i){
      df %>%
        mutate(
          prediction = final.model.predictions.list[[2]] %>%
            filter(mi.set == i) %>%
            pull(prediction)
        )
    }
  ) ->
    df.list
  
  # Define split points and time variable
  split.points <- split.points.2
  time.var <- "RRT.free.time"
  event.var <- "RRT.event"
  
  ## Confirm that function gives same result when pooling identical copies:
  # Apply Pooling to identical copies
  Pool.KM(
    df.list = lapply(
      1:10,
      function(x) df.list[[3]]
    ),
    time.var = time.var,
    event.var = event.var,
    split.points = split.points
  ) ->
    res
  
  # back transformed survival probabilities in first risk group
  cloglog.inv(res[[1]]$surv)
  
  # Normal Kaplan-Meier estimate in first risk group
  summary(
    survfit(
      Surv(RRT.free.time, RRT.event) ~ 1,
      data = df.list[[3]] %>% filter(prediction <= 5)
    )
  )$surv
  
  ## Confirm that pooled value at certain time point is correct
  # Surv pooled by hand:
  lapply(
    df.list[4:6],
    function(df){
      survfit(
        Surv(RRT.free.time, RRT.event) ~ 1,
        data = df %>% filter(prediction <= 10)
      )
    }
  ) ->
    model.list
  time.point <- 4
  sapply(
    model.list,
    function(res){
      summary(res, times = time.point)$surv
    }
  ) %>%
    cloglog() ->
    surv.vec.transf
  
  surv.vec.transf %>%
    mean() %>%
    cloglog.inv()
  
  # Pooled by function
  Pool.KM(
    df.list = df.list[4:6],
    time.var = time.var,
    event.var = event.var,
    split.points = c(10)
  ) ->
    pool.res
  cloglog.inv(pool.res[[1]]$surv.transf)[between(pool.res[[1]]$time, 3.76, 4)]
  
  # SE pooled by hand:
  sapply(
    model.list,
    function(res){
      x <- summary(res, times = time.point)
      (x$std.err/(log(1-x$surv)*(1-x$surv)))^2
    }
  ) %>%
    mean() ->
    within.var
  
  var(surv.vec.transf) ->
    between.var
  
  total.var <- within.var + (1+1/3)*between.var
  
  sqrt(total.var)
  
  pool.res[[1]]$se[between(pool.res[[1]]$time, 3.76, 4)]
}

# Plot Kaplan-Meier curves
Plot.KM <- function(
  km.list                 # list over all groups. Every Element is a list with vectors time, surv.transf, and se.transf
){
  # Create data frame for plotting
  plot.df <- data.frame()
  for(i in 1:length(km.list)){
    plot.df %>%
      bind_rows(
        data.frame(
          time = km.list[[i]]$time,
          surv.transf = km.list[[i]]$surv.transf,
          se.transf = km.list[[i]]$se.transf,
          group = factor(i)
        )
      ) ->
      plot.df
  }
  
  # Calculated back transformed survival probabilities and confidence intervals
  plot.df %>%
    mutate(
      surv = cloglog.inv(surv.transf),
      ci.lb = cloglog.inv(surv.transf-1.96*se.transf),
      ci.ub = cloglog.inv(surv.transf+1.96*se.transf)
    ) ->
    plot.df.2
  
  # Plot Kaplan-Meier curves
  plot.df.2 %>%
    ggplot() +
    geom_step(
      aes(
        x = time,
        y = surv,
        color = group
      )
    ) +
    utile.visuals::geom_stepconfint(
      aes(
        x = time,
        ymin = ci.lb,
        ymax = ci.ub,
        fill = group
      ),
      alpha = 0.3
    ) +
    geom_segment(                 # add lines until maximum follow up
      data = plot.df.2 %>%
        group_by(group) %>%
        arrange(group, time) %>%
        summarise(
          xstart =  time[length(time)-1],
          xend = time[length(time)],
          y = surv[length(time)-1]
        ),
        aes(
          x = xstart,
          xend = xend,
          y = y,
          yend = y,
          color = group
        )
    ) +
    theme_bw()
}

# Functions for development of prediction models ####
# Fit Cox regression and return prediction function
Fit.cox.regr <- function(
    df,                 # data set
    time.var,           # name of time variable
    event.var,          # name of event variable
    pred.set            # vector of predictor names
){
  # Fit Cox regression
  coxph(
    formula = as.formula(sprintf(
      "Surv(%s, %s) ~ %s",
      time.var,
      event.var,
      paste(pred.set, collapse = " + ")
    )),
    data = df
  ) ->
    cox.res
  
  # Extract coefficients
  cox.coefs <- coefficients(cox.res)
  
  # Define prediction function
  cox.fit <- function(new.data, pred.set, cox.coefs){
    if(!all(pred.set %in% names(new.data))) error("Needed variables in new.data missing!")
    
    # Calculate risk score dependent on variable type
    sapply(
      pred.set,
      function(pred.name){
        if(is.factor(new.data[[pred.name]])){
          score.add <- ifelse(
            new.data[[pred.name]] == levels(new.data[[pred.name]])[1],
            0,
            unname(cox.coefs[paste0(pred.name, new.data[[pred.name]])])
          )
        }
        if(is.numeric(new.data[[pred.name]])){
          score.add <- new.data[[pred.name]]*unname(cox.coefs[pred.name])
        }
        if(is.logical(new.data[[pred.name]])){
          score.add <- ifelse(
            new.data[[pred.name]],
            unname(cox.coefs[paste0(pred.name, "TRUE")]),
            0
          )
        }
        return(score.add)
      }
    ) %>%
      rowSums() %>%
      return()
  }
  
  return(list(
    fun = cox.fit,
    args = list(
      pred.set = pred.set,
      cox.coefs = cox.coefs
    )
  ))
}

# Fit RSF and return prediction function
Fit.RSF <- function(
  df,                 # data set
  time.var,           # name of time variable
  event.var,          # name of event variable
  pred.set,            # vector of predictor names
  ntree = 500,          # number of trees fitted
  mtry = ceiling(sqrt(length(pred.set))),       # number of variables considered at each node splot
  nodesize = 15,      # minimum size of terminal nodes
  nodedepth = NULL,   # maximumdepth of trees
  nsplit = 10,        # number of (randomly chosen) splitting points considered for split
  splitrule = "logrank",      # splitting rule
  fun = rfsrc,                # function to use (rfsrc or faster approximation rfsrc.fast)
  samptype = "swr"            # bootstrap sample sampled with replacement
){
  # Fit RSF
  fun(
    formula = as.formula(sprintf(
      "Surv(%s, %s) ~ %s",
      time.var,
      event.var,
      paste(pred.set, collapse = " + ")
    )),
    data = as.data.frame(df),
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    nsplit = nsplit,
    splitrule = splitrule,
    samptype = samptype,
    perf.type = "none",
    nodedepth = nodedepth
  ) ->
    rsf.res
  
  # Define prediction function
  rsf.fit <- function(new.data, pred.set, rsf.res){
    if(!all(pred.set %in% names(new.data))) error("Needed variables in new.data missing!")
    
    predict(
      object = rsf.res,
      newdata = as.data.frame(new.data)
    )$predicted %>%
      return()
  }
  
  return(list(
    fun = rsf.fit,
    args = list(
      pred.set = pred.set,
      rsf.res = rsf.res
    )
  ))
}

# Functions to estimate performance ####
# Calculate different performance metrics
Calc.performance <- function(
    prediction.model,             # prediction object
    new.data,                     # data to predict
    prediction.supplied = F,        # indicator if prediction is supplied as variable "prediction" in new.data. Skips application of prediction model. 
    time.var,           # name of time variable
    event.var,          # name of event variable
    method = "harrells.c",    # performance metric
    stand.err = F,           # give out standard error
    ...
){
  other.args <- list(...)
  
  # Calculate predictions if not supplied
  if(!prediction.supplied){
    new.data %>%
      mutate(
        prediction = do.call(
          what = prediction.model$fun,
          args = c(
            list(new.data = new.data),
            prediction.model$args
          )
        )
      ) ->
      new.data
  }
  if(method == "somers.d"){
    # Returns Somers' D for using the predicted score as predictor of endpoint
    rcorr.cens(
      x = -new.data$prediction,
      S = Surv(new.data[[time.var]], new.data[[event.var]])
    )[["Dxy"]] ->
      result
  }
  
  if(method == "harrells.c"){
    # Returns Harrell's C for using the predicted score as predictor of endpoint
    rcorr.cens(
      x = -new.data$prediction,
      S = Surv(new.data[[time.var]], new.data[[event.var]])
    ) ->
      result.all
    result <- result.all[["C Index"]]
    se <- result.all[["S.D."]] 
  }

  if(method == "unos.c"){
    # Returns Uno's C
    survC1::Est.Cval(
      mydata = new.data %>%
        select(
          all_of(c(
            time.var,
            event.var,
            "prediction"
          ))
        ),
      tau = other.args$uno.tau,
      nofit = TRUE
    )$Dhat ->
      result
    
    if(stand.err){
      survC1::Inf.Cval(
        mydata = new.data %>%
          select(
            all_of(c(
              time.var,
              event.var,
              "prediction"
            ))
          ),
        tau = other.args$uno.tau,
        itr = 1000,
        seed = other.args$unos.c.seed
      )$se ->
        se
    }
  }

  if(method == "roystons.d"){
    # Returns Royston and Sauerbrei's D
    coxph(
      formula = as.formula(
        sprintf(
          "Surv(%s, %s) ~ prediction",
          time.var,
          event.var
        )
      ),
      data = new.data
      ) ->
        fit

    royston(fit = fit) ->
      result.all
    result <- result.all[["D"]]
    se <- result.all[["se(D)"]]
  }

  if(method == "cox.coef"){
    # Returns coefficient of Cox regression
    coxph(
      formula = as.formula(
        sprintf(
          "Surv(%s, %s) ~ prediction",
          time.var,
          event.var
        )
      ),
      data = new.data
    ) ->
      fit

    coefficients(summary(fit))[,1] ->
      result
    coefficients(summary(fit))[,3] ->
      se
  }
  
  if(method == "surv.conc"){
    # Returns concordance of package survival
    concordance(
      object = as.formula(
        sprintf(
          "Surv(%s, %s) ~ prediction",
          time.var,
          event.var
        )
      ),
      data = new.data,
      ymax = other.args$conc.ymax
    ) ->
      result.all
    result <- result.all$concordance
    se <- sqrt(result.all$var)
  }
  if(!stand.err) return(result) else return(list(result, se))
}

# Function to aggregate multiple performances
Aggregate.performances <- function(
    performances,
    method,
    na.rm = F
){
  if(method == "mean") return(mean(performances, na.rm = na.rm))
}

# Function to aggregate multiple performances and theier standard errors
Aggregate.performances.incl.se <- function(
  perf.ls.ls            # List of performance objects, each being a list with two entries: point estimate and standard error
){
  sapply(
    perf.ls.ls,
    function(x) x[[1]]
  ) ->
    point.estimates
  
  sapply(
    perf.ls.ls,
    function(x) x[[2]]
  ) ->
    se.estimates
  
  mean(point.estimates, na.rm = T) ->
    estimate
  
  var.between <- var(point.estimates, na.rm = T)
  var.within <- mean(se.estimates^2, na.rm = T)
  var.total <- var.within + (1+1/length(perf.ls.ls))*var.between
  
  return(list(estimate = estimate, se = sqrt(var.total)))
}

# Calculate confidence intervals from estimate and standard error
Est.se.to.CI <- function(
  est.se.list             # list with two entries: estimate and se
){
  return(
    list(
      estimate = est.se.list$estimate,
      conf.int = est.se.list$estimate + c(-1, 1)*1.96*est.se.list$se
    )
  )
}

# Split data set randomly with seed in k parts for cross validation
Split.data.k.fold <- function(
    n,                     # number of data points
    k = 10,
    seed = 11052022
){
  set.seed(seed)
  sample(
    x = rep(1:k, each = ceiling(n/k)),
    size = n,
    replace = FALSE,
  ) %>%
    return()
}

# Estimate performance by cross validation
Eval.cv.model <- function(
    df,                  # data frame
    fold.vec,            # vector of fold indices for cross-validation
    model.function,
    performance.object,   # object to calculate performance from model and data
    aggr.perf.method = "mean",          # method to aggregate performances
    na.rm = F                           # should NAs be removed within aggregation?
  ){
  # performances for each of the k splits
  performances <- c()
  for (i in 1:max(fold.vec)) {
    df.train <- df[fold.vec != i,]
    df.test <- df[fold.vec == i,]
    do.call(
      what = model.function$fun,
      args = c(
        list(df = df.train),
        model.function$args
      )
    ) ->
      fitted.model
    do.call(
      what = performance.object$fun,
      args = c(
        list(
          new.data = df.test,
          prediction.model = fitted.model
        ),
        performance.object$args
      )
    ) ->
      performance
    performances <- c(
      performances,
      performance
    )
  }
  
  # Return aggregated performance
  return(
    Aggregate.performances(
      performances = performances,
      method = aggr.perf.method,
      na.rm = na.rm
    )
  )
}

# Estimate performance by cross validation after multiple imputation
Eval.cv.model.imp <- function(
    mi.sets,                  # list of multiple imputed data sets
    fold.vec,            # vector of fold indices for cross-validation
    model.function,
    performance.object,   # object to calculate performance from model and data
    aggr.perf.method = "mean",          # method to aggregate performances
    na.rm = F                        # should NAs be removed when aggregating?
){
  # Estimate performance by cross validation on each imputed data set and aggregate
  lapply(
    mi.sets,
    Eval.cv.model,
    fold.vec = fold.vec,
    model.function = model.function,
    performance.object = performance.object,
    na.rm = na.rm
  ) %>%
    unlist() %>%
    Aggregate.performances(
      performances = .,
      method = aggr.perf.method
    ) %>%
    return()
}

# Estimate performance of RSF model by using OOB predictions
Eval.RSF.model <-  function(
    df,                  # data frame
    model.function,
    performance.object   # object to calculate performance from model and data
){
  # Fit RSF
  do.call(
    what = model.function$fun,
    args = c(
      list(df = df),
      model.function$args
    )
  ) ->
    fitted.model
  
  # Calculate performance by using OOB predictions
  do.call(
    what = performance.object$fun,
    args = c(
      list(
        new.data = df,
        prediction.model = list(
          fun = function(new.data, rsf.res) return(rsf.res$predicted.oob),
          args = list(
            rsf.res = fitted.model$args$rsf.res
          )
        )
      ),
      performance.object$args
    )
  ) ->
    performance

  return(
    performance
  )
}

# Estimate performance of RSF model after multiple imputation
Eval.RSF.model.imp <- function(
  mi.sets,                  # list of multiple imputed data sets
  model.function,
  performance.object,   # object to calculate performance from model and data
  aggr.perf.method = "mean"          # method to aggregate performances
){
  # Aggregate performances on each imputed data set
  lapply(
    mi.sets,
    Eval.RSF.model,
    model.function = model.function,
    performance.object = performance.object
  ) %>%
    unlist() %>%
    Aggregate.performances(
      performances = .,
      method = aggr.perf.method
    ) %>%
    return()
}

# Returns function that calculates the pooled predictions of multiple RSF models 
Pool.RSF.models <- function(
  model.list                 # List of RSF models
){
  return(
    list(
      fun = function(new.data, model.list){
        sapply(
          model.list,
          function(model){
            do.call(
              what = model$fun,
              args = c(
                list(new.data = new.data),
                model$args
              )
            )
          }
        ) %>%
          rowSums() %>%
          return()
      },
      args = list(
        model.list = model.list
      )
    )
  )
}

# Apply prediction model ####
apply.model <- function(model, new.data){
  do.call(
    what = model$fun,
    args = c(
      list(new.data = new.data),
      model$args
    )
  )
}

# Calculate risk score from Cox coefficients
Calc.cox.risk.score <- function(new.data, pred.set, cox.coefs){
  if(!all(pred.set %in% names(new.data))) error("Needed variables in new.data missing!")
  
  sapply(
    pred.set,
    function(pred.name){
      if(is.factor(new.data[[pred.name]])){
        score.add <- unname(cox.coefs[paste0(pred.name, new.data[[pred.name]])])
      }
      if(is.numeric(new.data[[pred.name]])){
        score.add <- new.data[[pred.name]]*unname(cox.coefs[pred.name])
      }
      if(is.logical(new.data[[pred.name]])){
        score.add <- new.data[[pred.name]]*unname(cox.coefs[paste0(pred.name, "TRUE")])
      }
      return(score.add)
    }
  ) %>%
    rowSums() %>%
    return()
}
