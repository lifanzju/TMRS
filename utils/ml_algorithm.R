library(glmnet)
library(survival)
library(dplyr)
library(randomForestSRC)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)




##### RSF
RSF_function <- function(TrainingDataset, ValidationDatasets, seed) {
  set.seed(seed)
  fit <- rfsrc(Surv(time, status) ~ .,
    data = TrainingDataset$Analysis_table,
    ntree = 10000, nodesize = 15,
    splitrule = "logrank",
    samptype = "swor",
    importance = T,
    proximity = T,
    forest = T,
    seed = seed
  )
  best <- which.min(fit$err.rate)
  set.seed(seed)
  fit <- rfsrc(Surv(time, status) ~ .,
    data = TrainingDataset$Analysis_table,
    ntree = best, nodesize = 15, ## 该值建议多调整
    splitrule = "logrank",
    importance = T,
    proximity = T,
    forest = T,
    seed = seed
  )
  variable_selected <- var.select.rfsrc(fit)$topvars
  RS <- lapply(ValidationDatasets, function(x) {
    prediction <- predict(fit, newdata = x$Analysis_table)
    risk_score <- as.numeric(prediction$predicted)
    names(risk_score) <- row.names(x$Analysis_table)
    rs <- cbind(x$Survival, risk_score)
    rs <- as.data.frame(rs)
    return(rs)
  })
  C_indexs <- sapply(RS, function(x) {
    cox_regression <- coxph(Surv(time, status) ~ risk_score, data = x)
    c_index <- cox_regression$concordance[6]
    return(c_index)
  })
  output <- list(RS = RS, C_Index = C_indexs, Selected_Variable = variable_selected)
  return(output)
}
##### SurvivalSVM
SurvivalSVM_function <- function(TrainingDataset, ValidationDatasets, seed) {
  set.seed(seed)
  fit <- survivalsvm(Surv(time, status) ~ ., data = TrainingDataset$Analysis_table, gamma.mu = 1)
  RS <- lapply(ValidationDatasets, function(x) {
    prediction <- predict(fit, newdata = x$Analysis_table)
    risk_score <- as.numeric(prediction$predicted)
    names(risk_score) <- row.names(x$Analysis_table)
    rs <- cbind(x$Survival, risk_score)
    rs <- as.data.frame(rs)
    return(rs)
  })
  C_indexs <- sapply(RS, function(x) {
    cox_regression <- coxph(Surv(time, status) ~ risk_score, data = x)
    c_index <- cox_regression$concordance[6]
    return(c_index)
  })
  output <- list(RS = RS, C_Index = C_indexs)
  return(output)
}
##### CoxBoost
Coxboost_function <- function(TrainingDataset, ValidationDatasets, seed) {
  set.seed(seed)
  pen <- optimCoxBoostPenalty(
    time = TrainingDataset$Analysis_table$time, status = TrainingDataset$Analysis_table$status, TrainingDataset$EXP,
    maxstepno = 500, trace = F, parallel = T
  )
  cv.res <- cv.CoxBoost(
    time = TrainingDataset$Analysis_table$time, status = TrainingDataset$Analysis_table$status, TrainingDataset$EXP,
    maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty
  )
  fit <- CoxBoost(
    time = TrainingDataset$Analysis_table$time, status = TrainingDataset$Analysis_table$status, TrainingDataset$EXP,
    stepno = cv.res$optimal.step, penalty = pen$penalty
  )
  variable_selected <- names(coef(fit))[coef(fit) != 0]
  RS <- lapply(ValidationDatasets, function(x) {
    risk_score <- as.numeric(predict(fit, newdata = x$EXP, type = "lp"))
    names(risk_score) <- row.names(x$Analysis_table)
    rs <- cbind(x$Survival, risk_score)
    rs <- as.data.frame(rs)
    return(rs)
  })
  C_indexs <- sapply(RS, function(x) {
    cox_regression <- coxph(Surv(time, status) ~ risk_score, data = x)
    c_index <- cox_regression$concordance[6]
    return(c_index)
  })
  output <- list(RS = RS, C_Index = C_indexs, Selected_Variable = variable_selected)
  return(output)
}
##### GBM
GBM_function <- function(TrainingDataset, ValidationDatasets, seed) {
  set.seed(seed)
  fit <- gbm(
    formula = Surv(time, status) ~ ., data = TrainingDataset$Analysis_table, distribution = "coxph",
    n.trees = 10000,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 10, n.cores = 6
  )
  # find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(
    formula = Surv(time, status) ~ ., data = TrainingDataset$Analysis_table, distribution = "coxph",
    n.trees = best,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 10, n.cores = 6
  )
  RS <- lapply(ValidationDatasets, function(x) {
    risk_score <- as.numeric(predict(fit, newdata = x$Analysis_table, type = "link", n.trees = best))
    names(risk_score) <- row.names(x$Analysis_table)
    rs <- cbind(x$Survival, risk_score)
    rs <- as.data.frame(rs)
    return(rs)
  })
  C_indexs <- sapply(RS, function(x) {
    cox_regression <- coxph(Surv(time, status) ~ risk_score, data = x)
    c_index <- cox_regression$concordance[6]
    return(c_index)
  })
  output <- list(RS = RS, C_Index = C_indexs)
  return(output)
}
##### plsRcox
plsRcox_function <- function(TrainingDataset, ValidationDatasets, seed) {
  set.seed(seed)
  cv.plsRcox.res <- cv.plsRcox(list(x = TrainingDataset$EXP, time = TrainingDataset$Analysis_table$time, status = TrainingDataset$Analysis_table$status),
    nt = ncol(TrainingDataset$EXP), nfold = 10, verbose = F
  )
  nt <- as.numeric(cv.plsRcox.res$lambda.min5)
  set.seed(seed)
  fit <- plsRcox(TrainingDataset$EXP, TrainingDataset$Analysis_table$time, event = TrainingDataset$Analysis_table$status, nt = nt)
  RS <- lapply(ValidationDatasets, function(x) {
    risk_score <- as.numeric(predict(fit, newdata = x$EXP, type = "lp"))
    names(risk_score) <- row.names(x$Analysis_table)
    rs <- cbind(x$Survival, risk_score)
    rs <- as.data.frame(rs)
    return(rs)
  })
  C_indexs <- sapply(RS, function(x) {
    cox_regression <- coxph(Surv(time, status) ~ risk_score, data = x)
    c_index <- cox_regression$concordance[6]
    return(c_index)
  })
  output <- list(RS = RS, C_Index = C_indexs)
  return(output)
}
##### superPC
superPC_function <- function(TrainingDataset, ValidationDatasets, seed) {
  set.seed(seed)
  data <- list(
    x = t(TrainingDataset$EXP),
    y = TrainingDataset$Analysis_table$time,
    censoring.status = TrainingDataset$Analysis_table$status,
    featurenames = colnames(TrainingDataset$EXP)
  )
  fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
  set.seed(seed)
  cv.fit <- superpc.cv(fit, data,
    n.threshold = 20,
    n.fold = 10,
    n.components = 1,
    min.features = 1,
    max.features = nrow(data$x),
    compute.fullcv = TRUE,
    compute.preval = TRUE
  )
  threshold_selected <- cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])]
  variable_selected <- names(which(abs(fit$feature.scores) > threshold_selected))
  RS <- lapply(ValidationDatasets, function(w) {
    newdata <- list(
      x = t(w$EXP),
      y = w$Analysis_table$time,
      censoring.status = w$Analysis_table$status,
      featurenames = colnames(w$EXP)
    )
    set.seed(seed)
    prediction <- superpc.predict(fit, data = data, newdata = newdata, threshold = threshold_selected, n.components = 1)
    risk_score <- as.numeric(prediction$v.pred)
    rs <- cbind(w$Survival, risk_score)
    rs <- as.data.frame(rs)
    return(rs)
  })
  C_indexs <- sapply(RS, function(x) {
    cox_regression <- coxph(Surv(time, status) ~ risk_score, data = x)
    c_index <- cox_regression$concordance[6]
    return(c_index)
  })
  output <- list(RS = RS, C_Index = C_indexs, Selected_Variable = variable_selected)
  return(output)
}
##### Lasso
Lasso_function <- function(TrainingDataset, ValidationDatasets, seed) {
  set.seed(seed)
  fit <- cv.glmnet(x = TrainingDataset$EXP, y = TrainingDataset$Survival, family = "cox", type.measure = "deviance", alpha = 1)
  coef <- as.matrix(coef(fit$glmnet.fit, s = fit$lambda.min))
  variable_selected <- rownames(coef)[coef != 0]
  RS <- lapply(ValidationDatasets, function(x) {
    risk_score <- as.numeric(predict(fit, newx = x$EXP, type = "link", s = fit$lambda.min))
    names(risk_score) <- row.names(x$Analysis_table)
    rs <- cbind(x$Survival, risk_score)
    rs <- as.data.frame(rs)
    return(rs)
  })
  C_indexs <- sapply(RS, function(x) {
    cox_regression <- coxph(Surv(time, status) ~ risk_score, data = x)
    c_index <- cox_regression$concordance[6]
    return(c_index)
  })
  output <- list(RS = RS, C_Index = C_indexs, Selected_Variable = variable_selected)
  return(output)
}
##### Ridge
Ridge_function <- function(TrainingDataset, ValidationDatasets, seed) {
  set.seed(seed)
  fit <- cv.glmnet(x = TrainingDataset$EXP, y = TrainingDataset$Survival, family = "cox", type.measure = "deviance", alpha = 0)
  RS <- lapply(ValidationDatasets, function(x) {
    risk_score <- as.numeric(predict(fit, newx = x$EXP, type = "link", s = fit$lambda.min))
    names(risk_score) <- row.names(x$Analysis_table)
    rs <- cbind(x$Survival, risk_score)
    rs <- as.data.frame(rs)
    return(rs)
  })
  C_indexs <- sapply(RS, function(x) {
    cox_regression <- coxph(Surv(time, status) ~ risk_score, data = x)
    c_index <- cox_regression$concordance[6]
    return(c_index)
  })
  output <- list(RS = RS, C_Index = C_indexs)
  return(output)
}
##### Elastic Net
Elastic_Net_function <- function(TrainingDataset, ValidationDatasets, seed) {
  output <- list()
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)
    fit <- cv.glmnet(x = TrainingDataset$EXP, y = TrainingDataset$Survival, family = "cox", type.measure = "deviance", alpha = alpha)
    coef <- as.matrix(coef(fit$glmnet.fit, s = fit$lambda.min))
    variable_selected <- rownames(coef)[coef != 0]
    RS <- lapply(ValidationDatasets, function(x) {
      risk_score <- as.numeric(predict(fit, newx = x$EXP, type = "link", s = fit$lambda.min))
      names(risk_score) <- row.names(x$Analysis_table)
      rs <- cbind(x$Survival, risk_score)
      rs <- as.data.frame(rs)
      return(rs)
    })
    C_indexs <- sapply(RS, function(x) {
      cox_regression <- coxph(Surv(time, status) ~ risk_score, data = x)
      c_index <- cox_regression$concordance[6]
      return(c_index)
    })
    output[[paste0("alpha=", alpha)]] <- list(RS = RS, C_Index = C_indexs, Selected_Variable = variable_selected)
  }
  return(output)
}
##### Stepwise cox
Stepwise_cox_function <- function(TrainingDataset, ValidationDatasets, seed) {
  output <- list()
  for (direction in c("forward", "backward", "both")) {
    set.seed(seed)
    fit <- step(coxph(Surv(time, status) ~ ., data = TrainingDataset$Analysis_table), direction = direction)
    variable_selected <- names(fit$coefficients)
    RS <- lapply(ValidationDatasets, function(x) {
      risk_score <- as.numeric(predict(fit, type = "lp", newdata = as.data.frame(x$EXP)))
      names(risk_score) <- row.names(x$Analysis_table)
      rs <- cbind(x$Survival, risk_score)
      rs <- as.data.frame(rs)
      return(rs)
    })
    C_indexs <- sapply(RS, function(x) {
      cox_regression <- coxph(Surv(time, status) ~ risk_score, data = x)
      c_index <- cox_regression$concordance[6]
      return(c_index)
    })
    output[[direction]] <- list(RS = RS, C_Index = C_indexs, Selected_Variable = variable_selected)
  }
  return(output)
}


ml_run <- function(Methods, args) {
  # args_example: args <- list(TrainingDataset=training_dataset,ValidationDatasets=all_datasets,seed=seed)
  if (length(Methods) == 1) {
    return(do.call(paste0(Methods, "_function"), args = args))
  } else {
    screen_method <- Methods[1]
    fit_method <- Methods[2]
    if (screen_method %in% c("Elastic_Net", "Stepwise_cox")) {
      screen_result <- do.call(paste0(screen_method, "_function"), args = args)
      output <- list()
      for (i in 1:length(screen_result)) {
        sub_screen_result <- screen_result[[i]]
        screened_variables <- sub_screen_result$Selected_Variable
        TrainingDataset_fit <- args$TrainingDataset
        ValidationDatasets_fit <- args$ValidationDatasets
        TrainingDataset_fit$EXP <- TrainingDataset_fit$EXP[, screened_variables]
        TrainingDataset_fit$Analysis_table <- cbind(TrainingDataset_fit$Analysis_table[, 1:2], TrainingDataset_fit$Analysis_table[, screened_variables])
        ValidationDatasets_fit <- lapply(ValidationDatasets_fit, function(x) {
          EXP <- x$EXP[, screened_variables]
          Analysis_table <- cbind(x$Analysis_table[, 1:2], x$Analysis_table[, screened_variables])
          return_list <- list(EXP = EXP, Analysis_table = Analysis_table, Clinical = x$Clinical, Survival = x$Survival)
          return(return_list)
        })
        args_fit <- list(TrainingDataset = TrainingDataset_fit, ValidationDatasets = ValidationDatasets_fit, seed = args$seed)
        fit_output <- do.call(paste0(fit_method, "_function"), args = args_fit)
        output[[names(screen_result)[i]]] <- fit_output
        print(i)
      }
      return(output)
    } else {
      screen_result <- do.call(paste0(screen_method, "_function"), args = args)
      screened_variables <- screen_result$Selected_Variable
      TrainingDataset_fit <- args$TrainingDataset
      ValidationDatasets_fit <- args$ValidationDatasets
      TrainingDataset_fit$EXP <- TrainingDataset_fit$EXP[, screened_variables]
      TrainingDataset_fit$Analysis_table <- cbind(TrainingDataset_fit$Analysis_table[, 1:2], TrainingDataset_fit$Analysis_table[, screened_variables])
      ValidationDatasets_fit <- lapply(ValidationDatasets_fit, function(x) {
        EXP <- x$EXP[, screened_variables]
        Analysis_table <- cbind(x$Analysis_table[, 1:2], x$Analysis_table[, screened_variables])
        return_list <- list(EXP = EXP, Analysis_table = Analysis_table, Clinical = x$Clinical, Survival = x$Survival)
        return(return_list)
      })
      args_fit <- list(TrainingDataset = TrainingDataset_fit, ValidationDatasets = ValidationDatasets_fit, seed = args$seed)
      output <- do.call(paste0(fit_method, "_function"), args = args_fit)
      return(output)
    }
  }
}
