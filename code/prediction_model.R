#######################################################
### Deep Brain Stimulation: When to go directional.
### Prediction Model
### 2022 Katrin Petermann
#######################################################


#######################################################
### Create new Dataframe with labels for prediction
#######################################################

dd <- data.frame(id = c(MDSTNdata$patient, MDSTNdata$patient, MDSTNdata$patient, MDSTNdata$patient),
                 level = rep(c("best", "worst"), each = 2*dim(MDSTNdata)[1]),
                 SET_level = c(MDSTNdata$`TW_BL_SET_C1-8`, MDSTNdata$`TW_BL_SET_C9-16`, MDSTNdata$`TW_WL_SET_C1-8`, MDSTNdata$`TW_WL_SET_C9-16`),
                 ET_level = c(MDSTNdata$`TW_BL_ET_C1-8`, MDSTNdata$`TW_BL_ET_C9-16`, MDSTNdata$`TW_WL_ET_C1-8`, MDSTNdata$`TW_WL_ET_C9-16`),
                 TW_dir = c(MDSTNdata$`TW_BL-I_C1-8`, MDSTNdata$`TW_BL-I_C9-16`, MDSTNdata$`TW_WL-I_C1-8`, MDSTNdata$`TW_WL-I_C9-16`),
                 TW_level = c(MDSTNdata$`TW_BL_C1-8`, MDSTNdata$`TW_BL_C9-16`, MDSTNdata$`TW_WL_C1-8`, MDSTNdata$ `TW_WL_C9-16`))

dd <- na.omit(dd)

dd$TW_diff <- dd$TW_dir - dd$TW_level
dd$TW_diff_perc <- 100 * (dd$TW_diff)/(dd$TW_level)

## label to predict --> im TW is at least 25% larger
dd$TW_larger <- ifelse(dd$TW_diff_perc > 25, 1, 0)
dd$TW_larger[is.na(dd$TW_larger)] <- 0
dd$TW_larger <- factor(dd$TW_larger, levels = c(0, 1), labels = c("no", "yes"))

dd <- dd[order(dd$TW_level), ]

## Transform dataframe to long format

dd_long <- rbind(as.matrix(cbind("TW", dd[, c("id", "TW_level", "TW_larger")])),
                 as.matrix(cbind("ET", dd[, c("id", "ET_level", "TW_larger")])),
                 as.matrix(cbind("SET", dd[, c("id", "SET_level", "TW_larger")])))
dd_long <- as.data.frame(dd_long)
colnames(dd_long) = c("measure", "pid", "amplitude", "TW_larger")
dd_long$measure <- as.factor(dd_long$measure)
dd_long$amplitude <- as.numeric(dd_long$amplitude)
dd_long$TW_larger <- as.factor(dd_long$TW_larger)

## Figure 2 --> compiled in different file figures.R

#######################################################
### ROC Models to test different predictors
#######################################################

## Combined predictor TW and ET 
dd$TW_ET_level <- dd$TW_level - dd$ET_level

roc_ET <- roc(TW_larger ~ ET_level, data = dd, direction = "<", quiet = TRUE)
roc_SET <- roc(TW_larger ~ SET_level, data = dd, direction = ">", quiet = TRUE)
roc_TW <- roc(TW_larger ~ TW_level, data = dd, direction = ">", quiet = TRUE)
roc_TW_ET <- roc(TW_larger ~ TW_ET_level, data = dd, direction = ">", quiet = TRUE)

## calculate confidence intervalls for ROC curves
roc_ci <- lapply(list(roc_ET, roc_SET, roc_TW, roc_TW_ET), ci.se, specificities = seq(0, 1, l = 25))


dd_roc <- lapply(roc_ci, function(x){
  data.frame(x = as.numeric(rownames(x)),
             lower = x[, 1],
             upper = x[, 3])
})

## Figure 3a --> compiled in different file figures.R

## statistical difference of roc curves

## ET vs. SET
roc_test_ET_SET <- roc.test(roc_SET, roc_ET, 
                            paired = TRUE, method = "bootstrap")
## TW vs. SET
roc_test_TW_SET <- roc.test(roc_SET, roc_TW, 
                            paired = TRUE, method = "bootstrap")

## TW vs. ET
roc_test_TW_ET <- roc.test(roc_TW, roc_ET, 
                           paired = TRUE, method = "bootstrap")

## TW vs. ET
roc_test_TW_ET_TW <- roc.test(roc_TW, roc_TW_ET, 
                              paired = TRUE, method = "bootstrap")

#######################################################
### Crossvalidation for Predictive performance
#######################################################

n_samples <- dim(dd)[1]
n_folds <- 5
sens_thresh <- 0.75 ## threshold sensitivity

results <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(results) <- c("fold", "TW_threshold", "ET_threshold", "SET_threshold")


## create random folds with equal distribution of labels 0/1
dd_POS <- dd[dd$TW_larger == "yes", ]
dd_NEG <- dd[dd$TW_larger == "no", ]

set.seed(2808)
index_POS <- sample(cut(1:dim(dd_POS)[1], breaks = n_folds, labels = FALSE))
index_NEG <- sample(cut(1:dim(dd_NEG)[1], breaks = n_folds, labels = FALSE))

dd$fold_index[dd$TW_larger == "yes"] <- index_POS
dd$fold_index[dd$TW_larger == "no"] <- index_NEG

for(fold in 1:n_folds){
  ## define test and training data
  dd_train <- dd[dd$fold_index != fold, ]
  dd_test <- dd[dd$fold_index == fold, ]
  
  ## ROC Model
  roc_TW <- roc(TW_larger ~ TW_level, data = dd_train, direction = ">", quiet = TRUE)
  roc_ET <- roc(TW_larger ~ ET_level, data = dd_train, direction = "<", quiet = TRUE)
  roc_SET <- roc(TW_larger ~ SET_level, data = dd_train, direction = ">", quiet = TRUE)

  ## get thresholds for sensitivity
  tab_TW <- coords(roc_TW)
  tab_ET <- coords(roc_ET)
  tab_SET <- coords(roc_SET)
  
  th_TW <- min(tab_TW$threshold[tab_TW$sensitivity >= sens_thresh])
  th_ET <- max(tab_ET$threshold[tab_ET$sensitivity >= sens_thresh])
  th_SET <- min(tab_SET$threshold[tab_SET$sensitivity >= sens_thresh])
  
  ## predict if the contact will be retested (below or above threshold)
  dd$pred_TW[dd$fold_index == fold] <- ifelse(dd$TW_level[dd$fold_index == fold] <= th_TW, "yes", "no")
  dd$pred_ET[dd$fold_index == fold] <- ifelse(dd$ET_level[dd$fold_index == fold] >= th_ET, "yes", "no")
  dd$pred_SET[dd$fold_index == fold] <- ifelse(dd$SET_level[dd$fold_index == fold] <= th_SET, "yes", "no")
  dd$pred_TW_ET[dd$fold_index == fold] <- ifelse(dd$pred_TW[dd$fold_index == fold] == "yes" & dd$pred_ET[dd$fold_index == fold] == "yes", "yes", "no")
  
  thresholds <- rbind(results, data.frame(fold = fold, 
                                       TW_threshold = th_TW, 
                                       ET_threshold = th_ET, 
                                       SET_threshold = th_SET))
}

## calculate CV sensitivity and specificity for predictions
P <- sum(dd$TW_larger == "yes")
N <- sum(dd$TW_larger == "no")
TP_TW <- sum((dd$TW_larger == dd$pred_TW) & (dd$TW_larger == "yes"))
TN_TW <- sum((dd$TW_larger == dd$pred_TW) & (dd$TW_larger == "no"))

TP_ET <- sum((dd$TW_larger == dd$pred_ET) & (dd$TW_larger == "yes"))
TN_ET <- sum((dd$TW_larger == dd$pred_ET) & (dd$TW_larger == "no"))

TP_SET <- sum((dd$TW_larger == dd$pred_SET) & (dd$TW_larger == "yes"))
TN_SET <- sum((dd$TW_larger == dd$pred_SET) & (dd$TW_larger == "no"))

TP_TW_ET <- sum((dd$TW_larger == dd$pred_TW_ET) & (dd$TW_larger == "yes"))
TN_TW_ET <- sum((dd$TW_larger == dd$pred_TW_ET) & (dd$TW_larger == "no"))

sens_CV_TW <- TP_TW/P
spec_CV_TW <- TN_TW/N
acc_CV_TW <- (TP_TW + TN_TW) / n_samples

sens_CV_ET <- TP_ET/P
spec_CV_ET <- TN_ET/N
acc_CV_ET <- (TP_ET + TN_ET) / n_samples

sens_CV_SET <- TP_SET/P
spec_CV_SET <- TN_SET/N
acc_CV_SET <- (TP_SET + TN_SET) / n_samples

sens_CV_TW_ET <- TP_TW_ET/P
spec_CV_TW_ET <- TN_TW_ET/N
acc_CV_TW_ET <- (TP_TW_ET + TN_TW_ET) / n_samples

get_ciP <- function(measure, result, nom, denom){
  x <- confIntProportion(nom, denom)
  return(data.frame(measure = measure, result = result, estimate = x$p, lower = x$CIs[2, 2], upper = x$CIs[2, 3]))
}

## dataframe with prediction results
results_CV <- rbind(get_ciP("TW", "Sensitivity", TP_TW, P),
                    get_ciP("TW", "Specificity", TN_TW, N),
                    get_ciP("TW", "Accuracy", TP_TW + TN_TW, n_samples),
                    get_ciP("ET", "Sensitivity", TP_ET, P),
                    get_ciP("ET", "Specificity", TN_ET, N),
                    get_ciP("ET", "Accuracy", TP_ET + TN_ET, n_samples),
                    get_ciP("SET", "Sensitivity", TP_SET, P),
                    get_ciP("SET", "Specificity", TN_SET, N),
                    get_ciP("SET", "Accuracy", TP_SET + TN_SET, n_samples),
                    get_ciP("TW ET", "Sensitivity", TP_TW_ET, P),
                    get_ciP("TW ET", "Specificity", TN_TW_ET, N),
                    get_ciP("TW ET", "Accuracy", TP_TW_ET + TN_TW_ET, n_samples))


## Figure 3b --> compiled in different file figures.R

## CV for TW_threshold 1.75 and ET_threshold 2.25
TW_th <- 1.5

dd$pred <- ifelse(dd$TW_level <= TW_th, "yes", "no")

TP <- sum((dd$TW_larger == dd$pred) & (dd$TW_larger == "yes"))
TN <- sum((dd$TW_larger == dd$pred) & (dd$TW_larger == "no"))

sens <- confIntProportion(TP, P)
spec <- confIntProportion(TN, N)
acc <- confIntProportion(TP + TN, n_samples)

tab_pred <- data.frame(pred = c("Therapeutic Window", "Effect Threshold", 
                                "Side Effect Threshold", "TW and ET"),
                       sens = c(results_CV[1, 3], results_CV[4, 3], results_CV[7, 3], results_CV[10, 3]),
                       ci_sens = c(formatCI(results_CV[1, 4:5]), 
                                   formatCI(results_CV[4, 4:5]), 
                                   formatCI(results_CV[7, 4:5]), 
                                   formatCI(results_CV[10, 4:5])),
                       spec = c(results_CV[2, 3], results_CV[5, 3], results_CV[8, 3], results_CV[11, 3]),
                       ci_spec = c(formatCI(results_CV[2, 4:5]), 
                                   formatCI(results_CV[5, 4:5]), 
                                   formatCI(results_CV[8, 4:5]), 
                                   formatCI(results_CV[11, 4:5])),
                       acc = c(results_CV[3, 3], results_CV[6, 3], results_CV[9, 3], results_CV[12, 3]),
                       ci_acc = c(formatCI(results_CV[3, 4:5]), 
                                  formatCI(results_CV[6, 4:5]), 
                                  formatCI(results_CV[9, 4:5]), 
                                  formatCI(results_CV[12, 4:5])))