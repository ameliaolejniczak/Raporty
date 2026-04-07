library(pracma)
library(glmnet)
library(ggplot2)
library(knitr)
library(dplyr)

set.seed(47282456)
#standaryzacja ma znaczenie, gdy lambda jest ustalona
#napisać co i jak optymalizujemy (przy lambda i knockoffach)

#zadanie 1.

X <- randortho(500, "orthonormal")
k <- c(5, 20, 100)

beta <- setNames(lapply(k, function(k) {
  beta1 <- c(rep(4, k), rep(0, 500 - k))
  beta2 <- c(rep(sqrt(80/k), k), rep(0, 500 - k))
  return(list(beta1 = beta1, beta2 = beta2))
}), paste0("k", k))

Y <- lapply(beta, function(b) {
  Y_beta1 <- X%*%b[["beta1"]] + rnorm(500, 0, 1)
  Y_beta2 <- X%*%b[["beta2"]] + rnorm(500, 0, 1)
  return(list(Y_beta1 = Y_beta1, Y_beta2 = Y_beta2))
})

gamma <- seq(0.01, 50, 0.01)

data1 <- setNames(lapply(names(Y), function(name_k) {
  y <- Y[[name_k]]
  
  beta_ridge_beta1 <- matrix(rep(0, 500*length(gamma)), nrow = 500)
  MSE_ridge_beta1 <- rep(0, length(gamma))
  beta_ridge_beta2 <- matrix(rep(0, 500*length(gamma)), nrow = 500)
  MSE_ridge_beta2 <- rep(0, length(gamma))
  
  for (i in 1:length(gamma)) {
    beta_ridge_beta1[,i] <- t(X)%*%y[["Y_beta1"]]/(1 + gamma[i])
    MSE_ridge_beta1[i] <- sum((beta[[name_k]][["beta1"]]-beta_ridge_beta1[,i])^2)
    
    beta_ridge_beta2[,i] <- t(X)%*%y[["Y_beta2"]]/(1 + gamma[i])
    MSE_ridge_beta2[i] <- sum((beta[[name_k]][["beta2"]]-beta_ridge_beta2[,i])^2)
  }
  
  gamma_opt_beta1 <- 500/sum(beta[[name_k]][["beta1"]]^2)
  beta_ridge_beta1_optimal <- t(X)%*%y[["Y_beta1"]]/(1 + gamma_opt_beta1)
  MSE_ridge_beta1_optimal <- sum((beta[[name_k]][["beta1"]]-beta_ridge_beta1_optimal)^2)
  
  gamma_opt_beta2 <- 500/sum(beta[[name_k]][["beta2"]]^2)
  beta_ridge_beta2_optimal <- t(X)%*%y[["Y_beta2"]]/(1 + gamma_opt_beta2)
  MSE_ridge_beta2_optimal <- sum((beta[[name_k]][["beta2"]]-beta_ridge_beta2_optimal)^2)
  
  lasso_fits_beta1 <- glmnet(X, y[["Y_beta1"]], intercept = F, alpha = 1) 
  fitted_lambda_beta1 <- lasso_fits_beta1$lambda
  MSE_lasso_beta1 <- apply(coef(lasso_fits_beta1)[-1, ], 2, function(beta_l) {
    sum((beta_l - beta[[name_k]][["beta1"]])^2)
  })
  
  index_lasso_beta1 <- which.min(MSE_lasso_beta1)
  lambda_beta1 <- fitted_lambda_beta1[index_lasso_beta1]
  MSE_lasso_beta1_optimal <- MSE_lasso_beta1[index_lasso_beta1]
  
  lasso_fits_beta2 <- glmnet(X, y[["Y_beta2"]], intercept = F, alpha = 1) 
  fitted_lambda_beta2 <- lasso_fits_beta2$lambda
  MSE_lasso_beta2 <- apply(coef(lasso_fits_beta2)[-1, ], 2, function(beta_l) {
    sum((beta_l - beta[[name_k]][["beta2"]])^2)
  })
  
  index_lasso_beta2 <- which.min(MSE_lasso_beta2)
  lambda_beta2 <- fitted_lambda_beta2[index_lasso_beta2]
  MSE_lasso_beta2_optimal <- MSE_lasso_beta2[index_lasso_beta2]
  
  beta_ls_beta1 <- t(X) %*% y[["Y_beta1"]]
  MSE_ls_beta1 <- sum((beta[[name_k]][["beta1"]] - beta_ls_beta1)^2)
  
  beta_ls_beta2 <- t(X) %*% y[["Y_beta2"]]
  MSE_ls_beta2 <- sum((beta[[name_k]][["beta2"]] - beta_ls_beta2)^2)
  
  return(list(
    MSE_ridge_beta1 = MSE_ridge_beta1, 
    gamma_opt_beta1 = gamma_opt_beta1, 
    MSE_ridge_beta1_optimal = MSE_ridge_beta1_optimal,
    beta_ridge_beta1 = beta_ridge_beta1_optimal,
    
    MSE_ridge_beta2 = MSE_ridge_beta2, 
    gamma_opt_beta2 = gamma_opt_beta2, 
    MSE_ridge_beta2_optimal = MSE_ridge_beta2_optimal,
    beta_ridge_beta2 = beta_ridge_beta2_optimal,
    
    fitted_lambda_beta1 = fitted_lambda_beta1,
    MSE_lasso_beta1 = MSE_lasso_beta1, 
    lambda_beta1 = lambda_beta1, 
    MSE_lasso_beta1_optimal = MSE_lasso_beta1_optimal,
    
    fitted_lambda_beta2 = fitted_lambda_beta2,
    MSE_lasso_beta2 = MSE_lasso_beta2, 
    lambda_beta2 = lambda_beta2, 
    MSE_lasso_beta2_optimal = MSE_lasso_beta2_optimal,
    
    MSE_ls_beta1 = MSE_ls_beta1,
    MSE_ls_beta2 = MSE_ls_beta2)
  )
}), names(beta))

setwd("C:\\Users\\ameli\\Documents\\R\\analiza duzych zbiorow danych")
saveRDS(data1, "data1.rds")

MSE_ridge <- data.frame(
  MSE_ridge_beta1 <- c(data1[["k5"]][["MSE_ridge_beta1_optimal"]], data1[["k20"]][["MSE_ridge_beta1_optimal"]], data1[["k100"]][["MSE_ridge_beta1_optimal"]]),
  MSE_ridge_beta2 <- c(data1[["k5"]][["MSE_ridge_beta2_optimal"]], data1[["k20"]][["MSE_ridge_beta2_optimal"]], data1[["k100"]][["MSE_ridge_beta2_optimal"]])
)
colnames(MSE_ridge) <- c("i", "ii")
rownames(MSE_ridge) <- c("k = 5", "k = 20", "k = 100")

kable(MSE_ridge)

MSE_lasso <- data.frame(
  MSE_lasso_beta1 <- c(data1[["k5"]][["MSE_lasso_beta1_optimal"]], data1[["k20"]][["MSE_lasso_beta1_optimal"]], data1[["k100"]][["MSE_lasso_beta1_optimal"]]),
  MSE_lasso_beta2 <- c(data1[["k5"]][["MSE_lasso_beta2_optimal"]], data1[["k20"]][["MSE_lasso_beta2_optimal"]], data1[["k100"]][["MSE_lasso_beta2_optimal"]])
)
colnames(MSE_lasso) <- c("i", "ii")
rownames(MSE_lasso) <- c("k = 5", "k = 20", "k = 100")

kable(MSE_lasso)

MSE_ls <- data.frame(
  MSE_ls_beta1 <- c(data1[["k5"]][["MSE_ls_beta1"]], data1[["k20"]][["MSE_ls_beta1"]], data1[["k100"]][["MSE_ls_beta1"]]),
  MSE_ls_beta2 <- c(data1[["k5"]][["MSE_ls_beta2"]], data1[["k20"]][["MSE_ls_beta2"]], data1[["k100"]][["MSE_ls_beta2"]])
)
colnames(MSE_ls) <- c("i", "ii")
rownames(MSE_ls) <- c("k = 5", "k = 20", "k = 100")

kable(MSE_ls)

fd_power <- lapply(k, function(k) {
  name <- paste0("k", k)
  beta1 <- beta[[name]]$beta1
  beta2 <- beta[[name]]$beta2
  
  lambda1 <- data1[[name]][["lambda_beta1"]]
  lambda2 <- data1[[name]][["lambda_beta2"]]
  
  fd_lasso_beta1 <- (500 - k) * 2 * (1 - pnorm(lambda1))
  fd_lasso_beta2 <- (500 - k) * 2 * (1 - pnorm(lambda2))
  
  power_beta1 <- (1 - pnorm(lambda1 - beta1[1])) + pnorm(- lambda1 - beta1[1])
  power_beta2 <- (1 - pnorm(lambda1 - beta2[1])) + pnorm(- lambda1 - beta2[1])
  
  return(list(
    FD_lasso_beta1 = fd_lasso_beta1,
    FD_lasso_beta2 = fd_lasso_beta2,
    Power_lasso_beta1 = power_beta1,
    Power_lasso_beta2 = power_beta2
  ))
})

saveRDS(fd_power, "fd_power.rds")

#zadanie 2.
X <- matrix(rnorm(500 * 500, 0, 1/sqrt(500)), nrow = 500)

Y <- lapply(beta, function(b) {
  Y_beta1 <- X%*%b[["beta1"]] + rnorm(500, 0, 1)
  Y_beta2 <- X%*%b[["beta2"]] + rnorm(500, 0, 1)
  return(list(Y_beta1 = Y_beta1, Y_beta2 = Y_beta2))
})

data2 <- setNames(lapply(names(Y), function(name_k) {
  y <- Y[[name_k]]
  
  ridge_fit_2_beta1 <- cv.glmnet(X, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 0)
  hbeta_ridge_2_beta1 <- coef(ridge_fit_2_beta1, s = 'lambda.min')[2:501]
  ridge_fit_2_beta2 <- cv.glmnet(X, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 0)
  hbeta_ridge_2_beta2 <- coef(ridge_fit_2_beta2, s = 'lambda.min')[2:501]
  
  lasso_fit_2_beta1 <- cv.glmnet(X, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbeta_lasso_2_beta1 <- coef(lasso_fit_2_beta1, s = 'lambda.min')[2:501]
  lasso_fit_2_beta2 <- cv.glmnet(X, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbeta_lasso_2_beta2 <- coef(lasso_fit_2_beta2, s = 'lambda.min')[2:501]
  
  W1_beta1 <- 1/abs(hbeta_ridge_2_beta1)
  X_temp_beta1 <- X %*% diag(1/W1_beta1)
  adlasso1_fit_beta1 <- cv.glmnet(X_temp_beta1, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbetatemp_adlasso1_beta1 <- coef(adlasso1_fit_beta1, s = 'lambda.min')[2:501]
  hbeta_adlasso1_beta1 <- hbetatemp_adlasso1_beta1/W1_beta1
  W1_beta2 <- 1/(abs(hbeta_ridge_2_beta2))
  X_temp_beta2 <- X %*% diag(1/W1_beta2)
  adlasso1_fit_beta2 <- cv.glmnet(X_temp_beta2, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbetatemp_adlasso1_beta2 <- coef(adlasso1_fit_beta2, s = 'lambda.min')[2:501]
  hbeta_adlasso1_beta2 <- hbetatemp_adlasso1_beta2/W1_beta2
  
  RSS_beta1 <- sum((y[["Y_beta1"]] - X %*% hbeta_lasso_2_beta1)^2)
  ind_lasso_beta1 <- which(hbeta_lasso_2_beta1 > 0)
  if (length(ind_lasso_beta1) == 0) {
    hbeta_adlasso2_beta1 = rep(0, 500)
  } else {
    hsigma_beta1 <- sqrt(RSS_beta1/(500 - length(ind_lasso_beta1)))
    W2_beta1 <- hsigma_beta1/abs(hbeta_lasso_2_beta1[ind_lasso_beta1])
    lambda1_beta1 <- hsigma_beta1 * qnorm(1 - 0.2/(2 * 500))
    X_small_beta1 <- X[, ind_lasso_beta1]
    X_temp_beta1 <- X_small_beta1 %*% diag(1/W2_beta1)
    
    adlasso2_fit_beta1 <- glmnet(X_temp_beta1, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1, lambda = lambda1_beta1/500)
    hbetatemp_adlasso2_beta1 <- coef(adlasso2_fit_beta1)[2:(length(ind_lasso_beta1) +1)]
    hbeta_adlasso2_beta1 <- rep(0, 500)
    hbeta_adlasso2_beta1[ind_lasso_beta1] <- hbetatemp_adlasso2_beta1/W2_beta1
  }
  RSS_beta2 <- sum((y[["Y_beta2"]] - X %*% hbeta_lasso_2_beta2)^2)
  ind_lasso_beta2 <- which(hbeta_lasso_2_beta2 > 0)
  if (length(ind_lasso_beta2) == 0) {
    hbeta_adlasso2_beta2 = rep(0, 500)
  } else {
    hsigma_beta2 <- sqrt(RSS_beta2/(500 - length(ind_lasso_beta2)))
    W2_beta2 <- hsigma_beta2/abs(hbeta_lasso_2_beta2[ind_lasso_beta2])
    lambda1_beta2 <- hsigma_beta2 * qnorm(1 - 0.2/(2 * 500))
    X_small_beta2 <- X[, ind_lasso_beta2]
    X_temp_beta2 <- X_small_beta2 %*% diag(1/W2_beta2)
    
    adlasso2_fit_beta2 <- glmnet(X_temp_beta2, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1, lambda = lambda1_beta2/500)
    hbetatemp_adlasso2_beta2 <- coef(adlasso2_fit_beta2)[2:(length(ind_lasso_beta2) +1)]
    hbeta_adlasso2_beta2 <- rep(0, 500)
    hbeta_adlasso2_beta2[ind_lasso_beta2] <- hbetatemp_adlasso2_beta2/W2_beta2
  }
  
  X_copy <- matrix(rnorm(500 * 500, 0, 1/sqrt(500)), nrow = 500)
  X_large <- cbind(X, X_copy)
  large_fit_ridge_beta1 <- cv.glmnet(X_large, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 0)
  hbeta_large_ridge_beta1 <- coef(large_fit_ridge_beta1, s='lambda.min')[2:(2*500 + 1)]
  w <- abs(hbeta_large_ridge_beta1[1:500]) - abs(hbeta_large_ridge_beta1[501:1000])
  
  kn <- function(t) {
    ((1 + sum(w <= -t))/max(sum(w >= t), 1))
  }
  
  w_kn <- sapply(c(0, abs(w)), kn)
  thres <- min(c(0, abs(w))[w_kn <= 0.2])
  hbeta_knridge_beta1 <- hbeta_ridge_2_beta1 * (w >= thres)
  
  large_fit_ridge_beta2 <- cv.glmnet(X_large, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 0)
  hbeta_large_ridge_beta2 <- coef(large_fit_ridge_beta2, s='lambda.min')[2:(2*500 + 1)]
  w <- abs(hbeta_large_ridge_beta2[1:500]) - abs(hbeta_large_ridge_beta2[501:1000])
  w_kn <- sapply(c(0, abs(w)), kn)
  thres <- min(c(0, abs(w))[w_kn <= 0.2])
  hbeta_knridge_beta2 <- hbeta_ridge_2_beta2 * (w >= thres)
  
  large_fit_lasso_beta1 <- cv.glmnet(X_large, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbeta_large_lasso_beta1 <- coef(large_fit_lasso_beta1, s = 'lambda.min')[2:1001]
  w <- abs(hbeta_large_lasso_beta1[1:500]) - abs(hbeta_large_lasso_beta1[501:1000])
  w_kn_lasso <- sapply(c(0, abs(w)), kn)
  thres_lasso <- min(c(0, abs(w))[w_kn <= 0.2])
  hbeta_knlasso_beta1 <- hbeta_lasso_2_beta1 * (w >= thres)
  large_fit_lasso_beta2 <- cv.glmnet(X_large, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbeta_large_lasso_beta2 <- coef(large_fit_lasso_beta2, s = 'lambda.min')[2:1001]
  w <- abs(hbeta_large_lasso_beta2[1:500]) - abs(hbeta_large_lasso_beta2[501:1000])
  w_kn_lasso <- sapply(c(0, abs(w)), kn)
  thres_lasso <- min(c(0, abs(w))[w_kn <= 0.2])
  hbeta_knlasso_beta2 <- hbeta_lasso_2_beta2 * (w >= thres)
  
  SE_ridge_beta1 <- sum((hbeta_ridge_2_beta1 - beta[[name_k]][["beta1"]])^2)
  SE_ridge_beta2 <- sum((hbeta_ridge_2_beta2 - beta[[name_k]][["beta2"]])^2)
  SE_lasso_beta1 <- sum((hbeta_lasso_2_beta1 - beta[[name_k]][["beta1"]])^2)
  SE_lasso_beta2 <- sum((hbeta_lasso_2_beta2 - beta[[name_k]][["beta2"]])^2)
  SE_adlasso1_beta1 <- sum((hbeta_adlasso1_beta1 - beta[[name_k]][["beta1"]])^2)
  SE_adlasso1_beta2 <- sum((hbeta_adlasso1_beta2 - beta[[name_k]][["beta2"]])^2)
  SE_adlasso2_beta1 <- sum((hbeta_adlasso2_beta1 - beta[[name_k]][["beta1"]])^2)
  SE_adlasso2_beta2 <- sum((hbeta_adlasso2_beta2 - beta[[name_k]][["beta2"]])^2)
  
  true_beta <- beta[[name_k]][["beta1"]] != 0
  k <- sum(true_beta)
  
  lasso_selected_beta1 <- hbeta_lasso_2_beta1 != 0
  TP_lasso_beta1 <- sum(true_beta & lasso_selected_beta1)
  FD_lasso_beta1 <- sum(!true_beta & lasso_selected_beta1)
  R_lasso_beta1 <- sum(lasso_selected_beta1)
  power_lasso_beta1 <- TP_lasso_beta1 / k
  FDR_lasso_beta1 <- if (R_lasso_beta1 == 0) 0 else FD_lasso_beta1 / R_lasso_beta1
  lasso_selected_beta2 <- hbeta_lasso_2_beta2 != 0
  TP_lasso_beta2 <- sum(true_beta & lasso_selected_beta2)
  FD_lasso_beta2 <- sum(!true_beta & lasso_selected_beta2)
  R_lasso_beta2 <- sum(lasso_selected_beta2)
  power_lasso_beta2 <- TP_lasso_beta2 / k
  FDR_lasso_beta2 <- if (R_lasso_beta2 == 0) 0 else FD_lasso_beta2 / R_lasso_beta2
  
  adlasso1_selected_beta1 <- hbeta_adlasso1_beta1 != 0
  TP_adlasso1_beta1 <- sum(true_beta & adlasso1_selected_beta1)
  FD_adlasso1_beta1 <- sum(!true_beta & adlasso1_selected_beta1)
  R_adlasso1_beta1 <- sum(adlasso1_selected_beta1)
  power_adlasso1_beta1 <- TP_adlasso1_beta1 / k
  FDR_adlasso1_beta1 <- if (R_adlasso1_beta1 == 0) 0 else FD_adlasso1_beta1 / R_adlasso1_beta1
  adlasso1_selected_beta2 <- hbeta_adlasso1_beta2 != 0
  TP_adlasso1_beta2 <- sum(true_beta & adlasso1_selected_beta2)
  FD_adlasso1_beta2 <- sum(!true_beta & adlasso1_selected_beta2)
  R_adlasso1_beta2 <- sum(adlasso1_selected_beta2)
  power_adlasso1_beta2 <- TP_adlasso1_beta2 / k
  FDR_adlasso1_beta2 <- if (R_adlasso1_beta2 == 0) 0 else FD_adlasso1_beta2 / R_adlasso1_beta2
  
  adlasso2_selected_beta1 <- hbeta_adlasso2_beta1 != 0
  TP_adlasso2_beta1 <- sum(true_beta & adlasso2_selected_beta1)
  FD_adlasso2_beta1 <- sum(!true_beta & adlasso2_selected_beta1)
  R_adlasso2_beta1 <- sum(adlasso2_selected_beta1)
  power_adlasso2_beta1 <- TP_adlasso2_beta1 / k
  FDR_adlasso2_beta1 <- if (R_adlasso2_beta1 == 0) 0 else FD_adlasso2_beta1 / R_adlasso2_beta1
  adlasso2_selected_beta2 <- hbeta_adlasso2_beta2 != 0
  TP_adlasso2_beta2 <- sum(true_beta & adlasso2_selected_beta2)
  FD_adlasso2_beta2 <- sum(!true_beta & adlasso2_selected_beta2)
  R_adlasso2_beta2 <- sum(adlasso2_selected_beta2)
  power_adlasso2_beta2 <- TP_adlasso2_beta2 / k
  FDR_adlasso2_beta2 <- if (R_adlasso2_beta2 == 0) 0 else FD_adlasso2_beta2 / R_adlasso2_beta2
  
  ridgekn_selected_beta1 <- hbeta_knridge_beta1 != 0
  TP_ridgekn_beta1 <- sum(true_beta & ridgekn_selected_beta1)
  FD_ridgekn_beta1 <- sum(!true_beta & ridgekn_selected_beta1)
  R_ridgekn_beta1 <- sum(ridgekn_selected_beta1)
  power_ridgekn_beta1 <- TP_ridgekn_beta1 / k
  FDR_ridgekn_beta1 <- if (R_ridgekn_beta1 == 0) 0 else FD_ridgekn_beta1 / R_ridgekn_beta1
  ridgekn_selected_beta2 <- hbeta_knridge_beta2 != 0
  TP_ridgekn_beta2 <- sum(true_beta & ridgekn_selected_beta2)
  FD_ridgekn_beta2 <- sum(!true_beta & ridgekn_selected_beta2)
  R_ridgekn_beta2 <- sum(ridgekn_selected_beta2)
  power_ridgekn_beta2 <- TP_ridgekn_beta2 / k
  FDR_ridgekn_beta2 <- if (R_ridgekn_beta2 == 0) 0 else FD_ridgekn_beta2 / R_ridgekn_beta2
  
  lassokn_selected_beta1 <- hbeta_knlasso_beta1 != 0
  TP_lassokn_beta1 <- sum(true_beta & lassokn_selected_beta1)
  FD_lassokn_beta1 <- sum(!true_beta & lassokn_selected_beta1)
  R_lassokn_beta1 <- sum(lassokn_selected_beta1)
  power_lassokn_beta1 <- TP_lassokn_beta1 / k
  FDR_lassokn_beta1 <- if (R_lassokn_beta1 == 0) 0 else FD_lassokn_beta1 / R_lassokn_beta1
  lassokn_selected_beta2 <- hbeta_knlasso_beta2 != 0
  TP_lassokn_beta2 <- sum(true_beta & lassokn_selected_beta2)
  FD_lassokn_beta2 <- sum(!true_beta & lassokn_selected_beta2)
  R_lassokn_beta2 <- sum(lassokn_selected_beta2)
  power_lassokn_beta2 <- TP_lassokn_beta2 / k
  FDR_lassokn_beta2 <- if (R_lassokn_beta2 == 0) 0 else FD_lassokn_beta2 / R_lassokn_beta2
  
  return(list(
    hbeta_ridge_beta1 = hbeta_ridge_2_beta1,
    SE_ridge_beta1 = SE_ridge_beta1,
    hbeta_ridge_beta2 = hbeta_ridge_2_beta2,
    SE_ridge_beta2 = SE_ridge_beta2,
    hbeta_lasso_beta1 = hbeta_lasso_2_beta1,
    SE_lasso_beta1 = SE_lasso_beta1,
    power_lasso_beta1 = power_lasso_beta1,
    FDR_lasso_beta1 = FDR_lasso_beta1,
    hbeta_lasso_beta2 = hbeta_lasso_2_beta2,
    SE_lasso_beta2 = SE_lasso_beta2,
    power_lasso_beta2 = power_lasso_beta2,
    FDR_lasso_beta2 = FDR_lasso_beta2,
    hbeta_adlasso1_beta1 = hbeta_adlasso1_beta1,
    SE_adlasso1_beta1 = SE_adlasso1_beta1,
    power_adlasso1_beta1 = power_adlasso1_beta1,
    FDR_adlasso1_beta1 = FDR_adlasso1_beta1,
    hbeta_adlasso1_beta2 = hbeta_adlasso1_beta2,
    SE_adlasso1_beta2 = SE_adlasso1_beta2,
    power_adlasso1_beta2 = power_adlasso1_beta2,
    FDR_adlasso1_beta2 = FDR_adlasso1_beta2,
    hbeta_adlasso2_beta1 = hbeta_adlasso2_beta1,
    SE_adlasso2_beta1 = SE_adlasso2_beta1,
    power_adlasso2_beta1 = power_adlasso2_beta1,
    FDR_adlasso2_beta1 = FDR_adlasso2_beta1,
    hbeta_adlasso2_beta2 = hbeta_adlasso2_beta2,
    SE_adlasso2_beta2 = SE_adlasso2_beta2,
    power_adlasso2_beta2 = power_adlasso2_beta2,
    FDR_adlasso2_beta2 = FDR_adlasso2_beta2,
    hbeta_knridge_beta1 = hbeta_knridge_beta1,
    power_ridgekn_beta1 = power_ridgekn_beta1,
    FDR_ridgekn_beta1 = FDR_ridgekn_beta1,
    hbeta_knridge_beta2 = hbeta_knridge_beta2,
    power_ridgekn_beta2 = power_ridgekn_beta2,
    FDR_ridgekn_beta2 = FDR_ridgekn_beta2,
    hbeta_knlasso_beta1 = hbeta_knlasso_beta1,
    power_lassokn_beta1 = power_lassokn_beta1,
    FDR_lassokn_beta1 = FDR_lassokn_beta1,
    hbeta_knlasso_beta2 = hbeta_knlasso_beta2,
    power_lassokn_beta2 = power_lassokn_beta2,
    FDR_lassokn_beta2 = FDR_lassokn_beta2)
  )
}), names(beta))

saveRDS(data2, "data2.rds")

SE_beta1 <- data.frame(
  SE_ridge_beta1 <- c(data2[["k5"]][["SE_ridge_beta1"]], data2[["k20"]][["SE_ridge_beta1"]], data2[["k100"]][["SE_ridge_beta1"]]),
  SE_lasso_beta1 <- c(data2[["k5"]][["SE_lasso_beta1"]], data2[["k20"]][["SE_lasso_beta1"]], data2[["k100"]][["SE_lasso_beta1"]]),
  SE_adlasso1_beta1 <- c(data2[["k5"]][["SE_adlasso1_beta1"]], data2[["k20"]][["SE_adlasso1_beta1"]], data2[["k100"]][["SE_adlasso1_beta1"]]),
  SE_adlasso2_beta1 <- c(data2[["k5"]][["SE_adlasso2_beta1"]], data2[["k20"]][["SE_adlasso2_beta1"]], data2[["k100"]][["SE_adlasso2_beta1"]])
)
colnames(SE_beta1) <- c("Ridge", "LASSO", "AdLASSO1", "AdLASSO2")
rownames(SE_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(SE_beta1)

SE_beta2<- data.frame(
  SE_ridge_beta2 <- c(data2[["k5"]][["SE_ridge_beta2"]], data2[["k20"]][["SE_ridge_beta2"]], data2[["k100"]][["SE_ridge_beta2"]]),
  SE_lasso_beta2 <- c(data2[["k5"]][["SE_lasso_beta2"]], data2[["k20"]][["SE_lasso_beta2"]], data2[["k100"]][["SE_lasso_beta2"]]),
  SE_adlasso1_beta2 <- c(data2[["k5"]][["SE_adlasso1_beta2"]], data2[["k20"]][["SE_adlasso1_beta2"]], data2[["k100"]][["SE_adlasso1_beta2"]]),
  SE_adlasso2_beta2 <- c(data2[["k5"]][["SE_adlasso2_beta2"]], data2[["k20"]][["SE_adlasso2_beta2"]], data2[["k100"]][["SE_adlasso2_beta2"]])
)
colnames(SE_beta2) <- c("Ridge", "LASSO", "AdLASSO1", "AdLASSO2")
rownames(SE_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(SE_beta2)

power_beta1 <- data.frame(
  power_lasso_beta1 <- c(data2[["k5"]][["power_lasso_beta1"]], data2[["k20"]][["power_lasso_beta1"]], data2[["k100"]][["power_lasso_beta1"]]),
  power_adlasso1_beta1 <- c(data2[["k5"]][["power_adlasso1_beta1"]], data2[["k20"]][["power_adlasso1_beta1"]], data2[["k100"]][["power_adlasso1_beta1"]]),
  power_adlasso2_beta1 <- c(data2[["k5"]][["power_adlasso2_beta1"]], data2[["k20"]][["power_adlasso2_beta1"]], data2[["k100"]][["power_adlasso2_beta1"]]),
  power_ridgekn_beta1 <- c(data2[["k5"]][["power_ridgekn_beta1"]], data2[["k20"]][["power_ridgekn_beta1"]], data2[["k100"]][["power_ridgekn_beta1"]]),
  power_lassokn_beta1 <- c(data2[["k5"]][["power_lassokn_beta1"]], data2[["k20"]][["power_lassokn_beta1"]], data2[["k100"]][["power_lassokn_beta1"]])
) 

colnames(power_beta1) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(power_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(power_beta1)

power_beta2 <- data.frame(
  power_lasso_beta2 <- c(data2[["k5"]][["power_lasso_beta2"]], data2[["k20"]][["power_lasso_beta2"]], data2[["k100"]][["power_lasso_beta2"]]),
  power_adlasso1_beta2 <- c(data2[["k5"]][["power_adlasso1_beta2"]], data2[["k20"]][["power_adlasso1_beta2"]], data2[["k100"]][["power_adlasso1_beta2"]]),
  power_adlasso2_beta2 <- c(data2[["k5"]][["power_adlasso2_beta2"]], data2[["k20"]][["power_adlasso2_beta2"]], data2[["k100"]][["power_adlasso2_beta2"]]),
  power_ridgekn_beta2 <- c(data2[["k5"]][["power_ridgekn_beta2"]], data2[["k20"]][["power_ridgekn_beta2"]], data2[["k100"]][["power_ridgekn_beta2"]]),
  power_lassokn_beta2 <- c(data2[["k5"]][["power_lassokn_beta2"]], data2[["k20"]][["power_lassokn_beta2"]], data2[["k100"]][["power_lassokn_beta2"]])
) 

colnames(power_beta2) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(power_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(power_beta2)

FDR_beta1 <- data.frame(
  FDR_lasso_beta1 <- c(data2[["k5"]][["FDR_lasso_beta1"]], data2[["k20"]][["FDR_lasso_beta1"]], data2[["k100"]][["FDR_lasso_beta1"]]),
  FDR_adlasso1_beta1 <- c(data2[["k5"]][["FDR_adlasso1_beta1"]], data2[["k20"]][["FDR_adlasso1_beta1"]], data2[["k100"]][["FDR_adlasso1_beta1"]]),
  FDR_adlasso2_beta1 <- c(data2[["k5"]][["FDR_adlasso2_beta1"]], data2[["k20"]][["FDR_adlasso2_beta1"]], data2[["k100"]][["FDR_adlasso2_beta1"]]),
  FDR_ridgekn_beta1 <- c(data2[["k5"]][["FDR_ridgekn_beta1"]], data2[["k20"]][["FDR_ridgekn_beta1"]], data2[["k100"]][["FDR_ridgekn_beta1"]]),
  FDR_lassokn_beta1 <- c(data2[["k5"]][["FDR_lassokn_beta1"]], data2[["k20"]][["FDR_lassokn_beta1"]], data2[["k100"]][["FDR_lassokn_beta1"]])
)

colnames(FDR_beta1) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(FDR_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(FDR_beta1)

FDR_beta2 <- data.frame(
  FDR_lasso_beta2 <- c(data2[["k5"]][["FDR_lasso_beta2"]], data2[["k20"]][["FDR_lasso_beta2"]], data2[["k100"]][["FDR_lasso_beta2"]]),
  FDR_adlasso1_beta2 <- c(data2[["k5"]][["FDR_adlasso1_beta2"]], data2[["k20"]][["FDR_adlasso1_beta2"]], data2[["k100"]][["FDR_adlasso1_beta2"]]),
  FDR_adlasso2_beta2 <- c(data2[["k5"]][["FDR_adlasso2_beta2"]], data2[["k20"]][["FDR_adlasso2_beta2"]], data2[["k100"]][["FDR_adlasso2_beta2"]]),
  FDR_ridgekn_beta2 <- c(data2[["k5"]][["FDR_ridgekn_beta2"]], data2[["k20"]][["FDR_ridgekn_beta2"]], data2[["k100"]][["FDR_ridgekn_beta2"]]),
  FDR_lassokn_beta2 <- c(data2[["k5"]][["FDR_lassokn_beta2"]], data2[["k20"]][["FDR_lassokn_beta2"]], data2[["k100"]][["FDR_lassokn_beta2"]])
)

colnames(FDR_beta2) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(FDR_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(FDR_beta2)

rep <- 100

results <- replicate(rep, {
  Y <- lapply(beta, function(b) {
    Y_beta1 <- X%*%b[["beta1"]] + rnorm(500, 0, 1)
    Y_beta2 <- X%*%b[["beta2"]] + rnorm(500, 0, 1)
    return(list(Y_beta1 = Y_beta1, Y_beta2 = Y_beta2))
  })
  
  data2 <- setNames(lapply(names(Y), function(name_k) {
    y <- Y[[name_k]]
    
    ridge_fit_2_beta1 <- cv.glmnet(X, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 0)
    hbeta_ridge_2_beta1 <- coef(ridge_fit_2_beta1, s = 'lambda.min')[2:501]
    ridge_fit_2_beta2 <- cv.glmnet(X, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 0)
    hbeta_ridge_2_beta2 <- coef(ridge_fit_2_beta2, s = 'lambda.min')[2:501]
    
    lasso_fit_2_beta1 <- cv.glmnet(X, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbeta_lasso_2_beta1 <- coef(lasso_fit_2_beta1, s = 'lambda.min')[2:501]
    lasso_fit_2_beta2 <- cv.glmnet(X, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbeta_lasso_2_beta2 <- coef(lasso_fit_2_beta2, s = 'lambda.min')[2:501]
    
    W1_beta1 <- 1/abs(hbeta_ridge_2_beta1)
    X_temp_beta1 <- X %*% diag(1/W1_beta1)
    adlasso1_fit_beta1 <- cv.glmnet(X_temp_beta1, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbetatemp_adlasso1_beta1 <- coef(adlasso1_fit_beta1, s = 'lambda.min')[2:501]
    hbeta_adlasso1_beta1 <- hbetatemp_adlasso1_beta1/W1_beta1
    W1_beta2 <- 1/(abs(hbeta_ridge_2_beta2))
    X_temp_beta2 <- X %*% diag(1/W1_beta2)
    adlasso1_fit_beta2 <- cv.glmnet(X_temp_beta2, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbetatemp_adlasso1_beta2 <- coef(adlasso1_fit_beta2, s = 'lambda.min')[2:501]
    hbeta_adlasso1_beta2 <- hbetatemp_adlasso1_beta2/W1_beta2
    
    RSS_beta1 <- sum((y[["Y_beta1"]] - X %*% hbeta_lasso_2_beta1)^2)
    ind_lasso_beta1 <- which(hbeta_lasso_2_beta1 > 0)
    if (length(ind_lasso_beta1) == 0) {
      hbeta_adlasso2_beta1 = rep(0, 500)
    } else {
      hsigma_beta1 <- sqrt(RSS_beta1/(500 - length(ind_lasso_beta1)))
      W2_beta1 <- hsigma_beta1/abs(hbeta_lasso_2_beta1[ind_lasso_beta1])
      lambda1_beta1 <- hsigma_beta1 * qnorm(1 - 0.2/(2 * 500))
      X_small_beta1 <- X[, ind_lasso_beta1]
      if (length(ind_lasso_beta1) == 1) {
        X_temp_beta1 <- matrix(X_small_beta1 / W2_beta1, ncol = 1)
        adlasso2_fit_beta1 <- lm(y[["Y_beta1"]] ~ X_temp_beta1 + 0)
        hbetatemp_adlasso2_beta1 <- coef(adlasso2_fit_beta1)
      } else {
        X_temp_beta1 <- X_small_beta1 %*% diag(1 / W2_beta1)
        adlasso2_fit_beta1 <- glmnet(X_temp_beta1, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1, lambda = lambda1_beta1/500)
        hbetatemp_adlasso2_beta1 <- coef(adlasso2_fit_beta1)[2:(length(ind_lasso_beta1) + 1)]
      }
      hbeta_adlasso2_beta1 <- rep(0, 500)
      hbeta_adlasso2_beta1[ind_lasso_beta1] <- hbetatemp_adlasso2_beta1/W2_beta1
    }
    RSS_beta2 <- sum((y[["Y_beta2"]] - X %*% hbeta_lasso_2_beta2)^2)
    ind_lasso_beta2 <- which(hbeta_lasso_2_beta2 > 0)
    if (length(ind_lasso_beta2) == 0) {
      hbeta_adlasso2_beta2 = rep(0, 500)
    } else {
      hsigma_beta2 <- sqrt(RSS_beta2/(500 - length(ind_lasso_beta2)))
      W2_beta2 <- hsigma_beta2/abs(hbeta_lasso_2_beta2[ind_lasso_beta2])
      lambda1_beta2 <- hsigma_beta2 * qnorm(1 - 0.2/(2 * 500))
      X_small_beta2 <- X[, ind_lasso_beta2]
      if (length(ind_lasso_beta2) == 1) {
        X_temp_beta2 <- matrix(X_small_beta2 / W2_beta2, ncol = 1)
        adlasso2_fit_beta2 <- lm(y[["Y_beta2"]] ~ X_temp_beta2 + 0)
        hbetatemp_adlasso2_beta2 <- coef(adlasso2_fit_beta2)
      } else {
        X_temp_beta2 <- X_small_beta2 %*% diag(1 / W2_beta2)
        adlasso2_fit_beta2 <- glmnet(X_temp_beta2, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1, lambda = lambda1_beta2/500)
        hbetatemp_adlasso2_beta2 <- coef(adlasso2_fit_beta2)[2:(length(ind_lasso_beta2) +1)]
      }
      
      hbeta_adlasso2_beta2 <- rep(0, 500)
      hbeta_adlasso2_beta2[ind_lasso_beta2] <- hbetatemp_adlasso2_beta2/W2_beta2
    }
    
    X_copy <- matrix(rnorm(500 * 500, 0, 1/sqrt(500)), nrow = 500)
    X_large <- cbind(X, X_copy)
    large_fit_ridge_beta1 <- cv.glmnet(X_large, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 0)
    hbeta_large_ridge_beta1 <- coef(large_fit_ridge_beta1, s='lambda.min')[2:(2*500 + 1)]
    w <- abs(hbeta_large_ridge_beta1[1:500]) - abs(hbeta_large_ridge_beta1[501:1000])
    
    kn <- function(t) {
      ((1 + sum(w <= -t))/max(sum(w >= t), 1))
    }
    
    w_kn <- sapply(c(0, abs(w)), kn)
    thres <- min(c(0, abs(w))[w_kn <= 0.2])
    hbeta_knridge_beta1 <- hbeta_ridge_2_beta1 * (w >= thres)
    
    large_fit_ridge_beta2 <- cv.glmnet(X_large, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 0)
    hbeta_large_ridge_beta2 <- coef(large_fit_ridge_beta2, s='lambda.min')[2:(2*500 + 1)]
    w <- abs(hbeta_large_ridge_beta2[1:500]) - abs(hbeta_large_ridge_beta2[501:1000])
    w_kn <- sapply(c(0, abs(w)), kn)
    thres <- min(c(0, abs(w))[w_kn <= 0.2])
    hbeta_knridge_beta2 <- hbeta_ridge_2_beta2 * (w >= thres)
    
    large_fit_lasso_beta1 <- cv.glmnet(X_large, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbeta_large_lasso_beta1 <- coef(large_fit_lasso_beta1, s = 'lambda.min')[2:1001]
    w <- abs(hbeta_large_lasso_beta1[1:500]) - abs(hbeta_large_lasso_beta1[501:1000])
    w_kn_lasso <- sapply(c(0, abs(w)), kn)
    thres_lasso <- min(c(0, abs(w))[w_kn <= 0.2])
    hbeta_knlasso_beta1 <- hbeta_lasso_2_beta1 * (w >= thres)
    large_fit_lasso_beta2 <- cv.glmnet(X_large, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbeta_large_lasso_beta2 <- coef(large_fit_lasso_beta2, s = 'lambda.min')[2:1001]
    w <- abs(hbeta_large_lasso_beta2[1:500]) - abs(hbeta_large_lasso_beta2[501:1000])
    w_kn_lasso <- sapply(c(0, abs(w)), kn)
    thres_lasso <- min(c(0, abs(w))[w_kn <= 0.2])
    hbeta_knlasso_beta2 <- hbeta_lasso_2_beta2 * (w >= thres)
    
    SE_ridge_beta1 <- sum((hbeta_ridge_2_beta1 - beta[[name_k]][["beta1"]])^2)
    SE_ridge_beta2 <- sum((hbeta_ridge_2_beta2 - beta[[name_k]][["beta2"]])^2)
    SE_lasso_beta1 <- sum((hbeta_lasso_2_beta1 - beta[[name_k]][["beta1"]])^2)
    SE_lasso_beta2 <- sum((hbeta_lasso_2_beta2 - beta[[name_k]][["beta2"]])^2)
    SE_adlasso1_beta1 <- sum((hbeta_adlasso1_beta1 - beta[[name_k]][["beta1"]])^2)
    SE_adlasso1_beta2 <- sum((hbeta_adlasso1_beta2 - beta[[name_k]][["beta2"]])^2)
    SE_adlasso2_beta1 <- sum((hbeta_adlasso2_beta1 - beta[[name_k]][["beta1"]])^2)
    SE_adlasso2_beta2 <- sum((hbeta_adlasso2_beta2 - beta[[name_k]][["beta2"]])^2)
    
    true_beta <- beta[[name_k]][["beta1"]] != 0
    k <- sum(true_beta)
    
    lasso_selected_beta1 <- hbeta_lasso_2_beta1 != 0
    TP_lasso_beta1 <- sum(true_beta & lasso_selected_beta1)
    FD_lasso_beta1 <- sum(!true_beta & lasso_selected_beta1)
    R_lasso_beta1 <- sum(lasso_selected_beta1)
    power_lasso_beta1 <- TP_lasso_beta1 / k
    FDR_lasso_beta1 <- if (R_lasso_beta1 == 0) 0 else FD_lasso_beta1 / R_lasso_beta1
    lasso_selected_beta2 <- hbeta_lasso_2_beta2 != 0
    TP_lasso_beta2 <- sum(true_beta & lasso_selected_beta2)
    FD_lasso_beta2 <- sum(!true_beta & lasso_selected_beta2)
    R_lasso_beta2 <- sum(lasso_selected_beta2)
    power_lasso_beta2 <- TP_lasso_beta2 / k
    FDR_lasso_beta2 <- if (R_lasso_beta2 == 0) 0 else FD_lasso_beta2 / R_lasso_beta2
    
    adlasso1_selected_beta1 <- hbeta_adlasso1_beta1 != 0
    TP_adlasso1_beta1 <- sum(true_beta & adlasso1_selected_beta1)
    FD_adlasso1_beta1 <- sum(!true_beta & adlasso1_selected_beta1)
    R_adlasso1_beta1 <- sum(adlasso1_selected_beta1)
    power_adlasso1_beta1 <- TP_adlasso1_beta1 / k
    FDR_adlasso1_beta1 <- if (R_adlasso1_beta1 == 0) 0 else FD_adlasso1_beta1 / R_adlasso1_beta1
    adlasso1_selected_beta2 <- hbeta_adlasso1_beta2 != 0
    TP_adlasso1_beta2 <- sum(true_beta & adlasso1_selected_beta2)
    FD_adlasso1_beta2 <- sum(!true_beta & adlasso1_selected_beta2)
    R_adlasso1_beta2 <- sum(adlasso1_selected_beta2)
    power_adlasso1_beta2 <- TP_adlasso1_beta2 / k
    FDR_adlasso1_beta2 <- if (R_adlasso1_beta2 == 0) 0 else FD_adlasso1_beta2 / R_adlasso1_beta2
    
    adlasso2_selected_beta1 <- hbeta_adlasso2_beta1 != 0
    TP_adlasso2_beta1 <- sum(true_beta & adlasso2_selected_beta1)
    FD_adlasso2_beta1 <- sum(!true_beta & adlasso2_selected_beta1)
    R_adlasso2_beta1 <- sum(adlasso2_selected_beta1)
    power_adlasso2_beta1 <- TP_adlasso2_beta1 / k
    FDR_adlasso2_beta1 <- if (R_adlasso2_beta1 == 0) 0 else FD_adlasso2_beta1 / R_adlasso2_beta1
    adlasso2_selected_beta2 <- hbeta_adlasso2_beta2 != 0
    TP_adlasso2_beta2 <- sum(true_beta & adlasso2_selected_beta2)
    FD_adlasso2_beta2 <- sum(!true_beta & adlasso2_selected_beta2)
    R_adlasso2_beta2 <- sum(adlasso2_selected_beta2)
    power_adlasso2_beta2 <- TP_adlasso2_beta2 / k
    FDR_adlasso2_beta2 <- if (R_adlasso2_beta2 == 0) 0 else FD_adlasso2_beta2 / R_adlasso2_beta2
    
    ridgekn_selected_beta1 <- hbeta_knridge_beta1 != 0
    TP_ridgekn_beta1 <- sum(true_beta & ridgekn_selected_beta1)
    FD_ridgekn_beta1 <- sum(!true_beta & ridgekn_selected_beta1)
    R_ridgekn_beta1 <- sum(ridgekn_selected_beta1)
    power_ridgekn_beta1 <- TP_ridgekn_beta1 / k
    FDR_ridgekn_beta1 <- if (R_ridgekn_beta1 == 0) 0 else FD_ridgekn_beta1 / R_ridgekn_beta1
    ridgekn_selected_beta2 <- hbeta_knridge_beta2 != 0
    TP_ridgekn_beta2 <- sum(true_beta & ridgekn_selected_beta2)
    FD_ridgekn_beta2 <- sum(!true_beta & ridgekn_selected_beta2)
    R_ridgekn_beta2 <- sum(ridgekn_selected_beta2)
    power_ridgekn_beta2 <- TP_ridgekn_beta2 / k
    FDR_ridgekn_beta2 <- if (R_ridgekn_beta2 == 0) 0 else FD_ridgekn_beta2 / R_ridgekn_beta2
    
    lassokn_selected_beta1 <- hbeta_knlasso_beta1 != 0
    TP_lassokn_beta1 <- sum(true_beta & lassokn_selected_beta1)
    FD_lassokn_beta1 <- sum(!true_beta & lassokn_selected_beta1)
    R_lassokn_beta1 <- sum(lassokn_selected_beta1)
    power_lassokn_beta1 <- TP_lassokn_beta1 / k
    FDR_lassokn_beta1 <- if (R_lassokn_beta1 == 0) 0 else FD_lassokn_beta1 / R_lassokn_beta1
    lassokn_selected_beta2 <- hbeta_knlasso_beta2 != 0
    TP_lassokn_beta2 <- sum(true_beta & lassokn_selected_beta2)
    FD_lassokn_beta2 <- sum(!true_beta & lassokn_selected_beta2)
    R_lassokn_beta2 <- sum(lassokn_selected_beta2)
    power_lassokn_beta2 <- TP_lassokn_beta2 / k
    FDR_lassokn_beta2 <- if (R_lassokn_beta2 == 0) 0 else FD_lassokn_beta2 / R_lassokn_beta2
    
    return(list(
      hbeta_ridge_beta1 = hbeta_ridge_2_beta1,
      SE_ridge_beta1 = SE_ridge_beta1,
      hbeta_ridge_beta2 = hbeta_ridge_2_beta2,
      SE_ridge_beta2 = SE_ridge_beta2,
      hbeta_lasso_beta1 = hbeta_lasso_2_beta1,
      SE_lasso_beta1 = SE_lasso_beta1,
      power_lasso_beta1 = power_lasso_beta1,
      FDR_lasso_beta1 = FDR_lasso_beta1,
      hbeta_lasso_beta2 = hbeta_lasso_2_beta2,
      SE_lasso_beta2 = SE_lasso_beta2,
      power_lasso_beta2 = power_lasso_beta2,
      FDR_lasso_beta2 = FDR_lasso_beta2,
      hbeta_adlasso1_beta1 = hbeta_adlasso1_beta1,
      SE_adlasso1_beta1 = SE_adlasso1_beta1,
      power_adlasso1_beta1 = power_adlasso1_beta1,
      FDR_adlasso1_beta1 = FDR_adlasso1_beta1,
      hbeta_adlasso1_beta2 = hbeta_adlasso1_beta2,
      SE_adlasso1_beta2 = SE_adlasso1_beta2,
      power_adlasso1_beta2 = power_adlasso1_beta2,
      FDR_adlasso1_beta2 = FDR_adlasso1_beta2,
      hbeta_adlasso2_beta1 = hbeta_adlasso2_beta1,
      SE_adlasso2_beta1 = SE_adlasso2_beta1,
      power_adlasso2_beta1 = power_adlasso2_beta1,
      FDR_adlasso2_beta1 = FDR_adlasso2_beta1,
      hbeta_adlasso2_beta2 = hbeta_adlasso2_beta2,
      SE_adlasso2_beta2 = SE_adlasso2_beta2,
      power_adlasso2_beta2 = power_adlasso2_beta2,
      FDR_adlasso2_beta2 = FDR_adlasso2_beta2,
      hbeta_knridge_beta1 = hbeta_knridge_beta1,
      power_ridgekn_beta1 = power_ridgekn_beta1,
      FDR_ridgekn_beta1 = FDR_ridgekn_beta1,
      hbeta_knridge_beta2 = hbeta_knridge_beta2,
      power_ridgekn_beta2 = power_ridgekn_beta2,
      FDR_ridgekn_beta2 = FDR_ridgekn_beta2,
      hbeta_knlasso_beta1 = hbeta_knlasso_beta1,
      power_lassokn_beta1 = power_lassokn_beta1,
      FDR_lassokn_beta1 = FDR_lassokn_beta1,
      hbeta_knlasso_beta2 = hbeta_knlasso_beta2,
      power_lassokn_beta2 = power_lassokn_beta2,
      FDR_lassokn_beta2 = FDR_lassokn_beta2)
    )
  }), names(beta))
  
  return(data2)
}, simplify = FALSE)

saveRDS(results, "results.rds")

library(dplyr)

extract_metrics <- function(results) {
  k_values <- c("k5", "k20", "k100")
  all_k_metrics <- list()
  
  for (k in k_values) {
    metrics_list <- lapply(results, function(res) {
      res_k <- res[[k]]
      data.frame(
        SE_ridge_beta1 = res_k$SE_ridge_beta1,
        SE_ridge_beta2 = res_k$SE_ridge_beta2,
        SE_lasso_beta1 = res_k$SE_lasso_beta1,
        SE_lasso_beta2 = res_k$SE_lasso_beta2,
        SE_adlasso1_beta1 = res_k$SE_adlasso1_beta1,
        SE_adlasso1_beta2 = res_k$SE_adlasso1_beta2,
        SE_adlasso2_beta1 = res_k$SE_adlasso2_beta1,
        SE_adlasso2_beta2 = res_k$SE_adlasso2_beta2,
        FDR_lasso_beta1 = res_k$FDR_lasso_beta1,
        FDR_lasso_beta2 = res_k$FDR_lasso_beta2,
        power_lasso_beta1 = res_k$power_lasso_beta1,
        power_lasso_beta2 = res_k$power_lasso_beta2,
        FDR_adlasso1_beta1 = res_k$FDR_adlasso1_beta1,
        FDR_adlasso1_beta2 = res_k$FDR_adlasso1_beta2,
        power_adlasso1_beta1 = res_k$power_adlasso1_beta1,
        power_adlasso1_beta2 = res_k$power_adlasso1_beta2,
        FDR_adlasso2_beta1 = res_k$FDR_adlasso2_beta1,
        FDR_adlasso2_beta2 = res_k$FDR_adlasso2_beta2,
        power_adlasso2_beta1 = res_k$power_adlasso2_beta1,
        power_adlasso2_beta2 = res_k$power_adlasso2_beta2,
        FDR_ridgekn_beta1 = res_k$FDR_ridgekn_beta1,
        FDR_ridgekn_beta2 = res_k$FDR_ridgekn_beta2,
        power_ridgekn_beta1 = res_k$power_ridgekn_beta1,
        power_ridgekn_beta2 = res_k$power_ridgekn_beta2,
        FDR_lassokn_beta1 = res_k$FDR_lassokn_beta1,
        FDR_lassokn_beta2 = res_k$FDR_lassokn_beta2,
        power_lassokn_beta1 = res_k$power_lassokn_beta1,
        power_lassokn_beta2 = res_k$power_lassokn_beta2
      )
    })
    
    k_df <- bind_rows(metrics_list)
    k_avg <- colMeans(k_df, na.rm = TRUE)
    all_k_metrics[[k]] <- k_avg
  }
  
  return(all_k_metrics)
}

averaged_metrics_by_k <- extract_metrics(results)

MSE_beta1 <- data.frame(
  SE_ridge_beta1 <- c(averaged_metrics_by_k[["k5"]][["SE_ridge_beta1"]], averaged_metrics_by_k[["k20"]][["SE_ridge_beta1"]], averaged_metrics_by_k[["k100"]][["SE_ridge_beta1"]]),
  SE_lasso_beta1 <- c(averaged_metrics_by_k[["k5"]][["SE_lasso_beta1"]], averaged_metrics_by_k[["k20"]][["SE_lasso_beta1"]], averaged_metrics_by_k[["k100"]][["SE_lasso_beta1"]]),
  SE_adlasso1_beta1 <- c(averaged_metrics_by_k[["k5"]][["SE_adlasso1_beta1"]], averaged_metrics_by_k[["k20"]][["SE_adlasso1_beta1"]], averaged_metrics_by_k[["k100"]][["SE_adlasso1_beta1"]]),
  SE_adlasso2_beta1 <- c(averaged_metrics_by_k[["k5"]][["SE_adlasso2_beta1"]], averaged_metrics_by_k[["k20"]][["SE_adlasso2_beta1"]], averaged_metrics_by_k[["k100"]][["SE_adlasso2_beta1"]])
)
colnames(MSE_beta1) <- c("Ridge", "LASSO", "AdLASSO1", "AdLASSO2")
rownames(MSE_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(MSE_beta1)

MSE_beta2<- data.frame(
  SE_ridge_beta2 <- c(averaged_metrics_by_k[["k5"]][["SE_ridge_beta2"]], averaged_metrics_by_k[["k20"]][["SE_ridge_beta2"]], averaged_metrics_by_k[["k100"]][["SE_ridge_beta2"]]),
  SE_lasso_beta2 <- c(averaged_metrics_by_k[["k5"]][["SE_lasso_beta2"]], averaged_metrics_by_k[["k20"]][["SE_lasso_beta2"]], averaged_metrics_by_k[["k100"]][["SE_lasso_beta2"]]),
  SE_adlasso1_beta2 <- c(averaged_metrics_by_k[["k5"]][["SE_adlasso1_beta2"]], averaged_metrics_by_k[["k20"]][["SE_adlasso1_beta2"]], averaged_metrics_by_k[["k100"]][["SE_adlasso1_beta2"]]),
  SE_adlasso2_beta2 <- c(averaged_metrics_by_k[["k5"]][["SE_adlasso2_beta2"]], averaged_metrics_by_k[["k20"]][["SE_adlasso2_beta2"]], averaged_metrics_by_k[["k100"]][["SE_adlasso2_beta2"]])
)
colnames(MSE_beta2) <- c("Ridge", "LASSO", "AdLASSO1", "AdLASSO2")
rownames(MSE_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(MSE_beta2)

power_beta1 <- data.frame(
  power_lasso_beta1 <- c(averaged_metrics_by_k[["k5"]][["power_lasso_beta1"]], averaged_metrics_by_k[["k20"]][["power_lasso_beta1"]], averaged_metrics_by_k[["k100"]][["power_lasso_beta1"]]),
  power_adlasso1_beta1 <- c(averaged_metrics_by_k[["k5"]][["power_adlasso1_beta1"]], averaged_metrics_by_k[["k20"]][["power_adlasso1_beta1"]], averaged_metrics_by_k[["k100"]][["power_adlasso1_beta1"]]),
  power_adlasso2_beta1 <- c(averaged_metrics_by_k[["k5"]][["power_adlasso2_beta1"]], averaged_metrics_by_k[["k20"]][["power_adlasso2_beta1"]], averaged_metrics_by_k[["k100"]][["power_adlasso2_beta1"]]),
  power_ridgekn_beta1 <- c(averaged_metrics_by_k[["k5"]][["power_ridgekn_beta1"]], averaged_metrics_by_k[["k20"]][["power_ridgekn_beta1"]], averaged_metrics_by_k[["k100"]][["power_ridgekn_beta1"]]),
  power_lassokn_beta1 <- c(averaged_metrics_by_k[["k5"]][["power_lassokn_beta1"]], averaged_metrics_by_k[["k20"]][["power_lassokn_beta1"]], averaged_metrics_by_k[["k100"]][["power_lassokn_beta1"]])
) 

colnames(power_beta1) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(power_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(power_beta1)

power_beta2 <- data.frame(
  power_lasso_beta2 <- c(averaged_metrics_by_k[["k5"]][["power_lasso_beta2"]], averaged_metrics_by_k[["k20"]][["power_lasso_beta2"]], averaged_metrics_by_k[["k100"]][["power_lasso_beta2"]]),
  power_adlasso1_beta2 <- c(averaged_metrics_by_k[["k5"]][["power_adlasso1_beta2"]], averaged_metrics_by_k[["k20"]][["power_adlasso1_beta2"]], averaged_metrics_by_k[["k100"]][["power_adlasso1_beta2"]]),
  power_adlasso2_beta2 <- c(averaged_metrics_by_k[["k5"]][["power_adlasso2_beta2"]], averaged_metrics_by_k[["k20"]][["power_adlasso2_beta2"]], averaged_metrics_by_k[["k100"]][["power_adlasso2_beta2"]]),
  power_ridgekn_beta2 <- c(averaged_metrics_by_k[["k5"]][["power_ridgekn_beta2"]], averaged_metrics_by_k[["k20"]][["power_ridgekn_beta2"]], averaged_metrics_by_k[["k100"]][["power_ridgekn_beta2"]]),
  power_lassokn_beta2 <- c(averaged_metrics_by_k[["k5"]][["power_lassokn_beta2"]], averaged_metrics_by_k[["k20"]][["power_lassokn_beta2"]], averaged_metrics_by_k[["k100"]][["power_lassokn_beta2"]])
) 

colnames(power_beta2) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(power_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(power_beta2)

FDR_beta1 <- data.frame(
  FDR_lasso_beta1 <- c(averaged_metrics_by_k[["k5"]][["FDR_lasso_beta1"]], averaged_metrics_by_k[["k20"]][["FDR_lasso_beta1"]], averaged_metrics_by_k[["k100"]][["FDR_lasso_beta1"]]),
  FDR_adlasso1_beta1 <- c(averaged_metrics_by_k[["k5"]][["FDR_adlasso1_beta1"]], averaged_metrics_by_k[["k20"]][["FDR_adlasso1_beta1"]], averaged_metrics_by_k[["k100"]][["FDR_adlasso1_beta1"]]),
  FDR_adlasso2_beta1 <- c(averaged_metrics_by_k[["k5"]][["FDR_adlasso2_beta1"]], averaged_metrics_by_k[["k20"]][["FDR_adlasso2_beta1"]], averaged_metrics_by_k[["k100"]][["FDR_adlasso2_beta1"]]),
  FDR_ridgekn_beta1 <- c(averaged_metrics_by_k[["k5"]][["FDR_ridgekn_beta1"]], averaged_metrics_by_k[["k20"]][["FDR_ridgekn_beta1"]], averaged_metrics_by_k[["k100"]][["FDR_ridgekn_beta1"]]),
  FDR_lassokn_beta1 <- c(averaged_metrics_by_k[["k5"]][["FDR_lassokn_beta1"]], averaged_metrics_by_k[["k20"]][["FDR_lassokn_beta1"]], averaged_metrics_by_k[["k100"]][["FDR_lassokn_beta1"]])
)

colnames(FDR_beta1) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(FDR_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(FDR_beta1)

FDR_beta2 <- data.frame(
  FDR_lasso_beta2 <- c(averaged_metrics_by_k[["k5"]][["FDR_lasso_beta2"]], averaged_metrics_by_k[["k20"]][["FDR_lasso_beta2"]], averaged_metrics_by_k[["k100"]][["FDR_lasso_beta2"]]),
  FDR_adlasso1_beta2 <- c(averaged_metrics_by_k[["k5"]][["FDR_adlasso1_beta2"]], averaged_metrics_by_k[["k20"]][["FDR_adlasso1_beta2"]], averaged_metrics_by_k[["k100"]][["FDR_adlasso1_beta2"]]),
  FDR_adlasso2_beta2 <- c(averaged_metrics_by_k[["k5"]][["FDR_adlasso2_beta2"]], averaged_metrics_by_k[["k20"]][["FDR_adlasso2_beta2"]], averaged_metrics_by_k[["k100"]][["FDR_adlasso2_beta2"]]),
  FDR_ridgekn_beta2 <- c(averaged_metrics_by_k[["k5"]][["FDR_ridgekn_beta2"]], averaged_metrics_by_k[["k20"]][["FDR_ridgekn_beta2"]], averaged_metrics_by_k[["k100"]][["FDR_ridgekn_beta2"]]),
  FDR_lassokn_beta2 <- c(averaged_metrics_by_k[["k5"]][["FDR_lassokn_beta2"]], averaged_metrics_by_k[["k20"]][["FDR_lassokn_beta2"]], averaged_metrics_by_k[["k100"]][["FDR_lassokn_beta2"]])
)

colnames(FDR_beta2) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(FDR_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(FDR_beta2)

#zadanie 3.
library(mvtnorm)
library(MASS)
n <- 500
p <- 500
Sigma <- matrix(0.5, nrow = n, ncol = p)
diag(Sigma) <- 1
X <- rmvnorm(n, mean = rep(0, 500), sigma = Sigma/n)

Y <- lapply(beta, function(b) {
  Y_beta1 <- X%*%b[["beta1"]] + rnorm(500, 0, 1)
  Y_beta2 <- X%*%b[["beta2"]] + rnorm(500, 0, 1)
  return(list(Y_beta1 = Y_beta1, Y_beta2 = Y_beta2))
})

data3 <- setNames(lapply(names(Y), function(name_k) {
  y <- Y[[name_k]]
  
  ridge_fit_2_beta1 <- cv.glmnet(X, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 0)
  hbeta_ridge_2_beta1 <- coef(ridge_fit_2_beta1, s = 'lambda.min')[2:501]
  ridge_fit_2_beta2 <- cv.glmnet(X, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 0)
  hbeta_ridge_2_beta2 <- coef(ridge_fit_2_beta2, s = 'lambda.min')[2:501]
  
  lasso_fit_2_beta1 <- cv.glmnet(X, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbeta_lasso_2_beta1 <- coef(lasso_fit_2_beta1, s = 'lambda.min')[2:501]
  lasso_fit_2_beta2 <- cv.glmnet(X, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbeta_lasso_2_beta2 <- coef(lasso_fit_2_beta2, s = 'lambda.min')[2:501]
  
  W1_beta1 <- 1/abs(hbeta_ridge_2_beta1)
  X_temp_beta1 <- X %*% diag(1/W1_beta1)
  adlasso1_fit_beta1 <- cv.glmnet(X_temp_beta1, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbetatemp_adlasso1_beta1 <- coef(adlasso1_fit_beta1, s = 'lambda.min')[2:501]
  hbeta_adlasso1_beta1 <- hbetatemp_adlasso1_beta1/W1_beta1
  W1_beta2 <- 1/(abs(hbeta_ridge_2_beta2))
  X_temp_beta2 <- X %*% diag(1/W1_beta2)
  adlasso1_fit_beta2 <- cv.glmnet(X_temp_beta2, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbetatemp_adlasso1_beta2 <- coef(adlasso1_fit_beta2, s = 'lambda.min')[2:501]
  hbeta_adlasso1_beta2 <- hbetatemp_adlasso1_beta2/W1_beta2
  
  RSS_beta1 <- sum((y[["Y_beta1"]] - X %*% hbeta_lasso_2_beta1)^2)
  ind_lasso_beta1 <- which(hbeta_lasso_2_beta1 > 0)
  if (length(ind_lasso_beta1) == 0) {
    hbeta_adlasso2_beta1 = rep(0, 500)
  } else {
    hsigma_beta1 <- sqrt(RSS_beta1/(500 - length(ind_lasso_beta1)))
    W2_beta1 <- hsigma_beta1/abs(hbeta_lasso_2_beta1[ind_lasso_beta1])
    lambda1_beta1 <- hsigma_beta1 * qnorm(1 - 0.2/(2 * 500))
    X_small_beta1 <- X[, ind_lasso_beta1]
    X_temp_beta1 <- X_small_beta1 %*% diag(1/W2_beta1)
    
    adlasso2_fit_beta1 <- glmnet(X_temp_beta1, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1, lambda = lambda1_beta1/500)
    hbetatemp_adlasso2_beta1 <- coef(adlasso2_fit_beta1)[2:(length(ind_lasso_beta1) +1)]
    hbeta_adlasso2_beta1 <- rep(0, 500)
    hbeta_adlasso2_beta1[ind_lasso_beta1] <- hbetatemp_adlasso2_beta1/W2_beta1
  }
  RSS_beta2 <- sum((y[["Y_beta2"]] - X %*% hbeta_lasso_2_beta2)^2)
  ind_lasso_beta2 <- which(hbeta_lasso_2_beta2 > 0)
  if (length(ind_lasso_beta2) == 0) {
    hbeta_adlasso2_beta2 = rep(0, 500)
  } else {
    hsigma_beta2 <- sqrt(RSS_beta2/(500 - length(ind_lasso_beta2)))
    W2_beta2 <- hsigma_beta2/abs(hbeta_lasso_2_beta2[ind_lasso_beta2])
    lambda1_beta2 <- hsigma_beta2 * qnorm(1 - 0.2/(2 * 500))
    X_small_beta2 <- X[, ind_lasso_beta2]
    X_temp_beta2 <- X_small_beta2 %*% diag(1/W2_beta2)
    
    adlasso2_fit_beta2 <- glmnet(X_temp_beta2, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1, lambda = lambda1_beta2/500)
    hbetatemp_adlasso2_beta2 <- coef(adlasso2_fit_beta2)[2:(length(ind_lasso_beta2) +1)]
    hbeta_adlasso2_beta2 <- rep(0, 500)
    hbeta_adlasso2_beta2[ind_lasso_beta2] <- hbetatemp_adlasso2_beta2/W2_beta2
  }
  
  lambda_min <- min(eigen(Sigma)$values)
  s <- rep(2 * lambda_min, p)
  s <- pmin(s, 1)  
  diag_s <- diag(s)
  
  G_11 <- Sigma
  G_12 <- Sigma - diag_s
  G <- rbind(cbind(G_11, G_12), cbind(G_12, G_11))
  X_tilde <- matrix(0, n, p)
  for (i in 1:n) {
    mu <- X[i, ] - X[i, ] %*% solve(Sigma) %*% diag_s
    V <- 2 * diag_s - diag_s %*% solve(Sigma) %*% diag_s
    
    X_tilde[i, ] <- mvrnorm(1, mu = mu, Sigma = V)
  }
  
  X_large <- cbind(X, X_tilde)
  
  large_fit_ridge_beta1 <- cv.glmnet(X_large, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 0)
  hbeta_large_ridge_beta1 <- coef(large_fit_ridge_beta1, s='lambda.min')[2:(2*500 + 1)]
  w <- abs(hbeta_large_ridge_beta1[1:500]) - abs(hbeta_large_ridge_beta1[501:1000])
  
  kn <- function(t) {
    ((1 + sum(w <= -t))/max(sum(w >= t), 1))
  }
  
  w_kn <- sapply(c(0, abs(w)), kn)
  thres <- min(c(0, abs(w))[w_kn <= 0.2])
  hbeta_knridge_beta1 <- hbeta_ridge_2_beta1 * (w >= thres)
  
  large_fit_ridge_beta2 <- cv.glmnet(X_large, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 0)
  hbeta_large_ridge_beta2 <- coef(large_fit_ridge_beta2, s='lambda.min')[2:(2*500 + 1)]
  w <- abs(hbeta_large_ridge_beta2[1:500]) - abs(hbeta_large_ridge_beta2[501:1000])
  w_kn <- sapply(c(0, abs(w)), kn)
  thres <- min(c(0, abs(w))[w_kn <= 0.2])
  hbeta_knridge_beta2 <- hbeta_ridge_2_beta2 * (w >= thres)
  
  large_fit_lasso_beta1 <- cv.glmnet(X_large, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbeta_large_lasso_beta1 <- coef(large_fit_lasso_beta1, s = 'lambda.min')[2:1001]
  w <- abs(hbeta_large_lasso_beta1[1:500]) - abs(hbeta_large_lasso_beta1[501:1000])
  w_kn_lasso <- sapply(c(0, abs(w)), kn)
  thres_lasso <- min(c(0, abs(w))[w_kn <= 0.2])
  hbeta_knlasso_beta1 <- hbeta_lasso_2_beta1 * (w >= thres)
  large_fit_lasso_beta2 <- cv.glmnet(X_large, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
  hbeta_large_lasso_beta2 <- coef(large_fit_lasso_beta2, s = 'lambda.min')[2:1001]
  w <- abs(hbeta_large_lasso_beta2[1:500]) - abs(hbeta_large_lasso_beta2[501:1000])
  w_kn_lasso <- sapply(c(0, abs(w)), kn)
  thres_lasso <- min(c(0, abs(w))[w_kn <= 0.2])
  hbeta_knlasso_beta2 <- hbeta_lasso_2_beta2 * (w >= thres)
  
  SE_ridge_beta1 <- sum((hbeta_ridge_2_beta1 - beta[[name_k]][["beta1"]])^2)
  SE_ridge_beta2 <- sum((hbeta_ridge_2_beta2 - beta[[name_k]][["beta2"]])^2)
  SE_lasso_beta1 <- sum((hbeta_lasso_2_beta1 - beta[[name_k]][["beta1"]])^2)
  SE_lasso_beta2 <- sum((hbeta_lasso_2_beta2 - beta[[name_k]][["beta2"]])^2)
  SE_adlasso1_beta1 <- sum((hbeta_adlasso1_beta1 - beta[[name_k]][["beta1"]])^2)
  SE_adlasso1_beta2 <- sum((hbeta_adlasso1_beta2 - beta[[name_k]][["beta2"]])^2)
  SE_adlasso2_beta1 <- sum((hbeta_adlasso2_beta1 - beta[[name_k]][["beta1"]])^2)
  SE_adlasso2_beta2 <- sum((hbeta_adlasso2_beta2 - beta[[name_k]][["beta2"]])^2)
  
  true_beta <- beta[[name_k]][["beta1"]] != 0
  k <- sum(true_beta)
  
  lasso_selected_beta1 <- hbeta_lasso_2_beta1 != 0
  TP_lasso_beta1 <- sum(true_beta & lasso_selected_beta1)
  FD_lasso_beta1 <- sum(!true_beta & lasso_selected_beta1)
  R_lasso_beta1 <- sum(lasso_selected_beta1)
  power_lasso_beta1 <- TP_lasso_beta1 / k
  FDR_lasso_beta1 <- if (R_lasso_beta1 == 0) 0 else FD_lasso_beta1 / R_lasso_beta1
  lasso_selected_beta2 <- hbeta_lasso_2_beta2 != 0
  TP_lasso_beta2 <- sum(true_beta & lasso_selected_beta2)
  FD_lasso_beta2 <- sum(!true_beta & lasso_selected_beta2)
  R_lasso_beta2 <- sum(lasso_selected_beta2)
  power_lasso_beta2 <- TP_lasso_beta2 / k
  FDR_lasso_beta2 <- if (R_lasso_beta2 == 0) 0 else FD_lasso_beta2 / R_lasso_beta2
  
  adlasso1_selected_beta1 <- hbeta_adlasso1_beta1 != 0
  TP_adlasso1_beta1 <- sum(true_beta & adlasso1_selected_beta1)
  FD_adlasso1_beta1 <- sum(!true_beta & adlasso1_selected_beta1)
  R_adlasso1_beta1 <- sum(adlasso1_selected_beta1)
  power_adlasso1_beta1 <- TP_adlasso1_beta1 / k
  FDR_adlasso1_beta1 <- if (R_adlasso1_beta1 == 0) 0 else FD_adlasso1_beta1 / R_adlasso1_beta1
  adlasso1_selected_beta2 <- hbeta_adlasso1_beta2 != 0
  TP_adlasso1_beta2 <- sum(true_beta & adlasso1_selected_beta2)
  FD_adlasso1_beta2 <- sum(!true_beta & adlasso1_selected_beta2)
  R_adlasso1_beta2 <- sum(adlasso1_selected_beta2)
  power_adlasso1_beta2 <- TP_adlasso1_beta2 / k
  FDR_adlasso1_beta2 <- if (R_adlasso1_beta2 == 0) 0 else FD_adlasso1_beta2 / R_adlasso1_beta2
  
  adlasso2_selected_beta1 <- hbeta_adlasso2_beta1 != 0
  TP_adlasso2_beta1 <- sum(true_beta & adlasso2_selected_beta1)
  FD_adlasso2_beta1 <- sum(!true_beta & adlasso2_selected_beta1)
  R_adlasso2_beta1 <- sum(adlasso2_selected_beta1)
  power_adlasso2_beta1 <- TP_adlasso2_beta1 / k
  FDR_adlasso2_beta1 <- if (R_adlasso2_beta1 == 0) 0 else FD_adlasso2_beta1 / R_adlasso2_beta1
  adlasso2_selected_beta2 <- hbeta_adlasso2_beta2 != 0
  TP_adlasso2_beta2 <- sum(true_beta & adlasso2_selected_beta2)
  FD_adlasso2_beta2 <- sum(!true_beta & adlasso2_selected_beta2)
  R_adlasso2_beta2 <- sum(adlasso2_selected_beta2)
  power_adlasso2_beta2 <- TP_adlasso2_beta2 / k
  FDR_adlasso2_beta2 <- if (R_adlasso2_beta2 == 0) 0 else FD_adlasso2_beta2 / R_adlasso2_beta2
  
  ridgekn_selected_beta1 <- hbeta_knridge_beta1 != 0
  TP_ridgekn_beta1 <- sum(true_beta & ridgekn_selected_beta1)
  FD_ridgekn_beta1 <- sum(!true_beta & ridgekn_selected_beta1)
  R_ridgekn_beta1 <- sum(ridgekn_selected_beta1)
  power_ridgekn_beta1 <- TP_ridgekn_beta1 / k
  FDR_ridgekn_beta1 <- if (R_ridgekn_beta1 == 0) 0 else FD_ridgekn_beta1 / R_ridgekn_beta1
  ridgekn_selected_beta2 <- hbeta_knridge_beta2 != 0
  TP_ridgekn_beta2 <- sum(true_beta & ridgekn_selected_beta2)
  FD_ridgekn_beta2 <- sum(!true_beta & ridgekn_selected_beta2)
  R_ridgekn_beta2 <- sum(ridgekn_selected_beta2)
  power_ridgekn_beta2 <- TP_ridgekn_beta2 / k
  FDR_ridgekn_beta2 <- if (R_ridgekn_beta2 == 0) 0 else FD_ridgekn_beta2 / R_ridgekn_beta2
  
  lassokn_selected_beta1 <- hbeta_knlasso_beta1 != 0
  TP_lassokn_beta1 <- sum(true_beta & lassokn_selected_beta1)
  FD_lassokn_beta1 <- sum(!true_beta & lassokn_selected_beta1)
  R_lassokn_beta1 <- sum(lassokn_selected_beta1)
  power_lassokn_beta1 <- TP_lassokn_beta1 / k
  FDR_lassokn_beta1 <- if (R_lassokn_beta1 == 0) 0 else FD_lassokn_beta1 / R_lassokn_beta1
  lassokn_selected_beta2 <- hbeta_knlasso_beta2 != 0
  TP_lassokn_beta2 <- sum(true_beta & lassokn_selected_beta2)
  FD_lassokn_beta2 <- sum(!true_beta & lassokn_selected_beta2)
  R_lassokn_beta2 <- sum(lassokn_selected_beta2)
  power_lassokn_beta2 <- TP_lassokn_beta2 / k
  FDR_lassokn_beta2 <- if (R_lassokn_beta2 == 0) 0 else FD_lassokn_beta2 / R_lassokn_beta2
  
  return(list(
    hbeta_ridge_beta1 = hbeta_ridge_2_beta1,
    SE_ridge_beta1 = SE_ridge_beta1,
    hbeta_ridge_beta2 = hbeta_ridge_2_beta2,
    SE_ridge_beta2 = SE_ridge_beta2,
    hbeta_lasso_beta1 = hbeta_lasso_2_beta1,
    SE_lasso_beta1 = SE_lasso_beta1,
    power_lasso_beta1 = power_lasso_beta1,
    FDR_lasso_beta1 = FDR_lasso_beta1,
    hbeta_lasso_beta2 = hbeta_lasso_2_beta2,
    SE_lasso_beta2 = SE_lasso_beta2,
    power_lasso_beta2 = power_lasso_beta2,
    FDR_lasso_beta2 = FDR_lasso_beta2,
    hbeta_adlasso1_beta1 = hbeta_adlasso1_beta1,
    SE_adlasso1_beta1 = SE_adlasso1_beta1,
    power_adlasso1_beta1 = power_adlasso1_beta1,
    FDR_adlasso1_beta1 = FDR_adlasso1_beta1,
    hbeta_adlasso1_beta2 = hbeta_adlasso1_beta2,
    SE_adlasso1_beta2 = SE_adlasso1_beta2,
    power_adlasso1_beta2 = power_adlasso1_beta2,
    FDR_adlasso1_beta2 = FDR_adlasso1_beta2,
    hbeta_adlasso2_beta1 = hbeta_adlasso2_beta1,
    SE_adlasso2_beta1 = SE_adlasso2_beta1,
    power_adlasso2_beta1 = power_adlasso2_beta1,
    FDR_adlasso2_beta1 = FDR_adlasso2_beta1,
    hbeta_adlasso2_beta2 = hbeta_adlasso2_beta2,
    SE_adlasso2_beta2 = SE_adlasso2_beta2,
    power_adlasso2_beta2 = power_adlasso2_beta2,
    FDR_adlasso2_beta2 = FDR_adlasso2_beta2,
    hbeta_knridge_beta1 = hbeta_knridge_beta1,
    power_ridgekn_beta1 = power_ridgekn_beta1,
    FDR_ridgekn_beta1 = FDR_ridgekn_beta1,
    hbeta_knridge_beta2 = hbeta_knridge_beta2,
    power_ridgekn_beta2 = power_ridgekn_beta2,
    FDR_ridgekn_beta2 = FDR_ridgekn_beta2,
    hbeta_knlasso_beta1 = hbeta_knlasso_beta1,
    power_lassokn_beta1 = power_lassokn_beta1,
    FDR_lassokn_beta1 = FDR_lassokn_beta1,
    hbeta_knlasso_beta2 = hbeta_knlasso_beta2,
    power_lassokn_beta2 = power_lassokn_beta2,
    FDR_lassokn_beta2 = FDR_lassokn_beta2)
  )
}), names(beta))

saveRDS(data3, "data3.rds")

SE_beta1 <- data.frame(
  SE_ridge_beta1 <- c(data3[["k5"]][["SE_ridge_beta1"]], data3[["k20"]][["SE_ridge_beta1"]], data3[["k100"]][["SE_ridge_beta1"]]),
  SE_lasso_beta1 <- c(data3[["k5"]][["SE_lasso_beta1"]], data3[["k20"]][["SE_lasso_beta1"]], data3[["k100"]][["SE_lasso_beta1"]]),
  SE_adlasso1_beta1 <- c(data3[["k5"]][["SE_adlasso1_beta1"]], data3[["k20"]][["SE_adlasso1_beta1"]], data3[["k100"]][["SE_adlasso1_beta1"]]),
  SE_adlasso2_beta1 <- c(data3[["k5"]][["SE_adlasso2_beta1"]], data3[["k20"]][["SE_adlasso2_beta1"]], data3[["k100"]][["SE_adlasso2_beta1"]])
)
colnames(SE_beta1) <- c("Ridge", "LASSO", "AdLASSO1", "AdLASSO2")
rownames(SE_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(SE_beta1)

SE_beta2 <- data.frame(
  SE_ridge_beta2 <- c(data3[["k5"]][["SE_ridge_beta2"]], data3[["k20"]][["SE_ridge_beta2"]], data3[["k100"]][["SE_ridge_beta2"]]),
  SE_lasso_beta2 <- c(data3[["k5"]][["SE_lasso_beta2"]], data3[["k20"]][["SE_lasso_beta2"]], data3[["k100"]][["SE_lasso_beta2"]]),
  SE_adlasso1_beta2 <- c(data3[["k5"]][["SE_adlasso1_beta2"]], data3[["k20"]][["SE_adlasso1_beta2"]], data3[["k100"]][["SE_adlasso1_beta2"]]),
  SE_adlasso2_beta2 <- c(data3[["k5"]][["SE_adlasso2_beta2"]], data3[["k20"]][["SE_adlasso2_beta2"]], data3[["k100"]][["SE_adlasso2_beta2"]])
)
colnames(SE_beta2) <- c("Ridge", "LASSO", "AdLASSO1", "AdLASSO2")
rownames(SE_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(SE_beta2)

power_beta1 <- data.frame(
  power_lasso_beta1 <- c(data3[["k5"]][["power_lasso_beta1"]], data3[["k20"]][["power_lasso_beta1"]], data3[["k100"]][["power_lasso_beta1"]]),
  power_adlasso1_beta1 <- c(data3[["k5"]][["power_adlasso1_beta1"]], data3[["k20"]][["power_adlasso1_beta1"]], data3[["k100"]][["power_adlasso1_beta1"]]),
  power_adlasso2_beta1 <- c(data3[["k5"]][["power_adlasso2_beta1"]], data3[["k20"]][["power_adlasso2_beta1"]], data3[["k100"]][["power_adlasso2_beta1"]]),
  power_ridgekn_beta1 <- c(data3[["k5"]][["power_ridgekn_beta1"]], data3[["k20"]][["power_ridgekn_beta1"]], data3[["k100"]][["power_ridgekn_beta1"]]),
  power_lassokn_beta1 <- c(data3[["k5"]][["power_lassokn_beta1"]], data3[["k20"]][["power_lassokn_beta1"]], data3[["k100"]][["power_lassokn_beta1"]])
) 

colnames(power_beta1) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(power_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(power_beta1)

power_beta2 <- data.frame(
  power_lasso_beta2 <- c(data3[["k5"]][["power_lasso_beta2"]], data3[["k20"]][["power_lasso_beta2"]], data3[["k100"]][["power_lasso_beta2"]]),
  power_adlasso1_beta2 <- c(data3[["k5"]][["power_adlasso1_beta2"]], data3[["k20"]][["power_adlasso1_beta2"]], data3[["k100"]][["power_adlasso1_beta2"]]),
  power_adlasso2_beta2 <- c(data3[["k5"]][["power_adlasso2_beta2"]], data3[["k20"]][["power_adlasso2_beta2"]], data3[["k100"]][["power_adlasso2_beta2"]]),
  power_ridgekn_beta2 <- c(data3[["k5"]][["power_ridgekn_beta2"]], data3[["k20"]][["power_ridgekn_beta2"]], data3[["k100"]][["power_ridgekn_beta2"]]),
  power_lassokn_beta2 <- c(data3[["k5"]][["power_lassokn_beta2"]], data3[["k20"]][["power_lassokn_beta2"]], data3[["k100"]][["power_lassokn_beta2"]])
) 

colnames(power_beta2) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(power_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(power_beta2)

FDR_beta1 <- data.frame(
  FDR_lasso_beta1 <- c(data3[["k5"]][["FDR_lasso_beta1"]], data3[["k20"]][["FDR_lasso_beta1"]], data3[["k100"]][["FDR_lasso_beta1"]]),
  FDR_adlasso1_beta1 <- c(data3[["k5"]][["FDR_adlasso1_beta1"]], data3[["k20"]][["FDR_adlasso1_beta1"]], data3[["k100"]][["FDR_adlasso1_beta1"]]),
  FDR_adlasso2_beta1 <- c(data3[["k5"]][["FDR_adlasso2_beta1"]], data3[["k20"]][["FDR_adlasso2_beta1"]], data3[["k100"]][["FDR_adlasso2_beta1"]]),
  FDR_ridgekn_beta1 <- c(data3[["k5"]][["FDR_ridgekn_beta1"]], data3[["k20"]][["FDR_ridgekn_beta1"]], data3[["k100"]][["FDR_ridgekn_beta1"]]),
  FDR_lassokn_beta1 <- c(data3[["k5"]][["FDR_lassokn_beta1"]], data3[["k20"]][["FDR_lassokn_beta1"]], data3[["k100"]][["FDR_lassokn_beta1"]])
)

colnames(FDR_beta1) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(FDR_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(FDR_beta1)

FDR_beta2 <- data.frame(
  FDR_lasso_beta2 <- c(data3[["k5"]][["FDR_lasso_beta2"]], data3[["k20"]][["FDR_lasso_beta2"]], data3[["k100"]][["FDR_lasso_beta2"]]),
  FDR_adlasso1_beta2 <- c(data3[["k5"]][["FDR_adlasso1_beta2"]], data3[["k20"]][["FDR_adlasso1_beta2"]], data3[["k100"]][["FDR_adlasso1_beta2"]]),
  FDR_adlasso2_beta2 <- c(data3[["k5"]][["FDR_adlasso2_beta2"]], data3[["k20"]][["FDR_adlasso2_beta2"]], data3[["k100"]][["FDR_adlasso2_beta2"]]),
  FDR_ridgekn_beta2 <- c(data3[["k5"]][["FDR_ridgekn_beta2"]], data3[["k20"]][["FDR_ridgekn_beta2"]], data3[["k100"]][["FDR_ridgekn_beta2"]]),
  FDR_lassokn_beta2 <- c(data3[["k5"]][["FDR_lassokn_beta2"]], data3[["k20"]][["FDR_lassokn_beta2"]], data3[["k100"]][["FDR_lassokn_beta2"]])
)

colnames(FDR_beta2) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(FDR_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(FDR_beta2)

rep <- 100

results2 <- replicate(rep, {
  Y <- lapply(beta, function(b) {
    Y_beta1 <- X%*%b[["beta1"]] + rnorm(500, 0, 1)
    Y_beta2 <- X%*%b[["beta2"]] + rnorm(500, 0, 1)
    return(list(Y_beta1 = Y_beta1, Y_beta2 = Y_beta2))
  })
  
  data3 <- setNames(lapply(names(Y), function(name_k) {
    y <- Y[[name_k]]
    
    ridge_fit_2_beta1 <- cv.glmnet(X, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 0)
    hbeta_ridge_2_beta1 <- coef(ridge_fit_2_beta1, s = 'lambda.min')[2:501]
    ridge_fit_2_beta2 <- cv.glmnet(X, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 0)
    hbeta_ridge_2_beta2 <- coef(ridge_fit_2_beta2, s = 'lambda.min')[2:501]
    
    lasso_fit_2_beta1 <- cv.glmnet(X, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbeta_lasso_2_beta1 <- coef(lasso_fit_2_beta1, s = 'lambda.min')[2:501]
    lasso_fit_2_beta2 <- cv.glmnet(X, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbeta_lasso_2_beta2 <- coef(lasso_fit_2_beta2, s = 'lambda.min')[2:501]
    
    W1_beta1 <- 1/abs(hbeta_ridge_2_beta1)
    X_temp_beta1 <- X %*% diag(1/W1_beta1)
    adlasso1_fit_beta1 <- cv.glmnet(X_temp_beta1, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbetatemp_adlasso1_beta1 <- coef(adlasso1_fit_beta1, s = 'lambda.min')[2:501]
    hbeta_adlasso1_beta1 <- hbetatemp_adlasso1_beta1/W1_beta1
    W1_beta2 <- 1/(abs(hbeta_ridge_2_beta2))
    X_temp_beta2 <- X %*% diag(1/W1_beta2)
    adlasso1_fit_beta2 <- cv.glmnet(X_temp_beta2, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbetatemp_adlasso1_beta2 <- coef(adlasso1_fit_beta2, s = 'lambda.min')[2:501]
    hbeta_adlasso1_beta2 <- hbetatemp_adlasso1_beta2/W1_beta2
    
    RSS_beta1 <- sum((y[["Y_beta1"]] - X %*% hbeta_lasso_2_beta1)^2)
    ind_lasso_beta1 <- which(hbeta_lasso_2_beta1 > 0)
    if (length(ind_lasso_beta1) == 0) {
      hbeta_adlasso2_beta1 = rep(0, 500)
    } else {
      hsigma_beta1 <- sqrt(RSS_beta1/(500 - length(ind_lasso_beta1)))
      W2_beta1 <- hsigma_beta1/abs(hbeta_lasso_2_beta1[ind_lasso_beta1])
      lambda1_beta1 <- hsigma_beta1 * qnorm(1 - 0.2/(2 * 500))
      X_small_beta1 <- X[, ind_lasso_beta1]
      if (length(ind_lasso_beta1) == 1) {
        X_temp_beta1 <- matrix(X_small_beta1 / W2_beta1, ncol = 1)
        adlasso2_fit_beta1 <- lm(y[["Y_beta1"]] ~ X_temp_beta1 + 0)
        hbetatemp_adlasso2_beta1 <- coef(adlasso2_fit_beta1)
      } else {
        X_temp_beta1 <- X_small_beta1 %*% diag(1 / W2_beta1)
        adlasso2_fit_beta1 <- glmnet(X_temp_beta1, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1, lambda = lambda1_beta1/500)
        hbetatemp_adlasso2_beta1 <- coef(adlasso2_fit_beta1)[2:(length(ind_lasso_beta1) + 1)]
      }
      hbeta_adlasso2_beta1 <- rep(0, 500)
      hbeta_adlasso2_beta1[ind_lasso_beta1] <- hbetatemp_adlasso2_beta1/W2_beta1
    }
    RSS_beta2 <- sum((y[["Y_beta2"]] - X %*% hbeta_lasso_2_beta2)^2)
    ind_lasso_beta2 <- which(hbeta_lasso_2_beta2 > 0)
    if (length(ind_lasso_beta2) == 0) {
      hbeta_adlasso2_beta2 = rep(0, 500)
    } else {
      hsigma_beta2 <- sqrt(RSS_beta2/(500 - length(ind_lasso_beta2)))
      W2_beta2 <- hsigma_beta2/abs(hbeta_lasso_2_beta2[ind_lasso_beta2])
      lambda1_beta2 <- hsigma_beta2 * qnorm(1 - 0.2/(2 * 500))
      X_small_beta2 <- X[, ind_lasso_beta2]
      if (length(ind_lasso_beta2) == 1) {
        X_temp_beta2 <- matrix(X_small_beta2 / W2_beta2, ncol = 1)
        adlasso2_fit_beta2 <- lm(y[["Y_beta2"]] ~ X_temp_beta2 + 0)
        hbetatemp_adlasso2_beta2 <- coef(adlasso2_fit_beta2)
      } else {
        X_temp_beta2 <- X_small_beta2 %*% diag(1 / W2_beta2)
        adlasso2_fit_beta2 <- glmnet(X_temp_beta2, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1, lambda = lambda1_beta2/500)
        hbetatemp_adlasso2_beta2 <- coef(adlasso2_fit_beta2)[2:(length(ind_lasso_beta2) +1)]
      }
      
      hbeta_adlasso2_beta2 <- rep(0, 500)
      hbeta_adlasso2_beta2[ind_lasso_beta2] <- hbetatemp_adlasso2_beta2/W2_beta2
    }
    lambda_min <- min(eigen(Sigma)$values)
    s <- rep(2 * lambda_min, p)
    s <- pmin(s, 1)  
    diag_s <- diag(s)
    
    G_11 <- Sigma
    G_12 <- Sigma - diag_s
    G <- rbind(cbind(G_11, G_12), cbind(G_12, G_11))
    X_tilde <- matrix(0, n, p)
    for (i in 1:n) {
      mu <- X[i, ] - X[i, ] %*% solve(Sigma) %*% diag_s
      V <- 2 * diag_s - diag_s %*% solve(Sigma) %*% diag_s
      
      X_tilde[i, ] <- mvrnorm(1, mu = mu, Sigma = V)
    }
    
    X_large <- cbind(X, X_tilde)
    large_fit_ridge_beta1 <- cv.glmnet(X_large, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 0)
    hbeta_large_ridge_beta1 <- coef(large_fit_ridge_beta1, s='lambda.min')[2:(2*500 + 1)]
    w <- abs(hbeta_large_ridge_beta1[1:500]) - abs(hbeta_large_ridge_beta1[501:1000])
    
    kn <- function(t) {
      ((1 + sum(w <= -t))/max(sum(w >= t), 1))
    }
    
    w_kn <- sapply(c(0, abs(w)), kn)
    thres <- min(c(0, abs(w))[w_kn <= 0.2])
    hbeta_knridge_beta1 <- hbeta_ridge_2_beta1 * (w >= thres)
    
    large_fit_ridge_beta2 <- cv.glmnet(X_large, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 0)
    hbeta_large_ridge_beta2 <- coef(large_fit_ridge_beta2, s='lambda.min')[2:(2*500 + 1)]
    w <- abs(hbeta_large_ridge_beta2[1:500]) - abs(hbeta_large_ridge_beta2[501:1000])
    w_kn <- sapply(c(0, abs(w)), kn)
    thres <- min(c(0, abs(w))[w_kn <= 0.2])
    hbeta_knridge_beta2 <- hbeta_ridge_2_beta2 * (w >= thres)
    
    large_fit_lasso_beta1 <- cv.glmnet(X_large, y[["Y_beta1"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbeta_large_lasso_beta1 <- coef(large_fit_lasso_beta1, s = 'lambda.min')[2:1001]
    w <- abs(hbeta_large_lasso_beta1[1:500]) - abs(hbeta_large_lasso_beta1[501:1000])
    w_kn_lasso <- sapply(c(0, abs(w)), kn)
    thres_lasso <- min(c(0, abs(w))[w_kn <= 0.2])
    hbeta_knlasso_beta1 <- hbeta_lasso_2_beta1 * (w >= thres)
    large_fit_lasso_beta2 <- cv.glmnet(X_large, y[["Y_beta2"]], standardize = FALSE, intercept = FALSE, alpha = 1)
    hbeta_large_lasso_beta2 <- coef(large_fit_lasso_beta2, s = 'lambda.min')[2:1001]
    w <- abs(hbeta_large_lasso_beta2[1:500]) - abs(hbeta_large_lasso_beta2[501:1000])
    w_kn_lasso <- sapply(c(0, abs(w)), kn)
    thres_lasso <- min(c(0, abs(w))[w_kn <= 0.2])
    hbeta_knlasso_beta2 <- hbeta_lasso_2_beta2 * (w >= thres)
    
    SE_ridge_beta1 <- sum((hbeta_ridge_2_beta1 - beta[[name_k]][["beta1"]])^2)
    SE_ridge_beta2 <- sum((hbeta_ridge_2_beta2 - beta[[name_k]][["beta2"]])^2)
    SE_lasso_beta1 <- sum((hbeta_lasso_2_beta1 - beta[[name_k]][["beta1"]])^2)
    SE_lasso_beta2 <- sum((hbeta_lasso_2_beta2 - beta[[name_k]][["beta2"]])^2)
    SE_adlasso1_beta1 <- sum((hbeta_adlasso1_beta1 - beta[[name_k]][["beta1"]])^2)
    SE_adlasso1_beta2 <- sum((hbeta_adlasso1_beta2 - beta[[name_k]][["beta2"]])^2)
    SE_adlasso2_beta1 <- sum((hbeta_adlasso2_beta1 - beta[[name_k]][["beta1"]])^2)
    SE_adlasso2_beta2 <- sum((hbeta_adlasso2_beta2 - beta[[name_k]][["beta2"]])^2)
    
    true_beta <- beta[[name_k]][["beta1"]] != 0
    k <- sum(true_beta)
    
    lasso_selected_beta1 <- hbeta_lasso_2_beta1 != 0
    TP_lasso_beta1 <- sum(true_beta & lasso_selected_beta1)
    FD_lasso_beta1 <- sum(!true_beta & lasso_selected_beta1)
    R_lasso_beta1 <- sum(lasso_selected_beta1)
    power_lasso_beta1 <- TP_lasso_beta1 / k
    FDR_lasso_beta1 <- if (R_lasso_beta1 == 0) 0 else FD_lasso_beta1 / R_lasso_beta1
    lasso_selected_beta2 <- hbeta_lasso_2_beta2 != 0
    TP_lasso_beta2 <- sum(true_beta & lasso_selected_beta2)
    FD_lasso_beta2 <- sum(!true_beta & lasso_selected_beta2)
    R_lasso_beta2 <- sum(lasso_selected_beta2)
    power_lasso_beta2 <- TP_lasso_beta2 / k
    FDR_lasso_beta2 <- if (R_lasso_beta2 == 0) 0 else FD_lasso_beta2 / R_lasso_beta2
    
    adlasso1_selected_beta1 <- hbeta_adlasso1_beta1 != 0
    TP_adlasso1_beta1 <- sum(true_beta & adlasso1_selected_beta1)
    FD_adlasso1_beta1 <- sum(!true_beta & adlasso1_selected_beta1)
    R_adlasso1_beta1 <- sum(adlasso1_selected_beta1)
    power_adlasso1_beta1 <- TP_adlasso1_beta1 / k
    FDR_adlasso1_beta1 <- if (R_adlasso1_beta1 == 0) 0 else FD_adlasso1_beta1 / R_adlasso1_beta1
    adlasso1_selected_beta2 <- hbeta_adlasso1_beta2 != 0
    TP_adlasso1_beta2 <- sum(true_beta & adlasso1_selected_beta2)
    FD_adlasso1_beta2 <- sum(!true_beta & adlasso1_selected_beta2)
    R_adlasso1_beta2 <- sum(adlasso1_selected_beta2)
    power_adlasso1_beta2 <- TP_adlasso1_beta2 / k
    FDR_adlasso1_beta2 <- if (R_adlasso1_beta2 == 0) 0 else FD_adlasso1_beta2 / R_adlasso1_beta2
    
    adlasso2_selected_beta1 <- hbeta_adlasso2_beta1 != 0
    TP_adlasso2_beta1 <- sum(true_beta & adlasso2_selected_beta1)
    FD_adlasso2_beta1 <- sum(!true_beta & adlasso2_selected_beta1)
    R_adlasso2_beta1 <- sum(adlasso2_selected_beta1)
    power_adlasso2_beta1 <- TP_adlasso2_beta1 / k
    FDR_adlasso2_beta1 <- if (R_adlasso2_beta1 == 0) 0 else FD_adlasso2_beta1 / R_adlasso2_beta1
    adlasso2_selected_beta2 <- hbeta_adlasso2_beta2 != 0
    TP_adlasso2_beta2 <- sum(true_beta & adlasso2_selected_beta2)
    FD_adlasso2_beta2 <- sum(!true_beta & adlasso2_selected_beta2)
    R_adlasso2_beta2 <- sum(adlasso2_selected_beta2)
    power_adlasso2_beta2 <- TP_adlasso2_beta2 / k
    FDR_adlasso2_beta2 <- if (R_adlasso2_beta2 == 0) 0 else FD_adlasso2_beta2 / R_adlasso2_beta2
    
    ridgekn_selected_beta1 <- hbeta_knridge_beta1 != 0
    TP_ridgekn_beta1 <- sum(true_beta & ridgekn_selected_beta1)
    FD_ridgekn_beta1 <- sum(!true_beta & ridgekn_selected_beta1)
    R_ridgekn_beta1 <- sum(ridgekn_selected_beta1)
    power_ridgekn_beta1 <- TP_ridgekn_beta1 / k
    FDR_ridgekn_beta1 <- if (R_ridgekn_beta1 == 0) 0 else FD_ridgekn_beta1 / R_ridgekn_beta1
    ridgekn_selected_beta2 <- hbeta_knridge_beta2 != 0
    TP_ridgekn_beta2 <- sum(true_beta & ridgekn_selected_beta2)
    FD_ridgekn_beta2 <- sum(!true_beta & ridgekn_selected_beta2)
    R_ridgekn_beta2 <- sum(ridgekn_selected_beta2)
    power_ridgekn_beta2 <- TP_ridgekn_beta2 / k
    FDR_ridgekn_beta2 <- if (R_ridgekn_beta2 == 0) 0 else FD_ridgekn_beta2 / R_ridgekn_beta2
    
    lassokn_selected_beta1 <- hbeta_knlasso_beta1 != 0
    TP_lassokn_beta1 <- sum(true_beta & lassokn_selected_beta1)
    FD_lassokn_beta1 <- sum(!true_beta & lassokn_selected_beta1)
    R_lassokn_beta1 <- sum(lassokn_selected_beta1)
    power_lassokn_beta1 <- TP_lassokn_beta1 / k
    FDR_lassokn_beta1 <- if (R_lassokn_beta1 == 0) 0 else FD_lassokn_beta1 / R_lassokn_beta1
    lassokn_selected_beta2 <- hbeta_knlasso_beta2 != 0
    TP_lassokn_beta2 <- sum(true_beta & lassokn_selected_beta2)
    FD_lassokn_beta2 <- sum(!true_beta & lassokn_selected_beta2)
    R_lassokn_beta2 <- sum(lassokn_selected_beta2)
    power_lassokn_beta2 <- TP_lassokn_beta2 / k
    FDR_lassokn_beta2 <- if (R_lassokn_beta2 == 0) 0 else FD_lassokn_beta2 / R_lassokn_beta2
    
    return(list(
      hbeta_ridge_beta1 = hbeta_ridge_2_beta1,
      SE_ridge_beta1 = SE_ridge_beta1,
      hbeta_ridge_beta2 = hbeta_ridge_2_beta2,
      SE_ridge_beta2 = SE_ridge_beta2,
      hbeta_lasso_beta1 = hbeta_lasso_2_beta1,
      SE_lasso_beta1 = SE_lasso_beta1,
      power_lasso_beta1 = power_lasso_beta1,
      FDR_lasso_beta1 = FDR_lasso_beta1,
      hbeta_lasso_beta2 = hbeta_lasso_2_beta2,
      SE_lasso_beta2 = SE_lasso_beta2,
      power_lasso_beta2 = power_lasso_beta2,
      FDR_lasso_beta2 = FDR_lasso_beta2,
      hbeta_adlasso1_beta1 = hbeta_adlasso1_beta1,
      SE_adlasso1_beta1 = SE_adlasso1_beta1,
      power_adlasso1_beta1 = power_adlasso1_beta1,
      FDR_adlasso1_beta1 = FDR_adlasso1_beta1,
      hbeta_adlasso1_beta2 = hbeta_adlasso1_beta2,
      SE_adlasso1_beta2 = SE_adlasso1_beta2,
      power_adlasso1_beta2 = power_adlasso1_beta2,
      FDR_adlasso1_beta2 = FDR_adlasso1_beta2,
      hbeta_adlasso2_beta1 = hbeta_adlasso2_beta1,
      SE_adlasso2_beta1 = SE_adlasso2_beta1,
      power_adlasso2_beta1 = power_adlasso2_beta1,
      FDR_adlasso2_beta1 = FDR_adlasso2_beta1,
      hbeta_adlasso2_beta2 = hbeta_adlasso2_beta2,
      SE_adlasso2_beta2 = SE_adlasso2_beta2,
      power_adlasso2_beta2 = power_adlasso2_beta2,
      FDR_adlasso2_beta2 = FDR_adlasso2_beta2,
      hbeta_knridge_beta1 = hbeta_knridge_beta1,
      power_ridgekn_beta1 = power_ridgekn_beta1,
      FDR_ridgekn_beta1 = FDR_ridgekn_beta1,
      hbeta_knridge_beta2 = hbeta_knridge_beta2,
      power_ridgekn_beta2 = power_ridgekn_beta2,
      FDR_ridgekn_beta2 = FDR_ridgekn_beta2,
      hbeta_knlasso_beta1 = hbeta_knlasso_beta1,
      power_lassokn_beta1 = power_lassokn_beta1,
      FDR_lassokn_beta1 = FDR_lassokn_beta1,
      hbeta_knlasso_beta2 = hbeta_knlasso_beta2,
      power_lassokn_beta2 = power_lassokn_beta2,
      FDR_lassokn_beta2 = FDR_lassokn_beta2)
    )
  }), names(beta))
  
  return(data3)
}, simplify = FALSE)

saveRDS(results2, "results2.rds")

library(dplyr)

averaged_metrics_by_k <- extract_metrics(results2)

MSE_beta1 <- data.frame(
  SE_ridge_beta1 <- c(averaged_metrics_by_k[["k5"]][["SE_ridge_beta1"]], averaged_metrics_by_k[["k20"]][["SE_ridge_beta1"]], averaged_metrics_by_k[["k100"]][["SE_ridge_beta1"]]),
  SE_lasso_beta1 <- c(averaged_metrics_by_k[["k5"]][["SE_lasso_beta1"]], averaged_metrics_by_k[["k20"]][["SE_lasso_beta1"]], averaged_metrics_by_k[["k100"]][["SE_lasso_beta1"]]),
  SE_adlasso1_beta1 <- c(averaged_metrics_by_k[["k5"]][["SE_adlasso1_beta1"]], averaged_metrics_by_k[["k20"]][["SE_adlasso1_beta1"]], averaged_metrics_by_k[["k100"]][["SE_adlasso1_beta1"]]),
  SE_adlasso2_beta1 <- c(averaged_metrics_by_k[["k5"]][["SE_adlasso2_beta1"]], averaged_metrics_by_k[["k20"]][["SE_adlasso2_beta1"]], averaged_metrics_by_k[["k100"]][["SE_adlasso2_beta1"]])
)
colnames(MSE_beta1) <- c("Ridge", "LASSO", "AdLASSO1", "AdLASSO2")
rownames(MSE_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(MSE_beta1)

MSE_beta2<- data.frame(
  SE_ridge_beta2 <- c(averaged_metrics_by_k[["k5"]][["SE_ridge_beta2"]], averaged_metrics_by_k[["k20"]][["SE_ridge_beta2"]], averaged_metrics_by_k[["k100"]][["SE_ridge_beta2"]]),
  SE_lasso_beta2 <- c(averaged_metrics_by_k[["k5"]][["SE_lasso_beta2"]], averaged_metrics_by_k[["k20"]][["SE_lasso_beta2"]], averaged_metrics_by_k[["k100"]][["SE_lasso_beta2"]]),
  SE_adlasso1_beta2 <- c(averaged_metrics_by_k[["k5"]][["SE_adlasso1_beta2"]], averaged_metrics_by_k[["k20"]][["SE_adlasso1_beta2"]], averaged_metrics_by_k[["k100"]][["SE_adlasso1_beta2"]]),
  SE_adlasso2_beta2 <- c(averaged_metrics_by_k[["k5"]][["SE_adlasso2_beta2"]], averaged_metrics_by_k[["k20"]][["SE_adlasso2_beta2"]], averaged_metrics_by_k[["k100"]][["SE_adlasso2_beta2"]])
)
colnames(MSE_beta2) <- c("Ridge", "LASSO", "AdLASSO1", "AdLASSO2")
rownames(MSE_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(MSE_beta2)

power_beta1 <- data.frame(
  power_lasso_beta1 <- c(averaged_metrics_by_k[["k5"]][["power_lasso_beta1"]], averaged_metrics_by_k[["k20"]][["power_lasso_beta1"]], averaged_metrics_by_k[["k100"]][["power_lasso_beta1"]]),
  power_adlasso1_beta1 <- c(averaged_metrics_by_k[["k5"]][["power_adlasso1_beta1"]], averaged_metrics_by_k[["k20"]][["power_adlasso1_beta1"]], averaged_metrics_by_k[["k100"]][["power_adlasso1_beta1"]]),
  power_adlasso2_beta1 <- c(averaged_metrics_by_k[["k5"]][["power_adlasso2_beta1"]], averaged_metrics_by_k[["k20"]][["power_adlasso2_beta1"]], averaged_metrics_by_k[["k100"]][["power_adlasso2_beta1"]]),
  power_ridgekn_beta1 <- c(averaged_metrics_by_k[["k5"]][["power_ridgekn_beta1"]], averaged_metrics_by_k[["k20"]][["power_ridgekn_beta1"]], averaged_metrics_by_k[["k100"]][["power_ridgekn_beta1"]]),
  power_lassokn_beta1 <- c(averaged_metrics_by_k[["k5"]][["power_lassokn_beta1"]], averaged_metrics_by_k[["k20"]][["power_lassokn_beta1"]], averaged_metrics_by_k[["k100"]][["power_lassokn_beta1"]])
) 

colnames(power_beta1) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(power_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(power_beta1)

power_beta2 <- data.frame(
  power_lasso_beta2 <- c(averaged_metrics_by_k[["k5"]][["power_lasso_beta2"]], averaged_metrics_by_k[["k20"]][["power_lasso_beta2"]], averaged_metrics_by_k[["k100"]][["power_lasso_beta2"]]),
  power_adlasso1_beta2 <- c(averaged_metrics_by_k[["k5"]][["power_adlasso1_beta2"]], averaged_metrics_by_k[["k20"]][["power_adlasso1_beta2"]], averaged_metrics_by_k[["k100"]][["power_adlasso1_beta2"]]),
  power_adlasso2_beta2 <- c(averaged_metrics_by_k[["k5"]][["power_adlasso2_beta2"]], averaged_metrics_by_k[["k20"]][["power_adlasso2_beta2"]], averaged_metrics_by_k[["k100"]][["power_adlasso2_beta2"]]),
  power_ridgekn_beta2 <- c(averaged_metrics_by_k[["k5"]][["power_ridgekn_beta2"]], averaged_metrics_by_k[["k20"]][["power_ridgekn_beta2"]], averaged_metrics_by_k[["k100"]][["power_ridgekn_beta2"]]),
  power_lassokn_beta2 <- c(averaged_metrics_by_k[["k5"]][["power_lassokn_beta2"]], averaged_metrics_by_k[["k20"]][["power_lassokn_beta2"]], averaged_metrics_by_k[["k100"]][["power_lassokn_beta2"]])
) 

colnames(power_beta2) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(power_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(power_beta2)

FDR_beta1 <- data.frame(
  FDR_lasso_beta1 <- c(averaged_metrics_by_k[["k5"]][["FDR_lasso_beta1"]], averaged_metrics_by_k[["k20"]][["FDR_lasso_beta1"]], averaged_metrics_by_k[["k100"]][["FDR_lasso_beta1"]]),
  FDR_adlasso1_beta1 <- c(averaged_metrics_by_k[["k5"]][["FDR_adlasso1_beta1"]], averaged_metrics_by_k[["k20"]][["FDR_adlasso1_beta1"]], averaged_metrics_by_k[["k100"]][["FDR_adlasso1_beta1"]]),
  FDR_adlasso2_beta1 <- c(averaged_metrics_by_k[["k5"]][["FDR_adlasso2_beta1"]], averaged_metrics_by_k[["k20"]][["FDR_adlasso2_beta1"]], averaged_metrics_by_k[["k100"]][["FDR_adlasso2_beta1"]]),
  FDR_ridgekn_beta1 <- c(averaged_metrics_by_k[["k5"]][["FDR_ridgekn_beta1"]], averaged_metrics_by_k[["k20"]][["FDR_ridgekn_beta1"]], averaged_metrics_by_k[["k100"]][["FDR_ridgekn_beta1"]]),
  FDR_lassokn_beta1 <- c(averaged_metrics_by_k[["k5"]][["FDR_lassokn_beta1"]], averaged_metrics_by_k[["k20"]][["FDR_lassokn_beta1"]], averaged_metrics_by_k[["k100"]][["FDR_lassokn_beta1"]])
)

colnames(FDR_beta1) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(FDR_beta1) <- c("k = 5", "k = 20", "k = 100")

kable(FDR_beta1)

FDR_beta2 <- data.frame(
  FDR_lasso_beta2 <- c(averaged_metrics_by_k[["k5"]][["FDR_lasso_beta2"]], averaged_metrics_by_k[["k20"]][["FDR_lasso_beta2"]], averaged_metrics_by_k[["k100"]][["FDR_lasso_beta2"]]),
  FDR_adlasso1_beta2 <- c(averaged_metrics_by_k[["k5"]][["FDR_adlasso1_beta2"]], averaged_metrics_by_k[["k20"]][["FDR_adlasso1_beta2"]], averaged_metrics_by_k[["k100"]][["FDR_adlasso1_beta2"]]),
  FDR_adlasso2_beta2 <- c(averaged_metrics_by_k[["k5"]][["FDR_adlasso2_beta2"]], averaged_metrics_by_k[["k20"]][["FDR_adlasso2_beta2"]], averaged_metrics_by_k[["k100"]][["FDR_adlasso2_beta2"]]),
  FDR_ridgekn_beta2 <- c(averaged_metrics_by_k[["k5"]][["FDR_ridgekn_beta2"]], averaged_metrics_by_k[["k20"]][["FDR_ridgekn_beta2"]], averaged_metrics_by_k[["k100"]][["FDR_ridgekn_beta2"]]),
  FDR_lassokn_beta2 <- c(averaged_metrics_by_k[["k5"]][["FDR_lassokn_beta2"]], averaged_metrics_by_k[["k20"]][["FDR_lassokn_beta2"]], averaged_metrics_by_k[["k100"]][["FDR_lassokn_beta2"]])
)

colnames(FDR_beta2) <- c("LASSO", "AdLASSO1", "AdLASSO2", "Ridge knockoffs", "LASSO knockoffs")
rownames(FDR_beta2) <- c("k = 5", "k = 20", "k = 100")

kable(FDR_beta2)

