n <- 1000
p <- 950
set.seed(674564587)
X <- matrix(rnorm(n * p, 0, 1/sqrt(n)), nrow = n, ncol = p, byrow = FALSE)
beta <- c(rep(3.5, 20), rep(0, p - 20))
epsilon <- rnorm(n, 0, 1)
Y <- X %*% beta + epsilon
k <- c(10, 20, 30, 50, 100, 500, 950)

fun <- function(X, beta, epsilon, k) {
  X_k <- X[, 1:k, drop = FALSE]
  beta_k <- beta[1:k, drop = FALSE]
  Y_k <- X_k %*% beta_k + epsilon
  
  return(list(X_k = X_k, beta_k = beta_k, Y_k = Y_k))
}

#zadanie 1 + 3

wyniki_i <- setNames(lapply(k, function(i) {
  dane <- fun(X, beta, epsilon, i)
  est_beta <- solve(t(dane$X_k) %*% dane$X_k) %*% t(dane$X_k) %*% dane$Y_k
  return(est_beta = est_beta)
}), paste0("Model_", k))

wyniki <- setNames(lapply(k, function(i) {
  dane <- fun(X, beta, epsilon, i)
  model <- lm(dane$Y_k ~ dane$X_k)
  model_name <- paste0("Model_", i)
  est_beta <- wyniki_i[[model_name]]
  RSS <- sum((dane$Y_k - (dane$X_k %*% est_beta))**2)
  epsilon_star <- rnorm(n, 0, 1)
  Y_star <- dane$X_k %*% dane$beta_k + epsilon_star
  PE <- sum((Y_star - (dane$X_k %*% est_beta))**2)
  sigma_est_2 <- (1/(n-i))*RSS
  PE_hat_true <- RSS + 2 * i
  PE_hat_est <- RSS + 2 * sigma_est_2 * i
  H_ii <- diag(dane$X_k %*% solve(t(dane$X_k) %*% dane$X_k) %*% t(dane$X_k))
  CV <- sum(((dane$Y_k - (dane$X_k %*% est_beta))/(1 - H_ii))**2)
  summ <- summary(model)
  return(c(RSS = RSS, PE = PE, PE_hat_true = PE_hat_true, PE_hat_est = PE_hat_est, CV = CV))
}), paste0("Model_", k))

results_c <- setNames(lapply(k, function(i) replicate(100, {
  dane <- fun(X, beta, epsilon, i)
  model_name <- paste0("Model_", i)
  est_beta <- wyniki_i[[model_name]]
  n <- nrow(dane$X_k)  
  epsilon_star <- rnorm(n, 0, 1)
  Y_star <- dane$X_k %*% dane$beta_k + epsilon_star
  PE <- sum((Y_star - (dane$X_k %*% est_beta))**2)
  return(PE)
})), paste0("Model_", k))

roznice <- setNames(lapply(k, function(i) {
  model_name <- paste0("Model_", i)
  PE <- results_c[[model_name]]
  PE_hat_true <- wyniki[[model_name]][["PE_hat_true"]]
  PE_hat_est <- wyniki[[model_name]][["PE_hat_est"]]
  CV <- wyniki[[model_name]][["CV"]]
  return(list(PE_hat_true = PE_hat_true - PE, PE_hat_est = PE_hat_est - PE, CV = CV - PE))
}), paste0("Model_", k))

library(ggplot2)

dane_boxplot <- do.call(rbind, lapply(names(roznice), function(model_name) {
  model_data <- roznice[[model_name]]
  data.frame(
    Model = model_name,
    Estymator = rep(c("PE_hat_true", "PE_hat_est", "CV"), each = length(model_data[[1]])),
    Roznica = c(model_data$PE_hat_true, model_data$PE_hat_est, model_data$CV)
  )
}))

ggplot(dane_boxplot, aes(x = Estymator, y = Roznica, fill = Estymator)) +
  geom_boxplot() +
  facet_wrap(~ Model, nrow = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Różnice estymatorów względem prawdziwego PE",
    y = "PE_hat - PE",
    x = "Estymator"
  ) +
  scale_fill_brewer(palette = "Set2")

dane_boxplot_filtered <- subset(dane_boxplot, Model != "Model_950")

ggplot(dane_boxplot_filtered, aes(x = Estymator, y = Roznica, fill = Estymator)) +
  geom_boxplot() +
  facet_wrap(~ Model, nrow = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Różnice estymatorów (bez Model_950)",
    y = "PE_hat - PE",
    x = "Estymator"
  ) +
  scale_fill_brewer(palette = "Set2")

#zadanie 2.

k1 <- c(50, 100, 200, 500, 950)

library(bigstep)

wyniki_ii <- setNames(lapply(k1, function(i) {
  dane <- fun(X, beta, epsilon, i)
  est_beta <- solve(t(dane$X_k) %*% dane$X_k) %*% t(dane$X_k) %*% dane$Y_k
  return(est_beta)
}), paste0("Model_", k1))

wynikiAIC <- setNames(lapply(k1, function(i) {
  dane <- fun(X, beta, epsilon, i)
  d <- prepare_data(dane$Y_k, dane$X_k)
  wyn1 <- fast_forward(d, crit='aic', max_variables = 1000)
  ind1 <- sort(as.numeric(wyn1$model))
  X_aic <- dane$X_k[, ind1]
  beta_hat_aic <- solve(t(X_aic)%*% X_aic)%*%t(X_aic)%*%dane$Y_k
  TD <- sum(ind1 <= 20)
  FD <- sum(ind1 > 20)
  model_name <- paste0("Model_", i)
  est_beta <- wyniki_ii[[model_name]]
  SE <- sum((dane$X_k%*%beta_hat_aic - dane$X_k %*% dane$beta_k)^2)
  list(TD = TD, FD = FD, SE = SE)
}), paste0("Model_", k1))

wynikiBIC <- setNames(lapply(k1, function(i) {
  dane <- fun(X, beta, epsilon, i)
  d <- prepare_data(dane$Y_k, dane$X_k)
  wyn1 <- fast_forward(d, crit='bic')
  ind1 <- sort(as.numeric(wyn1$model))
  X_bic <- dane$X_k[, ind1]
  beta_hat_bic <- solve(t(X_bic)%*% X_bic)%*%t(X_bic)%*%dane$Y_k
  TD <- sum(ind1 <= 20)
  FD <- sum(ind1 > 20)
  model_name <- paste0("Model_", i)
  est_beta <- wyniki_ii[[model_name]]
  SE <- sum((X_bic%*%beta_hat_bic - dane$X_k %*% dane$beta_k)^2)
  list(TD = TD, FD = FD, SE = SE)
}), paste0("Model_", k1))

wynikiRIC <- setNames(lapply(k1, function(i) {
  dane <- fun(X, beta, epsilon, i)
  d <- prepare_data(dane$Y_k, dane$X_k)
  wyn1 <- stepwise(d,crit='mbic',const=sqrt(n))
  ind1 <- sort(as.numeric(wyn1$model))
  X_ric <- dane$X_k[, ind1]
  beta_hat_ric <- solve(t(X_ric) %*% X_ric) %*% t(X_ric) %*% dane$Y_k
  TD <- sum(ind1 <= 20)
  FD <- sum(ind1 > 20)
  model_name <- paste0("Model_", i)
  est_beta <- wyniki_ii[[model_name]]
  SE <- sum((X_ric %*% beta_hat_ric - dane$X_k %*% dane$beta_k) ^2)
  list(TD = TD, FD = FD, SE = SE)
}), paste0("Model_", k1))

wynikimBIC <- setNames(lapply(k1, function(i) {
  dane <- fun(X, beta, epsilon, i)
  d <- prepare_data(dane$Y_k, dane$X_k)
  wyn1 <- stepwise(d,crit='mbic')
  ind1 <- sort(as.numeric(wyn1$model))
  X_mbic <- dane$X_k[, ind1]
  beta_hat_mbic <- solve(t(X_mbic) %*% X_mbic) %*% t(X_mbic) %*% dane$Y_k
  TD <- sum(ind1 <= 20)
  FD <- sum(ind1 > 20)
  model_name <- paste0("Model_", i)
  est_beta <- wyniki_ii[[model_name]]
  SE <- sum((X_mbic %*% beta_hat_mbic - dane$X_k %*% dane$beta_k) ^2)
  list(TD = TD, FD = FD, SE = SE)
}), paste0("Model_", k1))

wynikimBIC2 <- setNames(lapply(k1, function(i) {
  dane <- fun(X, beta, epsilon, i)
  d <- prepare_data(dane$Y_k, dane$X_k)
  wyn1 <- stepwise(d,crit='mbic2')
  ind1 <- sort(as.numeric(wyn1$model))
  X_mbic2 <- dane$X_k[, ind1]
  beta_hat_mbic2 <- solve(t(X_mbic2) %*% X_mbic2) %*% t(X_mbic2) %*% dane$Y_k
  TD <- sum(ind1 <= 20)
  FD <- sum(ind1 > 20)
  model_name <- paste0("Model_", i)
  est_beta <- wyniki_ii[[model_name]]
  SE <- sum((X_mbic2 %*% beta_hat_mbic2 - dane$X_k %*% dane$beta_k) ^2)
  list(TD = TD, FD = FD, SE = SE)
}), paste0("Model_", k1))

#3
n <- 5000
p <- 950
set.seed(47868579)
X <- matrix(rnorm(n * p, 0, 1/sqrt(n)), nrow = n, ncol = p, byrow = FALSE)
beta <- c(rep(3.5, 20), rep(0, p - 20))
epsilon1 <- rnorm(n, 0, 1)
Y <- X %*% beta + epsilon
k <- c(10, 20, 30, 50, 100, 500, 950)

wyniki_i <- setNames(lapply(k, function(i) {
  dane <- fun(X, beta, epsilon1, i)
  est_beta <- solve(t(dane$X_k) %*% dane$X_k) %*% t(dane$X_k) %*% dane$Y_k
  return(est_beta = est_beta)
}), paste0("Model_", k))

wyniki <- setNames(lapply(k, function(i) {
  dane <- fun(X, beta, epsilon1, i)
  model <- lm(dane$Y_k ~ dane$X_k)
  model_name <- paste0("Model_", i)
  est_beta <- wyniki_i[[model_name]]
  RSS <- sum((dane$Y_k - (dane$X_k %*% est_beta))**2)
  epsilon_star <- rnorm(n, 0, 1)
  Y_star <- dane$X_k %*% dane$beta_k + epsilon_star
  PE <- sum((Y_star - (dane$X_k %*% est_beta))**2)
  sigma_est_2 <- (1/(n-i))*RSS
  PE_hat_true <- RSS + 2 * i
  PE_hat_est <- RSS + 2 * sigma_est_2 * i
  H_ii <- diag(dane$X_k %*% solve(t(dane$X_k) %*% dane$X_k) %*% t(dane$X_k))
  CV <- sum(((dane$Y_k - (dane$X_k %*% est_beta))/(1 - H_ii))**2)
  summ <- summary(model)
  return(c(RSS = RSS, PE = PE, PE_hat_true = PE_hat_true, PE_hat_est = PE_hat_est, CV = CV))
}), paste0("Model_", k))

results_c <- setNames(lapply(k, function(i) replicate(100, {
  dane <- fun(X, beta, epsilon, i)
  model_name <- paste0("Model_", i)
  est_beta <- wyniki_i[[model_name]]
  n <- nrow(dane$X_k)  
  epsilon_star <- rnorm(n, 0, 1)
  Y_star <- dane$X_k %*% dane$beta_k + epsilon_star
  PE <- sum((Y_star - (dane$X_k %*% est_beta))**2)
  return(PE)
})), paste0("Model_", k))

roznice <- setNames(lapply(k, function(i) {
  model_name <- paste0("Model_", i)
  PE <- results_c[[model_name]]
  PE_hat_true <- wyniki[[model_name]][["PE_hat_true"]]
  PE_hat_est <- wyniki[[model_name]][["PE_hat_est"]]
  CV <- wyniki[[model_name]][["CV"]]
  return(list(PE_hat_true = PE_hat_true - PE, PE_hat_est = PE_hat_est - PE, CV = CV - PE))
}), paste0("Model_", k))

library(ggplot2)

dane_boxplot <- do.call(rbind, lapply(names(roznice), function(model_name) {
  model_data <- roznice[[model_name]]
  data.frame(
    Model = model_name,
    Estymator = rep(c("PE_hat_true", "PE_hat_est", "CV"), each = length(model_data[[1]])),
    Roznica = c(model_data$PE_hat_true, model_data$PE_hat_est, model_data$CV)
  )
}))

ggplot(dane_boxplot, aes(x = Estymator, y = Roznica, fill = Estymator)) +
  geom_boxplot() +
  facet_wrap(~ Model, nrow = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Różnice estymatorów względem prawdziwego PE",
    y = "PE_hat - PE",
    x = "Estymator"
  ) +
  scale_fill_brewer(palette = "Set2")