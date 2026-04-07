# zadanie 1.
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(411)

#a
calculating_pvalue <- function(sample) {
  n <- length(sample)
  sum_s <- sum(sample)
  1 - ppois(sum_s - 1, n * 5)
}

#b
p_values <- replicate(1000, {
  sample <- rpois(100, 5)
  calculating_pvalue(sample)
})

df <- data.frame(p_values = p_values)

ggplot(df, aes(x = p_values)) +
  geom_histogram(bins = 20, fill = "skyblue", alpha = 0.8, aes(y = ..density..)) +
  geom_hline(yintercept = 1, col = "black") + 
  theme_bw() +
  labs(title = "Histogram p-values",
       x = "p-value",
       y = "Density")

#c
N <- 1000
alpha <- 0.05
n <- 100

testing_c <- replicate(N, {
  pvals <- replicate(1000, {
    sample <- rpois(100, 5)
    calculating_pvalue(sample)
  })
  bonf <- min(pvals) <= 0.05/1000
  T_stat <- - sum(2 * log(pvals))
  fish <- T_stat > qchisq(1 - alpha, 2000)
  c(bonf = bonf, fish = fish)
})

setwd("C:\\Users\\ameli\\Documents\\R\\theoretical foundations of large data sets")
saveRDS(testing_c, file = "testing_c.rds")

mean(testing_c[1, ])
mean(testing_c[2, ])

#d
needle_in_the_hay <- replicate(N, {
  pvals <- c(
  calculating_pvalue(rpois(100, 7)), 
  replicate(999, {
    sample <- rpois(100, 5)
    calculating_pvalue(sample)
  }))
  bonf <- min(pvals) <= 0.05/1000
  T_stat <- - sum(2 * log(pvals))
  fish <- T_stat > qchisq(1 - alpha, 2000)
  c(bonf = bonf, fish = fish)
})

saveRDS(needle_in_the_hay, "needle_in_the_hay.rds")

mean(needle_in_the_hay[1, ])
mean(needle_in_the_hay[2, ])

many_small_eff <- replicate(N, {
  pvals <- c(replicate(100, {
    sample <- rpois(100, 5.2)
    calculating_pvalue(sample)
  }), replicate(900, {
    sample <- rpois(100, 5)
    calculating_pvalue(sample)
  }))
  bonf <- min(pvals) <= 0.05/1000
  T_stat <- - sum(2 * log(pvals))
  fish <- T_stat > qchisq(1 - alpha, 2000)
  c(bonf = bonf, fish = fish)
})

saveRDS(many_small_eff, "many_small_eff.rds")

mean(many_small_eff[1, ])
mean(many_small_eff[2, ])

#zadanie 2.
rn_fun <- function() {
  x <- rnorm(100000)
  N <- 2:100000
  sapply(N, function(n) {max(x[1:n])/sqrt(2 * log(n))})
}

Rs <- lapply(1:10, function(x) {
  rn_fun()
})

saveRDS(Rs, "Rs.rds")

df <- data.frame(N = 2:100000, do.call(cbind, Rs))
colnames(df)[-1] <- paste0("Run_", 1:10)

df_long <- df %>%
  pivot_longer(cols = starts_with("Run_"), names_to = "Run", values_to = "Value")

ggplot(df_long, aes(x = N, y = Value, color = Run)) +
  geom_line(size = 0.8) +
  labs(title = "rn trajectories",
       x = "N",
       y = "max(x[1:n]) / sqrt(2 log(n))") +
  theme_bw() +
  theme(legend.position = "right")

#zadanie 3.
set.seed(12)
n_list <- c(1000, 10000, 100000)

Ls <- lapply(n_list, function(n) {
  replicate(1000, {
    Y <- rnorm(n)
    epsilon <- 0.1
    gamma <- (1 - epsilon) * sqrt(2 * log(n))
    L <- mean(exp(gamma * Y - (gamma^2)/2))
    Y_hat <- Y[Y < sqrt(2 * log(n))]
    L_hat <- mean(exp(gamma * Y_hat - (gamma^2)/2))
    c(L = L, L_tilde = L_hat)  
    })
})

saveRDS(Ls, "Ls.rds")

#a

df_list <- lapply(seq_along(n_list), function(i) {
  mat <- t(Ls[[i]])
  df <- as.data.frame(mat)
  df$n <- n_list[i]
  df
})

df_all <- bind_rows(df_list) %>%
  pivot_longer(cols = c("L", "L_tilde"), names_to = "Type", values_to = "Value")

ggplot(df_all, aes(x = Value, fill = Type)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  facet_wrap(~n) +
  labs(title = "L and L_hat histograms",
       fill = "Type") +
  theme_bw()

#b

variance <- lapply(Ls, function(l) {
  var_L <- var(l[1, ])
  var_L_hat <- var(l[2, ])
  list(var_L = var_L, var_L_hat = var_L_hat)
})

#c 

theoretical_probabilities <- lapply(n_list, function(n) {
  value <- sqrt(2 * log(n))
  pnorm(value)^n
})

probabilities <- lapply(Ls, function(l) {
  mean(l[1, ] == l[2, ])
})

#zadanie 4.

n_values <- c(500, 5000, 50000)
B <- 10000
epsilon_values <- c(0.05, 0.2)
alpha <- 0.05

results <- lapply(n_values, function(n) {
  lapply(epsilon_values, function(eps) {
    gamma <- (1 + eps) * sqrt(2 * log(n))
    log_Ls <- replicate(B, {
      Y_h0 <- rnorm(n)
      L <- log(mean(exp(gamma * Y_h0 - (gamma^2)/2)))
    })
    critical_value <- quantile(log_Ls, probs = 1 - alpha)
    tests <- replicate(B, {
      Y_h1 <- rnorm(n)
      Y_h1[sample(1:n, 1)] <- rnorm(1, mean = gamma)
      L <- log(mean(exp(gamma * Y_h1 - (gamma^2)/2)))
      L_check <- L > critical_value
      Bonf_check <- max(abs(Y_h1)) > abs(qnorm(alpha/(2 * n)))
      c(L_check = L_check, Bonf_check = Bonf_check)
    }, simplify = TRUE)
    power_L <- mean(tests[1, ])
    power_bonf <- mean(tests[2, ])
    list(critical_value = critical_value, power_L = power_L, power_bonf = power_bonf)
  })
})

saveRDS(results, "results.rds")

df <- bind_rows(
  lapply(seq_along(n_values), function(i) {
    bind_rows(
      lapply(seq_along(epsilon_values), function(j) {
        x <- results[[i]][[j]]
        data.frame(
          n = n_values[i],
          epsilon = epsilon_values[j],
          critical_value = x$critical_value,
          power_L = x$power_L,
          power_bonf = x$power_bonf
        )
      })
    )
  })
)

df_long <- df %>%
  pivot_longer(cols = c(power_L, power_bonf),
               names_to = "method",
               values_to = "power")

ggplot(df_long, aes(x = n, y = power, color = method)) +
  geom_point(size = 4) +
  facet_wrap(~ epsilon) +
  theme_bw() +
  labs(title = "Porównanie mocy testów",
       x = "n",
       y = "Power",
       color = "Metoda")

#zadanie 5.

#a 
x_range <- c(-3, 3)

ggplot(data.frame(x = x_range), aes(x)) +
  stat_function(fun = function(x) pt(x, 1), aes(color = "t(1)"), size = 1) +
  stat_function(fun = function(x) pt(x, 3), aes(color = "t(3)"), size = 1) +
  stat_function(fun = function(x) pt(x, 5), aes(color = "t(5)"), size = 1) +
  stat_function(fun = function(x) pt(x, 10), aes(color = "t(10)"), size = 1) +
  stat_function(fun = function(x) pt(x, 50), aes(color = "t(50)"), size = 1) +
  stat_function(fun = function(x) pt(x, 100), aes(color = "t(100)"), size = 1) +
  stat_function(fun = pnorm, aes(color = "N(0,1)"), size = 1) +
  scale_color_manual(
    name = "dist",
    values = c(
      "N(0,1)" = "black",
      "t(1)" = "#FFDE8B",
      "t(3)" = "#FFA88B",
      "t(5)" = "#FF6A8B",
      "t(10)" = "#C874AA",
      "t(50)" = "#AD8CFF",
      "t(100)" = "#90CFFF"
    )
  ) +
  coord_cartesian(xlim = x_range) +
  labs(title = "normal and student",
       x = "x",
       y = "cdf") +
  theme_bw() +
  theme(legend.position = "right")


#b
stand_chisq_cdf <- function(x, df) {
  pchisq(df + x * sqrt(2 * df), df)
}

ggplot(data.frame(x = x_range), aes(x)) +
  stat_function(fun = function(x) stand_chisq_cdf(x, 1), aes(color = "standchisq(1)"), size = 1) +
  stat_function(fun = function(x) stand_chisq_cdf(x, 3), aes(color = "standchisq(3)"), size = 1) +
  stat_function(fun = function(x) stand_chisq_cdf(x, 5), aes(color = "standchisq(5)"), size = 1) +
  stat_function(fun = function(x) stand_chisq_cdf(x, 10), aes(color = "standchisq(10)"), size = 1) +
  stat_function(fun = function(x) stand_chisq_cdf(x, 50), aes(color = "standchisq(50)"), size = 1) +
  stat_function(fun = function(x) stand_chisq_cdf(x, 100), aes(color = "standchisq(100)"), size = 1) +
  stat_function(fun = pnorm, aes(color = "N(0,1)"), size = 1) +
  scale_color_manual(
    name = "dist",
    values = c(
      "N(0,1)" = "black",
      "standchisq(1)" = "#FFDE8B",
      "standchisq(3)" = "#FFA88B",
      "standchisq(5)" = "#FF6A8B",
      "standchisq(10)" = "#C874AA",
      "standchisq(50)" = "#AD8CFF",
      "standchisq(100)" = "#90CFFF"
    )
  ) +
  coord_cartesian(xlim = x_range) +
  labs(title = "normal and standardised chi square",
       x = "x",
       y = "cdf") +
  theme_bw() +
  theme(legend.position = "right")

