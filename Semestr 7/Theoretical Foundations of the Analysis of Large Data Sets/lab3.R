#zadanie 1.
library(goftest)
library(dplyr)
library(ggplot2)
setwd("C:\\Users\\ameli\\Documents\\R\\theoretical foundations of large data sets")

set.seed(12345)

hc_mod <- function(pvals) {
  n <- length(pvals)
  Fn <- ecdf(pvals)
  ts <- seq(1e-10, 1 - 1e-10, length.out = 10000)
  qt <- log(log(1/(ts * (1 - ts))))
  hc_mod_fun <- sqrt(n) * (Fn(ts) - ts)/sqrt(ts * (1 - ts) * qt)
  hc_mod <- max(hc_mod_fun)
  return(hc_mod)
}

n <- c(5000, 50000)

type_I_errors <- lapply(n, function(n) {
  rejected <- replicate(10000, {  
    pvals <- runif(n)
    hc_mod_vals <- hc_mod(pvals)
    hc_mod_vals > 4.14
  })
  mean(rejected)
})

saveRDS(type_I_errors, file = "type_I_errors.rds")
#dodatkowe wykresy

n <- 5000

p_vals <- runif(n)
Fn <- ecdf(p_vals)
ts <- p_vals
qt <- log(log(1/(ts * (1 - ts))))
fun_in_hc_mod <- sqrt(n) * (Fn(ts) - ts)/sqrt(ts * (1 - ts) * qt)
fun_in_hc <- sqrt(n) * ((Fn(ts) - ts)/(sqrt(ts * (1 - ts))))

df <- data.frame(
  p_vals = p_vals,
  fun_in_hc_mod = fun_in_hc_mod,
  fun_in_hc = fun_in_hc
)

ggplot(df, aes(x = p_vals, y = fun_in_hc)) +
  geom_line() +
  xlim(1/n, 1/2)

ggplot(df, aes(x = p_vals, y = fun_in_hc_mod)) +
  geom_line() +
  xlim(0, 1)

#zadanie 2.
n <- 5000
alpha <- 0.05

hc <- function(pvals) {
  n <- length(pvals)
  Fn <- ecdf(pvals)
  ts <- seq(1/n + 1e-10, 1/2 - 1e-10, length.out = 10000)
  hc_fun <- sqrt(n) * (Fn(ts) - ts)/sqrt(ts * (1 - ts))
  hc <- max(hc_fun)
  return(hc)
}

critical_values <- replicate(10000, {
  p_values <- runif(n)
  hc <- hc(p_values)
  hc_mod <- hc_mod(p_values)
  c(hc, hc_mod)
})

saveRDS(critical_values, "critical_values.rds")

critval_hc <- quantile(critical_values[1, ], probs = 1 - alpha)
critval_hc_mod <- quantile(critical_values[2, ], probs = 1 - alpha)

# zadanie 3. 
bonferroni <- function(pvals) {
  n <- length(pvals)
  return(min(pvals))
}

chi_square <- function(samp) {
  n <- length(samp)
  T_stat <- sum(samp^2)
  return(T_stat)
}

fisher <- function(pvals) {
  n <- length(pvals)
  T_stat <- - sum(2 * log(pvals))
  return(T_stat)
}

n <- 5000
mu <- c(1.2 * sqrt(2 * log(n)), 1.02 * sqrt(2 * log(n/2000)), 1.002 * sqrt(2 * log(n/4000)))

test_values <- replicate(10000, {
  sample_1 <- c(rnorm(1, mu[1], 1), rnorm(n - 1, 0, 1))
  sample_2 <- c(rnorm(100, mu[2], 1), rnorm(n - 100, 0, 1))
  sample_3 <- c(rnorm(1000, mu[3], 1), rnorm(n - 1000, 0, 1))
  
  pvals_1 <- 2 * (1 - pnorm(abs(sample_1)))
  pvals_2 <- 2 * (1 - pnorm(abs(sample_2)))
  pvals_3 <- 2 * (1 - pnorm(abs(sample_3)))
  
  hc_1 <- hc(pvals_1) > critval_hc
  hc_2 <- hc(pvals_2) > critval_hc
  hc_3 <- hc(pvals_3) > critval_hc
  
  hc_mod_1 <- hc_mod(pvals_1) > critval_hc_mod
  hc_mod_2 <- hc_mod(pvals_2) > critval_hc_mod
  hc_mod_3 <- hc_mod(pvals_3) > critval_hc_mod
  
  bonf_1 <- bonferroni(pvals_1) <= alpha/n
  bonf_2 <- bonferroni(pvals_2) <= alpha/n
  bonf_3 <- bonferroni(pvals_3) <= alpha/n
  
  chi_sq_1 <- chi_square(sample_1) > qchisq(1 - alpha, n)
  chi_sq_2 <- chi_square(sample_2) > qchisq(1 - alpha, n)
  chi_sq_3 <- chi_square(sample_3) > qchisq(1 - alpha, n)
  
  fish_1 <- fisher(pvals_1) > qchisq(1 - alpha, 2 * n)
  fish_2 <- fisher(pvals_2) > qchisq(1 - alpha, 2 * n)
  fish_3 <- fisher(pvals_3) > qchisq(1 - alpha, 2 * n)
  
  ks_1 <- ks.test(pvals_1, y = "punif")$p.value <= alpha
  ks_2 <- ks.test(pvals_2, y = "punif")$p.value <= alpha
  ks_3 <- ks.test(pvals_3, y = "punif")$p.value <= alpha
  
  ad_1 <- ad.test(pvals_1)$p.value <= alpha
  ad_2 <- ad.test(pvals_2)$p.value <= alpha
  ad_3 <- ad.test(pvals_3)$p.value <= alpha
  
  c(hc_1 = hc_1, hc_mod_1 = hc_mod_1, bonf_1 = bonf_1, chi_sq_1 = chi_sq_1, fish_1 = fish_1, ks_1 = ks_1, ad_1 = ad_1,
    hc_2 = hc_2, hc_mod_2 = hc_mod_2, bonf_2 = bonf_2, chi_sq_2 = chi_sq_2, fish_2 = fish_2, ks_2 = ks_2, ad_2 = ad_2,
    hc_3 = hc_3, hc_mod_3 = hc_mod_3, bonf_3 = bonf_3, chi_sq_3 = chi_sq_3, fish_3 = fish_3, ks_3 = ks_3, ad_3 = ad_3)
})

test_powers <- apply(test_values, 1, FUN = mean)
test_powers_df <- data.frame(
  hc = c(test_powers["hc_1.95%"], test_powers["hc_2.95%"], test_powers["hc_3.95%"]),
  hc_mod = c(test_powers["hc_mod_1.95%"], test_powers["hc_mod_2.95%"], test_powers["hc_mod_3.95%"]),
  bonf = c(test_powers["bonf_1"], test_powers["bonf_2"], test_powers["bonf_3"]),
  chi2 = c(test_powers["chi_sq_1"], test_powers["chi_sq_2"], test_powers["chi_sq_3"]),
  fisher = c(test_powers["fish_1"], test_powers["fish_2"], test_powers["fish_3"]),
  ks = c(test_powers["ks_1"], test_powers["ks_2"], test_powers["ks_3"]),
  ad = c(test_powers["ad_1"], test_powers["ad_2"], test_powers["ad_3"])
)

rownames(test_powers_df) <- c("mu_1", "mu_2", "mu_3")

saveRDS(test_powers_df, "test_powers_df.rds")

#zadanie 4.
library(sde)
n = 5000
N = 1000
t_seq <- seq(0, 1, length.out = 1000)

trajectories <- lapply(1:N, function(x) {
  X <- runif(n)
  Fn <- ecdf(X)
  Un <- sqrt(n) * (Fn(t_seq) - t_seq)
  B <- BBridge(N = length(t_seq) - 1)
  list(Un = Un, B = B)
})
saveRDS(trajectories, "trajectories.rds")

idx <- sample(1:N, 5)
subset_traj <- trajectories[idx]

df_Un <- lapply(1:5, function(i) {
  data.frame(
    t = t_seq,
    value = subset_traj[[i]]$Un,
    traj = paste0("traj_", i),
    type = "Un"
  )
}) %>% bind_rows()

df_B <- lapply(1:5, function(i) {
  data.frame(
    t = t_seq,
    value = subset_traj[[i]]$B,
    traj = paste0("traj_", i),
    type = "B"
  )
}) %>% bind_rows()

df_all <- bind_rows(df_Un, df_B)

ggplot(df_all, aes(x = t, y = value, color = traj)) +
  geom_line() +
  facet_wrap(~ type, scales = "free_y") +
  theme_minimal() +
  labs(title = "Trajektorie Un i mostu Browna B",
       x = "t", y = "") +
  theme(legend.position = "none")

K_S_vals <- sapply(trajectories, function(traj) max(abs(traj$Un)))
T_vals <- sapply(trajectories, function(traj) max(abs(traj$B)))

alpha_levels <- c(0.8, 0.9, 0.95)
KS_quantile <- quantile(K_S_vals, probs = alpha_levels)
T_quantile <- quantile(T_vals, probs = alpha_levels)
