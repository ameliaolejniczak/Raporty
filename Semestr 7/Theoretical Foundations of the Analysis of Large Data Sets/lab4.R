#zadanie 1.
#a
set.seed(2137)
n_list <- c(5000, 50000)
beta_list <- c(0.6, 0.8)
r_list <- c(0.1, 0.4)
alpha <- 0.05

critical_values_L <- setNames(lapply(n_list, function(n) {
  setNames(lapply(beta_list, function(beta) {
    setNames(lapply(r_list, function(r) {
      Ls <- replicate(10000, {
        X <- rnorm(n)
        epsilon <- n ^ -beta
        mi <- sqrt(2 * r * log(n))
        L <- prod((1 - epsilon) + epsilon * exp(mi * X - (mi ^ 2)/2))
        L
      })
      critval_L <- quantile(Ls, probs = 1 - alpha)
    }), c("r = 0.1", "r = 0.4"))
  }), c("beta = 0.6", "beta = 0.8"))
}), c("n = 5000", "n = 50000"))

setwd("C:\\Users\\ameli\\Documents\\R\\theoretical foundations of large data sets")
saveRDS(critical_values_L, file = "critical_values_l.rds")

#b
hc_mod <- function(pvals) {
  n <- length(pvals)
  Fn <- ecdf(pvals)
  ts <- seq(1e-10, 1 - 1e-10, length.out = 10000)
  qt <- log(log(1/(ts * (1 - ts))))
  hc_mod_fun <- sqrt(n) * (Fn(ts) - ts)/sqrt(ts * (1 - ts) * qt)
  hc_mod <- max(hc_mod_fun)
  return(hc_mod)
}

hc <- function(pvals) {
  n <- length(pvals)
  Fn <- ecdf(pvals)
  ts <- seq(1/n + 1e-10, 1/2 - 1e-10, length.out = 10000)
  hc_fun <- sqrt(n) * (Fn(ts) - ts)/sqrt(ts * (1 - ts))
  hc <- max(hc_fun)
  return(hc)
}

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

neyman_pearson <- function(sample, beta, r) {
  n <- length(sample)
  epsilon <- n ^ -beta
  mi <- sqrt(2 * r * log(n))
  L <- prod((1 - epsilon) + epsilon * exp(mi * sample - (mi ^ 2)/2))
  return(L)
}

critical_values <- setNames(lapply(n_list, function(n)
{
  critvals <- replicate(10000, {
    p_values <- runif(n)
    hc <- hc(p_values)
    hc_mod <- hc_mod(p_values)
    c(hc, hc_mod)
  })
  critval_hc <- quantile(critvals[1, ], probs = 1 - alpha)
  critval_hc_mod <- quantile(critvals[2, ], probs = 1 - alpha)
  list(critval_hc = critval_hc, critval_hc_mod = critval_hc_mod)
}),  c("n = 5000", "n = 50000"))

powers <- setNames(lapply(n_list, function(n) {
  setNames(lapply(beta_list, function(beta) {
    setNames(lapply(r_list, function(r) {
      epsilon <- n ^ -beta
      mi <- sqrt(2 * r * log(n))
      tests <- replicate(10000, {
        signal <- rbinom(n, 1, epsilon)
        sample <- rnorm(n, mean = mi * signal)
        p_vals <- 2 * (1 - pnorm(abs(sample)))
        hc_mod_test <- hc_mod(p_vals) > critical_values[[paste0("n = ", n)]][["critval_hc_mod"]]
        hc_test <- hc(p_vals) > critical_values[[paste0("n = ", n)]][["critval_hc"]]
        bonf_test <- bonferroni(p_vals) <= alpha/n
        chi_sq_test <- chi_square(sample) > qchisq(1 - alpha, n)
        fish_test <- fisher(p_vals) > qchisq(1 - alpha, 2 * n)
        neym_pear_test <- neyman_pearson(sample, beta, r) > critical_values_L[[paste0("n = ", n)]][[paste0("beta = ", beta)]][[paste0("r = ", r)]]
        c(hc_mod_test = hc_mod_test, hc_test = hc_test, bonf_test = bonf_test, chi_sq_test = chi_sq_test, fish_test = fish_test, neym_pear_test = neym_pear_test)
      })
      apply(tests, 1, FUN = mean)
    }), c("r = 0.1", "r = 0.4"))
  }), c("beta = 0.6", "beta = 0.8"))
}), c("n = 5000", "n = 50000"))

saveRDS(powers, "powers.rds")

powers_to_table <- function(powers) {
  rows <- list()
  
  for (n_name in names(powers)) {
    for (beta_name in names(powers[[n_name]])) {
      for (r_name in names(powers[[n_name]][[beta_name]])) {
        
        vec <- powers[[n_name]][[beta_name]][[r_name]]
        
        rows[[length(rows) + 1]] <- data.frame(
          n = gsub("n = ", "", n_name),
          beta = gsub("beta = ", "", beta_name),
          r = gsub("r = ", "", r_name),
          test = names(vec),
          power = as.numeric(vec),
          row.names = NULL
        )
      }
    }
  }
  
  do.call(rbind, rows)
}

power_table <- powers_to_table(powers)
power_table

library(dplyr)
library(ggplot2)

plot_data <- power_table %>%
  mutate(
    x = paste0("b=", beta, ", r=", r),
    x = factor(x, levels = unique(x)),
    n = factor(n),
    test = factor(test)
  )

ggplot(plot_data,
       aes(x = x, y = power, color = test, group = test)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ n, nrow = 1,
             labeller = labeller(n = function(x) paste("n =", x))) +
  labs(
    x = "",
    y = "Moc testu",
    color = "",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


library(tidyr)

power_table_wide <- pivot_wider(
  power_table,
  names_from = test,
  values_from = power
)

power_table_wide

#zadanie 2.
n <- 20
N <- 10000
alpha <- 0.05

mu1 <- c(1.5 * sqrt(2 * log(n)), rep(0, n - 1))
mu2 <- c(rep(1.05 * sqrt(2 * log(n)), 5), rep(0, n - 5))
mu3 <- c(1.5 * sqrt(2 * log(20/1:10)), rep(0, n - 10))
mu_list <- list(mu1 = mu1, mu2 = mu2, mu3 = mu3)

powers_fdr_fwer_results <- lapply(mu_list, function(mu) {
  test_results <- replicate(N, {
    X <- rnorm(n, mu, 1)
    pvals <- 2 * (1 - pnorm(abs(X)))
    alpha_sidak <- 1 - (1 - alpha)^(1/n)
    
    tests <- list(Bonferroni = p.adjust(pvals, method = "bonferroni") <= alpha,
    Sidak = pvals <= alpha_sidak,
    Holm = p.adjust(pvals, method = "holm") <= alpha,
    Hochberg = p.adjust(pvals, method = "hochberg") <= alpha,
    BH = p.adjust(pvals, method = "BH") <= alpha)
    
    lapply(tests, function(x) {
      true_alt <- mu != 0
      rejected <- x
      FP <- sum(rejected & !true_alt)
      TP <- sum(rejected & true_alt)
      R <- sum(rejected)
      mu1 <- sum(true_alt)
      
      FWER <- as.numeric(FP>0)
      FDR <- ifelse(R==0, 0, FP/R)
      Power <- ifelse(mu1==0, 0, TP/mu1)
      
      c(FWER=FWER, FDR=FDR, Power=Power)
    })
  })
})

library(dplyr)

summarize_results <- function(res) {
  methods <- rownames(res)
  metrics <- c("FWER", "FDR", "Power")
  
  out <- matrix(NA, nrow = length(methods), ncol = length(metrics),
                dimnames = list(methods, metrics))
  
  for (i in seq_along(methods)) {
    mat <- do.call(cbind, res[i, ])
    out[i, ] <- rowMeans(mat)
  }
  
  out
}

final_results <- lapply(powers_fdr_fwer_results, summarize_results)
final_results

saveRDS(final_results, "results_low.rds")

#zadanie 3.
n <- 5000
N <- 10000
alpha <- 0.05

mu1 <- c(1.5 * sqrt(2 * log(n)), rep(0, n - 1))
mu2 <- c(rep(1.05 * sqrt(2 * log(n)), 100), rep(0, n - 100))
mu3 <- c(rep(1.005 * sqrt(2 * log(n)), 100), rep(0, n - 100))
mu4 <- c(rep(1.005 * sqrt(2 * log(n/100)), 1000), rep(0, n - 1000))
mu_list <- list(mu1 = mu1, mu2 = mu2, mu3 = mu3, mu4 = mu4)

powers_fdr_fwer_results <- lapply(mu_list, function(mu) {
  test_results <- replicate(N, {
    X <- rnorm(n, mu, 1)
    pvals <- 2 * (1 - pnorm(abs(X)))
    alpha_sidak <- 1 - (1 - alpha)^(1/n)
    
    tests <- list(Bonferroni = p.adjust(pvals, method = "bonferroni") <= alpha,
                  Sidak = pvals <= alpha_sidak,
                  Holm = p.adjust(pvals, method = "holm") <= alpha,
                  Hochberg = p.adjust(pvals, method = "hochberg") <= alpha,
                  BH = p.adjust(pvals, method = "BH") <= alpha)
    
    lapply(tests, function(x) {
      true_alt <- mu != 0
      rejected <- x
      FP <- sum(rejected & !true_alt)
      TP <- sum(rejected & true_alt)
      R <- sum(rejected)
      mu1 <- sum(true_alt)
      
      FWER <- as.numeric(FP>0)
      FDR <- ifelse(R==0, 0, FP/R)
      Power <- ifelse(mu1==0, 0, TP/mu1)
      
      c(FWER=FWER, FDR=FDR, Power=Power)
    })
  })
})

library(dplyr)

summarize_results <- function(res) {
  methods <- rownames(res)
  metrics <- c("FWER", "FDR", "Power")
  
  out <- matrix(NA, nrow = length(methods), ncol = length(metrics),
                dimnames = list(methods, metrics))
  
  for (i in seq_along(methods)) {
    mat <- do.call(cbind, res[i, ])
    out[i, ] <- rowMeans(mat)
  }
  
  out
}

final_results <- lapply(powers_fdr_fwer_results, summarize_results)
final_results

saveRDS(final_results, "results_large.rds")
