# Zadanie 2.
m_list <- c(20, 100)
n_list <- c(200, 1000)
eps_list <- c(0.05, 0.1, 0.2)
c_0 <- 1
c_1 <- 1

results <- lapply(n_list, function(n) {
  lapply(m_list, function(m) {
    alpha_list <- c(0.1, 0.1 * sqrt(200/m))
    lapply(eps_list, function(eps) {
      lapply(alpha_list, function(alpha) {
        results <- replicate(1000, {
          epsilon <- rbinom(n, 1, eps)
          X <- rgamma(n, m, 1/3) * (1 - epsilon) + rgamma(n, m, 1/5.5)*epsilon
          p_vals <- 1 - pgamma(X, m, 1/3)
          
          #Bonferroni
          results_Bon <- p_vals < alpha/n
          FD_Bon <- results_Bon & (1 - epsilon)
          TD_Bon <- results_Bon & epsilon
          bl_II_Bon <- (1 - results_Bon) & epsilon
          FWER_Bon <- sum(FD_Bon) >= 1
          FDR_Bon <- sum(FD_Bon)/max(sum(results_Bon), 1)
          Power_Bon <- sum(TD_Bon)/sum(epsilon)
          E_Bon <- c_0 * sum(FD_Bon) + c_1 * sum(bl_II_Bon)
          
          #Benjamini-Hochberg
          BH_threshold <- alpha * (n:1)/n
          p_sorted <- sort(p_vals, decreasing = TRUE)
          results_BH <- p_vals <= p_sorted[which.max(p_sorted <= BH_threshold)]
          FD_BH <- results_BH & (1 - epsilon)
          TD_BH <- results_BH & epsilon
          bl_II_BH <- (1 - results_BH) & epsilon
          FWER_BH <- sum(FD_BH) >= 1
          FDR_BH <- sum(FD_BH)/max(sum(results_BH), 1)
          Power_BH <- sum(TD_BH)/sum(epsilon)
          E_BH <- c_0 * sum(FD_BH) + c_1 * sum(bl_II_BH)
          
          #BFDR
          BFDR <- function(x) {
            F0 <- pgamma(x, m, 1/3)
            F1 <- (1 - eps) * pgamma(x, m, 1/3) + eps * pgamma(x, m, 1/5.5) 
            return((1 - eps) * (1 - F0)/(1 - F1) - alpha)
            }
          if (m == 20) {g <- 150} else {g <- 500}
          cBFDR <- uniroot(BFDR, interval = c(0.0001,g))$root
          result_BFDR <- X >= cBFDR
          FD_BFDR <- result_BFDR & (1 - epsilon)
          TD_BFDR <- result_BFDR & epsilon
          bl_II_BFDR <- (1 - result_BFDR) & epsilon
          FWER_BFDR <- sum(FD_BFDR) >= 1
          FDR_BFDR <- sum(FD_BFDR)/max(sum(result_BFDR), 1)
          Power_BFDR <- sum(TD_BFDR)/sum(epsilon)
          E_BFDR <- c_0*sum(FD_BFDR)+c_1*sum(bl_II_BFDR)
          
          #Bayes_Class
          thres_BC <- log((1 - eps)/eps * (5.5/3)^m)/(1/3 - 1/5.5)
          result_BC <- X >= thres_BC
          FD_BC <- result_BC & (1 - epsilon)
          TD_BC <- result_BC & epsilon
          bl_II_BC <- (1 - result_BC) & epsilon
          FWER_BC <- sum(FD_BC) >= 1
          FDR_BC <- sum(FD_BC)/max(sum(result_BC), 1)
          Power_BC <- sum(TD_BC)/sum(epsilon)
          E_BC <- c_0 * sum(FD_BC) + c_1 * sum(bl_II_BC)
          
          c(FWER_Bon = FWER_Bon, 
            FDR_Bon = FDR_Bon, 
            Power_Bon = Power_Bon, 
            E_Bon = E_Bon, 
            FWER_BH = FWER_BH, 
            FDR_BH = FDR_BH, 
            Power_BH = Power_BH, 
            E_BH = E_BH, 
            FWER_BFDR = FWER_BFDR, 
            FDR_BFDR = FDR_BFDR, 
            Power_BFDR = Power_BFDR, 
            E_BFDR = E_BFDR, 
            FWER_BC = FWER_BC, 
            FDR_BC = FDR_BC, 
            Power_BC = Power_BC, 
            E_BC = E_BC)
        })
        apply(results, 1, mean)
      })
    })
  })
})

df_results <- data.frame()

for (i_n in seq_along(results)) {
  for (i_m in seq_along(results[[i_n]])) {
    for (i_eps in seq_along(results[[i_n]][[i_m]])) {
      for (i_alpha in seq_along(results[[i_n]][[i_m]][[i_eps]])) {
        
        res <- results[[i_n]][[i_m]][[i_eps]][[i_alpha]]
        
        df_results <- rbind(
          df_results,
          cbind(
            n     = n_list[i_n],
            m     = m_list[i_m],
            eps   = eps_list[i_eps],
            alpha = c(0.1, 0.1 * sqrt(200/m_list[i_m]))[i_alpha],
            as.data.frame(t(res))
          )
        )
      }
    }
  }
}

library(dplyr)
library(tidyr)

df_long <- df_results %>%
  pivot_longer(
    cols = -c(n, m, eps, alpha),
    names_to = c("Metric", "Test"),
    names_sep = "_",
    values_to = "Value"
)

setwd("C:\\Users\\ameli\\Documents\\R\\theoretical foundations of large data sets")
saveRDS(df_long, "df_long.rds")

library(ggplot2)

plot_metric <- function(metric_name) {
  ggplot(
    df_long %>% filter(Metric == metric_name),
    aes(
      x = eps,
      y = Value,
      color = Test,
      shape = factor(alpha),
      group = interaction(Test, alpha)
    )
  ) +
    geom_line() +
    geom_point(size = 2) +
    facet_grid(m ~ n, labeller = label_both) +
    labs(
      x = expression(epsilon),
      y = metric_name,
      shape = expression(alpha),
      color = "Test"
    ) +
    theme_bw() +
    theme(
      legend.position = "right"
    )
}

plot_metric("FWER")
plot_metric("FDR")
plot_metric("Power")
plot_metric("E")

theor_values <- lapply(n_list, function(n) {
  lapply(m_list, function(m) {
    alpha_list <- c(0.1, 0.1 * sqrt(200/m))
    lapply(eps_list, function(eps) {
      lapply(alpha_list, function(alpha) {
        #Bonferroni
        thres_bonf <- qgamma(1 - alpha/n, m, 1/3)
        Power_theor_bonf <- 1 - pgamma(thres_bonf, m, 1/5.5)
        cost_theor_bonf <- (c_0 * (1 - eps) * (1 - pgamma(thres_bonf, m, 1/3)) + c_1 * eps * pgamma(thres_bonf, m, 1/5.5)) * n
        
        #BFDR
        BFDR <- function(x) {
          F0 <- pgamma(x, m, 1/3)
          F1 <- (1 - eps) * pgamma(x, m, 1/3) + eps * pgamma(x, m, 1/5.5) 
          return((1 - eps) * (1 - F0)/(1 - F1) - alpha)
        }
        if (m == 20) {g <- 150} else {g <- 500}
        cBFDR <- uniroot(BFDR, interval = c(0.0001,g))$root
        Power_theor_BFDR <- 1 - pgamma(cBFDR, m, 1/5.5)
        cost_theor_BFDR <- (c_0 * (1 - eps) * (1 - pgamma(cBFDR, m, 1/3)) + c_1 * eps * pgamma(cBFDR, m, 1/5.5)) * n
        
        #BC
        thres_BC <- log((1 - eps)/eps * (5.5/3)^m)/(1/3 - 1/5.5)
        Power_theor_BC <- 1 - pgamma(thres_BC, m, 1/5.5)
        cost_theor_BC <- (c_0 * (1 - eps) * (1 - pgamma(thres_BC, m, 1/3)) + c_1 * eps * pgamma(thres_BC, m, 1/5.5)) * n
        
        list(Power_theor_bonf = Power_theor_bonf, 
             cost_theor_bonf = cost_theor_bonf, 
             Power_theor_BFDR = Power_theor_BFDR, 
             cost_theor_BFDR = cost_theor_BFDR, 
             Power_theor_BC = Power_theor_BC, 
             cost_theor_BC = cost_theor_BC)
      })
    })
  })
})       

df_theor_values <- do.call(rbind, lapply(seq_along(n_list), function(i_n) {
  do.call(rbind, lapply(seq_along(m_list), function(i_m) {
    do.call(rbind, lapply(seq_along(eps_list), function(i_eps) {
      do.call(rbind, lapply(1:2, function(i_alpha) {
        
        res <- theor_values[[i_n]][[i_m]][[i_eps]][[i_alpha]]
        
        alpha <- c(0.1, 0.1 * sqrt(200/m_list[i_m]))[i_alpha]
        
        data.frame(
          n = n_list[i_n],
          m = m_list[i_m],
          eps = eps_list[i_eps],
          alpha = alpha,
          Power_theor_bonf = res$Power_theor_bonf,
          cost_theor_bonf = res$cost_theor_bonf,
          Power_theor_BFDR = res$Power_theor_BFDR,
          cost_theor_BFDR = res$cost_theor_BFDR,
          Power_theor_BC = res$Power_theor_BC,
          cost_theor_BC = res$cost_theor_BC
        )
        
      }))
    }))
  }))
}))

library(knitr)

df_theor_long <- df_theor_values %>%
  pivot_longer(
    cols = -c(n, m, eps, alpha),
    names_to = c("Metric", "Test"),
    names_sep = "_theor_",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(Metric,
                    Power = "Power",
                    cost  = "E"),      # koszt traktujemy jak E
    Test = recode(Test,
                  bonf = "Bon",
                  BFDR = "BFDR",
                  BC   = "BC")
  )

saveRDS(df_theor_long, "df_theor_long.rds")

plot_metric_theor <- function(metric_name) {
  ggplot(
    df_theor_long %>% filter(Metric == metric_name),
    aes(
      x = eps,
      y = Value,
      color = Test,
      shape = factor(alpha),
      group = interaction(Test, alpha)
    )
  ) +
    geom_line() +
    geom_point(size = 2) +
    facet_grid(m ~ n, labeller = label_both) +
    labs(
      x = expression(epsilon),
      y = paste("Theoretical", metric_name),
      shape = expression(alpha),
      color = "Test"
    ) +
    theme_bw()
}

plot_metric_theor("Power")
plot_metric_theor("E")

# Zadanie 3.
EM_algorithm <- function(n, m, eps, k_pow) {
  epsi <- numeric(k_pow)
  mu_1 <- numeric(k_pow)
  mu_0 <- numeric(k_pow)
  for (i in 1:k_pow){
    altern <- rbinom(n, 1, eps)
    X <- rgamma(n, m, 1/3) * (1 - altern) + rgamma(n, m, 1/5.5) * altern
    eps_k <- rep(0, 500)
    mu1_k <- rep(0, 500)
    mu0_k <- rep(0, 500)
    eps_k[1] <- 0.5
    mu1_k[1] <- mean(X[X >= median(X)])
    mu0_k[1] <- mean(X[X <= median(X)])
    # mu1_k[1] <- mean(X[X <= median(X)])
    # mu0_k[1] <- mean(X[X >= median(X)])
    D <- 1
    k <- 1
    while ((D > 10^(-5)) & (k < 500)) {
      pi <- dgamma(X, m, 1/mu1_k[k]) * eps_k[k]/(dgamma(X, m, 1/mu1_k[k]) * eps_k[k] + dgamma(X, m, 1/mu0_k[k]) * (1 - eps_k[k]))
      eps_k[k + 1] <- sum(pi)/n
      mu1_k[k + 1] <- sum(pi * X)/sum(pi * m)
      mu0_k[k + 1] <- sum(X * (1 - pi))/sum(m * (1 - pi))
      D <- (mu1_k[k + 1] - mu1_k[k])^2 + (eps_k[k + 1] - eps_k[k])^2 + (mu0_k[k + 1] - mu0_k[k])^2
      k <- k + 1
    }
    epsi[i] <- eps_k[eps_k > 0][length(eps_k[eps_k > 0])]
    mu_1[i] <- mu1_k[mu1_k > 0][length(mu1_k[mu1_k > 0])]
    mu_0[i] <- mu0_k[mu0_k > 0][length(mu0_k[mu0_k > 0])]
  }
  return(list(p_1 = epsi, 
              p_2 = mu_1,
              p_3 = mu_0
              ))
}

param <- expand.grid(n = c(200,1000),
                     m = c(20,100), 
                     eps = c(0.05,0.1,0.2))
bias <- function(est, praw) abs(mean(est) - praw)
mse <- function(est, praw) mean((est - praw)^2)

wyniki_zad_2 <- list()
for(i in 1:nrow(param)){
  n_i <- param$n[i]
  m_i <- param$m[i]
  eps_i <- param$eps[i]
  wyn <- EM_algorithm(n = n_i, m = m_i, eps = eps_i, k_pow=1000)
  wyniki_zad_2[[i]] <- list(para = param[i, ], 
                            eps = list(bias = bias(wyn$p_1, eps_i),
                                       var = var(wyn$p_1), 
                                       mse = mse(wyn$p_1, eps_i), 
                                       val = wyn$p_1), 
                            mu1 = list(bias = bias(wyn$p_2, 5.5), 
                                       var = var(wyn$p_2), 
                                       mse = mse(wyn$p_2, 5.5), 
                                       val = wyn$p_2),
                            mu0 = list(bias = bias(wyn$p_3, 3), 
                                       var = var(wyn$p_3), 
                                       mse = mse(wyn$p_3, 3), 
                                       val = wyn$p_3)
                            )
}

df_zad_2 <- data.frame(n = integer(), 
                       m = integer(), 
                       eps_praw = numeric(),
                       eps_bias = numeric(), 
                       eps_var = numeric(), 
                       eps_mse = numeric(),
                       mu1_bias = numeric(), 
                       mu1_var = numeric(), 
                       mu1_mse = numeric(),
                       mu0_bias = numeric(),
                       mu0_var = numeric(),
                       mu0_mse = numeric()
                       )

for(i in seq_along(wyniki_zad_2)){
  df_zad_2[i,] <- c(wyniki_zad_2[[i]]$para$n,
                    wyniki_zad_2[[i]]$para$m,
                    wyniki_zad_2[[i]]$para$eps,
                    wyniki_zad_2[[i]]$eps$bias, 
                    wyniki_zad_2[[i]]$eps$var, 
                    wyniki_zad_2[[i]]$eps$mse,
                    wyniki_zad_2[[i]]$mu1$bias, 
                    wyniki_zad_2[[i]]$mu1$var, 
                    wyniki_zad_2[[i]]$mu1$mse,
                    wyniki_zad_2[[i]]$mu0$bias, 
                    wyniki_zad_2[[i]]$mu0$var, 
                    wyniki_zad_2[[i]]$mu0$mse
                    ) 
}

library(dplyr)
library(purrr)

hist_df <- map_dfr(wyniki_zad_2, function(w) {
  data.frame(
    n = w$para$n,
    m = w$para$m,
    eps_praw = w$para$eps,
    eps_hat = w$eps$val,
    mu1_hat = w$mu1$val,
    mu0_hat = w$mu0$val
  )
})

saveRDS(hist_df, "hist_df.rds")

metrics_df <- map_dfr(wyniki_zad_2, function(w) {
  data.frame(
    n = w$para$n,
    m = w$para$m,
    eps_praw = w$para$eps,
    eps_bias = w$eps$bias,
    eps_mse = w$eps$mse,
    eps_var = w$eps$var,
    mu1_bias = w$mu1$bias,
    mu1_mse = w$mu1$mse,
    mu1_var = w$mu1$var,
    mu0_bias = w$mu0$bias,
    mu0_mse = w$mu0$mse,
    mu0_var = w$mu0$var
  )
})

saveRDS(metrics_df, "metrics_df.rds")

ggplot(hist_df, aes(x = eps_hat)) +
  geom_histogram(bins = 30, fill = "lightblue3", color = "blue") +
  facet_grid(n + m ~ eps_praw) +
  labs(
    title = expression("Histogram estymatora " * hat(epsilon)),
    x = expression(hat(epsilon)),
    y = "Częstość"
  ) +
  theme_bw()

ggplot(hist_df, aes(x = mu1_hat)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "darkgreen") +
  facet_grid(n + m ~ eps_praw) +
  labs(
    title = expression("Histogram estymatora " * hat(mu)[1]),
    x = expression(hat(mu)[1]),
    y = "Częstość"
  ) +
  theme_bw()

ggplot(hist_df, aes(x = mu0_hat)) +
  geom_histogram(bins = 30, fill = "tomato", color = "red4") +
  facet_grid(n + m ~ eps_praw) +
  labs(
    title = expression("Histogram estymatora " * hat(mu)[0]),
    x = expression(hat(mu)[1]),
    y = "Częstość"
  ) +
  theme_bw()

library(tidyr)

metrics_long <- metrics_df %>%
  pivot_longer(
    cols = -c(n, m, eps_praw),
    names_to = c("param", "metric"),
    names_sep = "_",
    values_to = "value"
  )

ggplot(
  metrics_long %>% filter(metric == "bias"),
  aes(x = n, y = value, color = param)
) +
  geom_line() +
  geom_point() +
  facet_grid(m ~ eps_praw) +
  labs(
    title = "Bias estymatorów",
    x = "n",
    y = "Bias",
    color = "Parametr"
  ) +
  theme_bw()

ggplot(
  metrics_long %>% filter(metric == "var"),
  aes(x = n, y = value, color = param)
) +
  geom_line() +
  geom_point() +
  facet_grid(m ~ eps_praw) +
  labs(
    title = "Wariancja estymatorów",
    x = "n",
    y = "Wariancja",
    color = "Parametr"
  ) +
  theme_bw()

ggplot(
  metrics_long %>% filter(metric == "mse"),
  aes(x = n, y = value, color = param)
) +
  geom_line() +
  geom_point() +
  facet_grid(m ~ eps_praw) +
  labs(
    title = "MSE estymatorów",
    x = "n",
    y = "MSE",
    color = "Parametr"
  ) +
  theme_bw()

#Zadanie 4.

zad_4 <- function(m, n, e, a, sim){
  m_BH <- matrix(0, nrow = sim, ncol=3)
  m_BFDR <- matrix(0, nrow = sim, ncol=3)
  m_BC <- matrix(0, nrow = sim, ncol=3)
  for (i in 1:sim){
    eps <- rbinom(n, 1, e)
    X <- rgamma(n, m, 1/3) * (1 - eps) + rgamma(n, m, 1/5.5) * eps
    p_w <- 1 - pgamma(X, m, 1/3)
    c_0 <- 1
    c_1 <- 1
    eps_k <- rep(0, 500)
    mu1_k <- rep(0, 500)
    mu0_k <- rep(0, 500)
    eps_k[1] <- 0.5
    mu0_k[1] <- mean(X[X <= median(X)])
    mu1_k[1] <- mean(X[X >= median(X)])
    D <- 1
    k <- 1
    while (D > 10^(-5) & k < 500){
      pi <- dgamma(X, m, 1/mu1_k[k]) * eps_k[k]/(dgamma(X, m, 1/mu1_k[k]) * eps_k[k] + dgamma(X, m, 1/mu0_k[k]) * (1 - eps_k[k]))
      eps_k[k + 1] <- sum(pi)/n
      mu1_k[k + 1] <- sum(pi * X)/sum(pi * m)
      mu0_k[k + 1] <- sum(X * (1 - pi))/sum(m * (1 - pi))
      D <- (mu1_k[k + 1] - mu1_k[k])^2 + (eps_k[k + 1] - eps_k[k])^2 + (mu0_k[k + 1] - mu0_k[k])^2
      k <- k + 1
    }
    eps_pl_in <- eps_k[eps_k > 0][length(eps_k[eps_k > 0])]
    mu1_pl_in <-mu1_k[mu1_k > 0][length(mu1_k[mu1_k > 0])]
    mu0_pl_in <- mu0_k[mu0_k > 0][length(mu0_k[mu0_k > 0])]
    
    #BH
    progBH <- a * (n:1)/((1 - eps_pl_in) * n)
    psort <- sort(p_w, decreasing = TRUE)
    prog_BH <- psort[which.max(psort <= progBH)]
    result_BH <- (p_w <= prog_BH)
    FD_BH <- result_BH & (1 - eps)
    TD_BH <- result_BH & eps
    bl_II_BH <- (1 - result_BH) & eps
    m_BH[i, 1] <- sum(FD_BH)/max(sum(result_BH), 1)
    m_BH[i, 2]  <- sum(TD_BH)/sum(eps)
    m_BH[i, 3] <- c_0 * sum(FD_BH) + c_1 * sum(bl_II_BH)
    
    # BFDR
    BFDR <- function(x) {
      F0 <- pgamma(x, m, 1/mu0_pl_in)
      F1 <- (1 - eps_pl_in) * pgamma(x, m, 1/mu0_pl_in) + eps_pl_in * pgamma(x, m, 1/mu1_pl_in) 
      return((1 - eps_pl_in) * (1 - F0)/(1 - F1) - a)
      }
    if(m == 20){g <- 200} else {g <- 500}
    low <- BFDR(0.0001)
    high <- BFDR(g)
    if (is.nan(low) || is.nan(high) || low * high > 0) {
      cBFDR <- Inf
    } else {
      cBFDR <- uniroot(BFDR, interval = c(0.0001, g))$root
    }
    result_BFDR <- X >= cBFDR
    FD_BFDR <- result_BFDR & (1 - eps)
    TD_BFDR <- result_BFDR & eps
    bl_II_BFDR <- (1 - result_BFDR) & eps
    m_BFDR[i, 1] <- sum(FD_BFDR)/max(sum(result_BFDR), 1)
    m_BFDR[i, 2] <- sum(TD_BFDR)/sum(eps)
    m_BFDR[i, 3] <- c_0 * sum(FD_BFDR) + c_1 * sum(bl_II_BFDR)
    
    #Bayes_Class
    thres_BC <- log((1 - eps_pl_in)/eps_pl_in * (mu1_pl_in/mu0_pl_in)^m)/(1/mu0_pl_in-1/mu1_pl_in)
    result_BC <- X >= thres_BC
    FD_BC <- result_BC & (1 - eps)
    TD_BC <- result_BC & eps
    bl_II_BC <- (1 - result_BC) & eps
    m_BC[i, 1] <- sum(FD_BC)/max(sum(result_BC), 1)
    m_BC[i, 2] <- sum(TD_BC)/sum(eps)
    m_BC[i, 3] <- c_0 * sum(FD_BC) + c_1 * sum(bl_II_BC)
  }
  bh <- apply(m_BH, 2, mean)
  bfdr <- apply(m_BFDR, 2, mean)
  bc <- apply(m_BC, 2, mean)
  return(list(Bh = bh, Bfdr = bfdr, Bc = bc))
}

param_m_20_pl_in <- expand.grid(a = c(0.1, 0.1 * sqrt(10)), n = c(200, 1000), e = c(0.05, 0.1, 0.2))
wyniki_m_20_pl_in <- vector("list", 12)
for(i in seq_len(12)) {
  wyniki_m_20_pl_in[[i]] <- zad_4(m = 20, n = param_m_20_pl_in$n[i], e = param_m_20_pl_in$e[i], a = param_m_20_pl_in$a[i], sim = 1000)
}
names(wyniki_m_20_pl_in) <- apply(param_m_20_pl_in, 1, function(x) paste0("a=", x["a"], "_n=", x["n"], "_e=", x["e"]))

param_m_100_pl_in <- expand.grid(a = c(0.1, 0.1 * sqrt(2)), n = c(200, 1000), e = c(0.05, 0.1, 0.2))
wyniki_m_100_pl_in <- vector("list", 12)
for(i in seq_len(12)) {
  wyniki_m_100_pl_in[[i]] <- zad_4(m = 100, n = param_m_100_pl_in$n[i], e = param_m_100_pl_in$e[i], a = param_m_100_pl_in$a[i], sim = 1000)
}
names(wyniki_m_100_pl_in) <- apply(param_m_100_pl_in, 1, function(x) paste0("a=", x["a"], "_n=", x["n"], "_e=", x["e"]))

zrob_df_z_4 <- function(wynik) {
  df <- rbind(BH = wynik$Bh, BFDR = wynik$Bfdr, BC = wynik$Bc)
  colnames(df) <- c("FDR", "moc", "koszt")
  as.data.frame(df)
}

tabele_m_20_pl_in <- lapply(wyniki_m_20_pl_in, zrob_df_z_4)
tabele_m_100_pl_in <- lapply(wyniki_m_100_pl_in, zrob_df_z_4)

library(dplyr)
library(tidyr)
library(purrr)

df_all_zad4 <- imap_dfr(tabele_m_20_pl_in, function(x, name) {
  
  # wyciągamy parametry z nazwy
  params <- strsplit(name, "_")[[1]]
  a <- as.numeric(sub("a=", "", params[1]))
  n <- as.numeric(sub("n=", "", params[2]))
  e <- as.numeric(sub("e=", "", params[3]))
  
  # zamiana data.frame 3x3 na tidy format
  as.data.frame(x) %>%
    mutate(metryka = rownames(x)) %>%
    pivot_longer(
      cols = -metryka,
      names_to = "wariant",
      values_to = "wartosc"
    ) %>%
    mutate(a = a, n = n, e = e)
})

library(dplyr)

library(dplyr)
library(ggplot2)

#ZAPISANIE JAKO FACTOR IT I LEKKO ZMIANA KOLUMN

df_all_zad4 <- df_all_zad4 %>%
  mutate(
    a = factor(a),
    n = factor(n),
    procedura = factor(metryka, levels = c("BH", "BFDR", "BC")),
    wariant = factor(wariant, levels = c("FDR", "moc", "koszt"))
  )


library(ggplot2)

#WYKRESY DLA FDR DLA KAŻDEGO TESTU ITP

ggplot(
  df_all_zad4 %>% filter(wariant == "FDR"),
  aes(x = e, y = wartosc, color = n, group = n)
) +
  geom_line() +
  geom_point(size = 3) +
  facet_grid(
    procedura ~ a,   
    scales = "free_y"
  ) +
  labs(
    x = expression(epsilon),
    y = "FDR",
    title = "False Discovery Rate (plug-in procedures)",
    color = "n"
  ) +
  theme_bw()

#WYKRESY MOCY

ggplot(
  df_all_zad4 %>% filter(wariant == "moc"),
  aes(x = e, y = wartosc, color = n, group = n)
) +
  geom_line() +
  geom_point(size = 3) +
  facet_grid(
    procedura ~ a, 
    scales = "free_y"
  ) +
  labs(
    x = expression(epsilon),
    y = "FDR",
    title = "False Discovery Rate (plug-in procedures)",
    color = "n"
  ) +
  theme_bw()


#WYKRESY KOSZTU

ggplot(
  df_all_zad4 %>% filter(wariant == "koszt"),
  aes(x = e, y = wartosc, color = n, group = n)
) +
  geom_line() +
  geom_point(size = 3) +
  facet_grid(
    procedura ~ a,  
    scales = "free_y"
  ) +
  labs(
    x = expression(epsilon),
    y = "FDR",
    title = "False Discovery Rate (plug-in procedures)",
    color = "n"
  ) +
  theme_bw()

zrob_df_all <- function(tabele, m_val){
  imap_dfr(tabele, function(x, name) {
    
    params <- strsplit(name, "_")[[1]]
    a <- as.numeric(sub("a=", "", params[1]))
    n <- as.numeric(sub("n=", "", params[2]))
    e <- as.numeric(sub("e=", "", params[3]))
    
    as.data.frame(x) %>%
      mutate(metoda = rownames(x)) %>%
      pivot_longer(
        cols = -metoda,
        names_to = "miara",
        values_to = "wartosc"
      ) %>%
      mutate(
        a = a,
        n = n,
        e = e,
        m = m_val
      )
  })
}

df_plug_in <- bind_rows(
  zrob_df_all(tabele_m_20_pl_in, 20),
  zrob_df_all(tabele_m_100_pl_in, 100)
) %>% mutate(typ = "plug-in")

df_direct <- df_long %>%
  filter(
    Metric %in% c("FDR", "Power", "E"),
    Test %in% c("BH", "BFDR", "BC")
  ) %>%
  mutate(
    miara = recode(
      Metric,
      "Power" = "moc",
      "E"     = "koszt"
    ),
    metoda = Test,
    typ = "direct",
    wartosc = Value,
    e = eps,
    a = alpha
  ) %>%
  select(
    metoda, miara, wartosc,
    a, n, e, m, typ
  )

df_zad4_all <- bind_rows(df_direct, df_plug_in)

saveRDS(df_zad4_all, "df_zad4_all.rds")

plot_zad4_facet <- function(df, miara_nazwa, tytul) {
  
  ggplot(
    df %>% filter(miara == miara_nazwa),
    aes(
      x = e,
      y = wartosc,
      color = typ,
      shape = factor(a),
      group = interaction(typ, a)
    )
  ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.5) +
    facet_grid(
      n + m ~ metoda,
      labeller = labeller(
        n = function(x) paste0("n: ", x),
        m = function(x) paste0("m: ", x)
      )
    ) +
    labs(
      title = tytul,
      x = expression(epsilon),
      y = miara_nazwa,
      color = "Estimation method",
      shape = expression(alpha)
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5)
    )
}

plot_zad4_facet(df_zad4_all, "FDR", "Comparison of FDR: direct vs plug-in")

plot_zad4_facet(df_zad4_all, "moc", "Comparison of power: direct vs plug-in")

plot_zad4_facet(df_zad4_all, "koszt", "Comparison of cost: direct vs plug-in")


