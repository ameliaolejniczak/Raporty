set.seed(1)
n <- 1000
p <- 950
num_rep <- 1000

X <- matrix(rnorm(n * p, sd = 1 / sqrt(1000)), nrow = n, ncol = p)
beta <- c(rep(3, 5), rep(0, 945))
num <- c(5, 10, 20, 100, 500, 950)
epsilon <- rnorm(1000)

simulated_models <- lapply(seq_len(num_rep), function(iter) {
  Y <- X %*% beta + rnorm(1000)
  lapply(num, function(number) {
    lm(Y ~ 0 + X[, seq_len(number)])
  })
})

simulated_summary <- lapply(1:1000, function(x) {
  lapply(simulated_models[[x]], function(sim_mod) {
    summary(sim_mod)
  })
})

library(dplyr)
library(tidyr)
first_betas <- lapply(simulated_summary[[1]], function(s) s$coefficients)
#zadanie1.
#a
betas_info <- data.frame(
  model_size = rep(num, sapply(first_betas, nrow)),
  estimate = unlist(lapply(first_betas, function(mat) mat[, 1])),
  test = unlist(lapply(first_betas, function(mat) mat[, "Pr(>|t|)"]))
)

saveRDS(betas_info, "betas_info.rds")

library(gridExtra)
#b
intervals_1_model <- lapply(simulated_models[[1]], function(sim_mod) confint(sim_mod)[1, 2] - confint(sim_mod)[1, 1])
beta1_info <- data.frame(
  model_size = num,
  std = unlist(lapply(first_betas, function(mat) mat[1, "Std. Error"])),
  interval_length = unlist(intervals_1_model)
)

saveRDS(beta1_info, "beta1_info.rds")

#c
discoveries_one_try <- sapply(simulated_summary[[1]], function(sim_sum) {
  true_hypotheses = 1:5
  rejected_hypotheses = which(sim_sum$coefficients[, 4] < 0.05)
  
  true_discoveries = intersect(rejected_hypotheses, true_hypotheses)
  false_discoveries = setdiff(rejected_hypotheses, true_hypotheses)
  
  n = length(sim_sum$coefficients[, 4])
  rejected_hypotheses_bon = which(sim_sum$coefficients[, 4] < 0.05/n)
  
  true_discoveries_bon = intersect(rejected_hypotheses_bon, true_hypotheses)
  false_discoveries_bon = setdiff(rejected_hypotheses_bon, true_hypotheses)
  
  i = order(sim_sum$coefficients[, 4])
  condition = ((1:length(i))/n)*0.05
  rejected = which(sim_sum$coefficients[i, 4] <= condition)
  i_0 <- if (length(rejected) > 0) max(rejected) else 0
  rejected_hypotheses_bh = if (i_0 > 0) i[1:i_0] else integer(0)
  
  true_discoveries_bh = intersect(rejected_hypotheses_bh, true_hypotheses)
  false_discoveries_bh = setdiff(rejected_hypotheses_bh, true_hypotheses)
  
  
  c(length(true_discoveries), length(false_discoveries), length(true_discoveries_bon), length(false_discoveries_bon), length(true_discoveries_bh), length(false_discoveries_bh))
})

saveRDS(discoveries_one_try, "discoveries_one_try.rds")

discoveries_1_df <- data.frame(
  model_size = num,
  true_discoveries = discoveries_one_try[1,],
  false_discoveries = discoveries_one_try[2,],
  true_discoveries_bon = discoveries_one_try[3,],
  false_discoveries_bon = discoveries_one_try[4,],
  true_discoveries_bh = discoveries_one_try[5, ],
  false_discoveries_bh = discoveries_one_try[6, ]
)

saveRDS(discoveries_1_df, "discoveries.rds")

tables <- lapply(1:6, function(x) {
  tabela_df <- data.frame(
    "Aktualna sytuacja" = c("Prawdziwe", "Fałszywe", "Razem"),
    "Przyjęte" = c(num[x] - 5 - discoveries_one_try[2, x], 5 - discoveries_one_try[1, x], num[x] - 5 - discoveries_one_try[2, x] + 5 - discoveries_one_try[1, x]),
    "Odrzucone" = c(discoveries_one_try[2, x], discoveries_one_try[1, x], discoveries_one_try[2, x] + discoveries_one_try[1, x]),
    "Razem" = c(num[x] - 5, 5, num[x])
  )
  return(tabela_df)
})

tables_bon <- lapply(1:6, function(x) {
  tabela_df <- data.frame(
    "Aktualna sytuacja" = c("Prawdziwe", "Fałszywe", "Razem"),
    "Przyjęte" = c(num[x] - 5 - discoveries_one_try[4, x], 5 - discoveries_one_try[3, x], num[x] - 5 - discoveries_one_try[4, x] + 5 - discoveries_one_try[3, x]),
    "Odrzucone" = c(discoveries_one_try[4, x], discoveries_one_try[3, x], discoveries_one_try[3, x] + discoveries_one_try[4, x]),
    "Razem" = c(num[x] - 5, 5, num[x])
  )
  return(tabela_df)
})

library(knitr)
library(kableExtra)

for (i in 1:6) {
  cat("\n\n### Tabela dla modelu z", num[i], "zmiennymi\n\n")  # Nagłówek dla każdej tabeli
  print(
    kable(tables[[i]], align = "c", booktabs = TRUE, escape = FALSE) %>%
      kable_styling(full_width = FALSE, position = "center")
  )
}

for (i in 1:6) {
  cat("\n\n### Tabela dla modelu z", num[i], "zmiennymi\n\n")  # Nagłówek dla każdej tabeli
  print(
    kable(tables_bon[[i]], align = "c", booktabs = TRUE, escape = FALSE) %>%
      kable_styling(full_width = FALSE, position = "center")
  )
}

library(ggplot2)
library(gridExtra)

p1 <- ggplot(discoveries_1_df, aes(x = model_size, y = true_discoveries_bon)) +
  geom_point(fill = "lightgreen", color = "black", shape = 21, size = 2, stroke = 0.75) +
  labs(title = "Prawdziwe odkrycia",
       x = "Liczba zmiennych",
       y = "Prawdziwe odkrycia") +
  theme_minimal()

p2 <- ggplot(discoveries_1_df, aes(x = model_size, y = false_discoveries_bon)) +
  geom_point(fill = "lightgreen", color = "black", shape = 21, size = 2, stroke = 0.75) +
  labs(title = "Fałszywe odkrycia",
       x = "Liczba zmiennych",
       y = "Fałszywe odkrycia") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)

p1 <- ggplot(discoveries_1_df, aes(x = model_size, y = true_discoveries)) +
  geom_point(fill = "lightgreen", color = "black", shape = 21, size = 2, stroke = 0.75) +
  labs(title = "Prawdziwe odkrycia",
       x = "Liczba zmiennych",
       y = "Prawdziwe odkrycia") +
  theme_minimal()

p2 <- ggplot(discoveries_1_df, aes(x = model_size, y = false_discoveries)) +
  geom_point(fill = "lightgreen", color = "black", shape = 21, size = 2, stroke = 0.75) +
  labs(title = "Fałszywe odkrycia",
       x = "Liczba zmiennych",
       y = "Fałszywe odkrycia") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)

p1 <- ggplot(discoveries_1_df, aes(x = model_size, y = true_discoveries_bh)) +
  geom_point(fill = "lightgreen", color = "black", shape = 21, size = 2, stroke = 0.75) +
  labs(title = "Prawdziwe odkrycia",
       x = "Liczba zmiennych",
       y = "Prawdziwe odkrycia") +
  theme_minimal()

p2 <- ggplot(discoveries_1_df, aes(x = model_size, y = false_discoveries_bh)) +
  geom_point(fill = "lightgreen", color = "black", shape = 21, size = 2, stroke = 0.75) +
  labs(title = "Fałszywe odkrycia",
       x = "Liczba zmiennych",
       y = "Fałszywe odkrycia") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)

#zadanie2
#a i b

simulated_models_unl <- unlist(simulated_models, F, F)
variance <- lapply(simulated_models_unl, function(sim_mod) c(length(sim_mod$coefficients), diag(vcov(sim_mod))[1], confint(sim_mod)[1, 2] - confint(sim_mod)[1, 1]))
variance_df <- data.frame(
  model_size = sapply(variance, `[`, 1),
  variance = sapply(variance, `[`, 2),
  interval = sapply(variance, `[`, 3)
)

saveRDS(variance_df, "variance_df.rds")

library(dplyr)
variance_mean_table <- variance_df %>%
  group_by(model_size) %>%
  summarise(variance = mean(variance), interval = mean(interval))

#Cov(\hat{\beta}) = \sigma^2(\mathbb{X}'\mathbb{X})^{-1} ??

theoretical_results <- sapply(num, function(k) {
  variance_k <- 1000 / (1000 - k - 1)
  interval_k <- 2 * qt(0.975, 1000 - k) * sqrt(variance_k)
  c(variance_k, interval_k)
})

#c
simulated_summary_unlist <- unlist(simulated_summary, F, F)
true_hypotheses = 1:5

results_discoveries = lapply(simulated_summary_unlist, function(model) {
  num_vars_model = length(model$coefficients[, 4])
  
  rejected_hypotheses = unname(which(model$coefficients[, 4] < 0.05))
  true_discoveries = intersect(rejected_hypotheses, true_hypotheses)
  false_discoveries = setdiff(rejected_hypotheses, true_hypotheses)
  
  rejected_hypotheses_bonf = unname(which(model$coefficients[, 4] < (0.05 / num_vars_model)))
  true_discoveries_bonf = intersect(rejected_hypotheses_bonf, true_hypotheses)
  false_discoveries_bonf = setdiff(rejected_hypotheses_bonf, true_hypotheses)
  
  i = unname(order(model$coefficients[, 4]))
  condition = ((1:length(i))/num_vars_model)*0.05
  rejected = which(model$coefficients[i, 4] <= condition)
  i_0 <- if (length(rejected) > 0) max(rejected) else 0
  
  rejected_hypotheses_bh = if (i_0 > 0) i[1:i_0] else integer(0)
  true_discoveries_bh = intersect(rejected_hypotheses_bh, true_hypotheses)
  false_discoveries_bh = setdiff(rejected_hypotheses_bh, true_hypotheses)
  
  list(Correction = c("brak", "Bonferroni", "Benjaminiego-Hochberga"),
       NumVars = rep(num_vars_model, 3),
       NumTrueDiscoveries = c(length(true_discoveries), length(true_discoveries_bonf), length(true_discoveries_bh)),
       NumFalseDiscoveries = c(length(false_discoveries), length(false_discoveries_bonf), length(false_discoveries_bh)),
       FWER = c(length(false_discoveries) > 0, length(false_discoveries_bonf) > 0, length(false_discoveries_bh) > 0),
       FDP = c(length(false_discoveries) / max(length(true_discoveries) + length(false_discoveries), 1),
               length(false_discoveries_bonf) / max(length(true_discoveries_bonf) + length(false_discoveries_bonf), 1),
               length(false_discoveries_bh) / max(length(true_discoveries_bh) + length(false_discoveries_bh), 1)))
})

saveRDS(results_discoveries, "results_discoveries.rds")

library(dplyr)

dplyr::bind_rows(results_discoveries) %>%
  group_by(Correction, NumVars) %>%
  summarize(MeanTrueDiscoveries = mean(NumTrueDiscoveries),
            MeanFalseDiscoveries = mean(NumFalseDiscoveries),
            FWER = mean(FWER),
            FDR = mean(FDP))
