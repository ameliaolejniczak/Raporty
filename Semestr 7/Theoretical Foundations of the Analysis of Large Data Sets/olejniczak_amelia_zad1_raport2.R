set.seed(411)

#a przed poprawką
calculating_pvalue <- function(sample) {
  n <- length(sample)
  sum_s <- sum(sample)
  1 - ppois(sum_s, n * 5)
}

#a po poprawce
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

mean(many_small_eff[1, ])
mean(many_small_eff[2, ])
