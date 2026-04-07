#zasada 3 sigm dodać do raportu
library(car)
library(ggplot2)
setwd("C:\\Users\\ameli\\Documents\\R\\statistical learning")

WL <- read.delim("WeightLength.txt")

mu_w <- mean(WL$Weight)
mu_l <- mean(WL$Length)

#możesz zamiast tego colMeans

var_w <- var(WL$Weight)
var_l <- var(WL$Length)
cov_wl <- cov(WL$Weight, WL$Length)

mu <- c(mu_w, mu_l)
sigma <- matrix(c(var_w, cov_wl, cov_wl, var_l), nrow = 2, ncol = 2)

qqnorm(WL$Weight)
qqline(WL$Weight)
qqnorm(WL$Length)
qqline(WL$Length)

WL_mat <- as.matrix(WL)
fun <- mahalanobis(WL, mu, sigma)
alpha_95 <- qchisq(0.95, 2)
alpha_75 <- qchisq(0.75, 2)
score2 <- sum(fun <= alpha_75)
score1 <- sum(fun <= alpha_95 & fun > alpha_75)
score0 <- sum(fun > alpha_95)

ell_95 <- as.data.frame(ellipse(mu, sigma, sqrt(alpha_95)))
ell_75 <- as.data.frame(ellipse(mu, sigma, sqrt(alpha_75)))
colnames(ell_95) <- c("Weight", "Length")
colnames(ell_75) <- c("Weight", "Length")

score <- ifelse(fun <= alpha_75, 2,
                ifelse(fun <= alpha_95, 1, 0))

WL$score <- score

ggplot(WL, aes(Weight, Length, color = factor(score))) +
  geom_point(size = 2) +
  geom_path(data = ell_95, aes(Weight, Length), color = "plum4", linewidth = 1) +
  geom_path(data = ell_75, aes(Weight, Length), color = "plum1", linewidth = 1) +
  scale_color_manual(values = c("0"="skyblue1","1"="skyblue3","2"="skyblue4")) +
  labs(color = "Score") +
  theme_minimal()

P <- eigen(sigma)$vectors
len_P <- apply(P, 2, function(v) sqrt(sum(v^2))) #1 so it is ok

WL_new <- apply(WL_mat, 1, function(x) t(P) %*% x)

WL_new_df <- as.data.frame(t(WL_new))
colnames(WL_new_df) <- c("Z1","Z2")

mean_WL_new <- c(mean(WL_new_df$Z1), mean(WL_new_df$Z2))
sigma_WL_new <- cov(WL_new_df)

ell_95 <- as.data.frame(ellipse(mean_WL_new, sigma_WL_new, sqrt(alpha_95)))
ell_75 <- as.data.frame(ellipse(mean_WL_new, sigma_WL_new, sqrt(alpha_75)))
colnames(ell_95) <- c("Weight", "Length")
colnames(ell_75) <- c("Weight", "Length")
fun <- mahalanobis(WL_new_df, mean_WL_new, sigma_WL_new)

score <- ifelse(fun <= alpha_75, 2,
                ifelse(fun <= alpha_95, 1, 0))

WL_new_df$score <- score

ggplot(WL_new_df, aes(Z1, Z2)) +
  geom_point(size = 2) +
  geom_path(data = ell_95, aes(Weight, Length), color = "plum4", linewidth = 1) +
  geom_path(data = ell_75, aes(Weight, Length), color = "plum1", linewidth = 1) +
  scale_color_manual(values = c("0"="skyblue1","1"="skyblue3","2"="skyblue4")) +
  labs(color = "Score") +
  theme_minimal()

PWL <- read.delim("ParentsWeightLength.txt")

mu_fh <- mean(PWL$FatherHeight)
mu_mh <- mean(PWL$MotherHeight)
mu_w2 <- mean(PWL$Weight)
mu_l2 <- mean(PWL$Length)

mu <- c(mu_fh, mu_mh, mu_w2, mu_l2)
sigma2 <- cov(PWL)

qqnorm(PWL$FatherHeight)
qqline(PWL$FatherHeight)
qqnorm(PWL$MotherHeight)
qqline(PWL$MotherHeight)
qqnorm(PWL$Weight)
qqline(PWL$Weight)
qqnorm(PWL$Length)
qqline(PWL$Length)

child_mu <- c(mu_w2, mu_l2)
parent_mu <- c(mu_fh, mu_mh)
sigma22 <- sigma2[1:2, 1:2]
sigma11 <- sigma2[3:4, 3:4]
sigma21 <- sigma2[1:2, 3:4]
sigma12 <- sigma2[3:4, 1:2]

mu_cond <- t(apply(PWL[, 1:2],1,function(x){
  child_mu + sigma12 %*% solve(sigma22) %*% (x - parent_mu)
}))

sigma_cond <- sigma11 - sigma12 %*% solve(sigma22) %*% sigma21

fun2 <- numeric(nrow(PWL))

for(i in 1:nrow(PWL)){
  fun2[i] <- mahalanobis(PWL[i, 3:4], mu_cond[i,], sigma_cond)
}

ell_95 <- as.data.frame(ellipse(child_mu, sigma_cond, sqrt(alpha_95)))
ell_75 <- as.data.frame(ellipse(child_mu, sigma_cond, sqrt(alpha_75)))
colnames(ell_95) <- c("Weight", "Length")
colnames(ell_75) <- c("Weight", "Length")

score2 <- sum(fun2 <= alpha_75)
score1 <- sum(fun2 <= alpha_95 & fun2 > alpha_75)
score0 <- sum(fun2 > alpha_95)

score <- ifelse(fun2 <= alpha_75,2,
                ifelse(fun2 <= alpha_95,1,0))

PWL_mat <- as.matrix(PWL)
PWL$score <- score

ggplot(PWL, aes(Weight, Length, color=factor(score)))+
  geom_point(size=2)+
  geom_path(data=ell_95, aes(Weight,Length),
            color="purple", linewidth=1)+
  geom_path(data=ell_75, aes(Weight,Length),
            color="pink", linewidth=1)+
  theme_minimal()

mu_spec <- child_mu + sigma12 %*% solve(sigma22) %*% (c(185, 178) - parent_mu)
mu_spec_vector <- as.vector(mu_spec)

ell2_95 <- as.data.frame(ellipse(as.vector(mu_spec), sigma_cond, sqrt(alpha_95)))
ell2_75 <- as.data.frame(ellipse(as.vector(mu_spec), sigma_cond, sqrt(alpha_75)))
colnames(ell2_95) <- c("Weight", "Length")
colnames(ell2_75) <- c("Weight", "Length")

ggplot(PWL, aes(x = Weight, y = Length)) +
  geom_point(alpha = 0.3) + # stare dane w tle
  geom_path(data = ell2_95, color = "red", linewidth = 1) +  # Elipsa 95% dla nowych rodziców
  geom_path(data = ell2_75, color = "blue", linewidth = 1) + # Elipsa 75% dla nowych rodziców
  geom_point(aes(x = mu_spec_vector[1], y = mu_spec_vector[2]), 
             color = "black", shape = 4, size = 5) + # Środek (przewidywany wzrost dziecka)
  labs(title = "Classification Ellipsoids for Parents (185cm, 178cm)",
       subtitle = "Cross marks the predicted (conditional mean) weight and length") +
  theme_minimal()

sigma2_P <- eigen(sigma2)$vectors
len_sigma2_P <- apply(sigma2_P, 2, function(v) sqrt(sum(v^2))) #1 so it is ok


