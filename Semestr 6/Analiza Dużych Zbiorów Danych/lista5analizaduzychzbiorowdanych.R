#zadanie( 3. lasso z walidacją krzyżową, najmniejszych kwadratów, gdzie PESEL (ustalimy jeszcze, jak skrócić)
set.seed(56758654)
n <- c(150, 500, 2000)
p <- 100

X <- lapply(n, function(n) {
  f <- matrix(rnorm(n * 3), nrow = n)
  w <- rep((4 - 1:3), each = p) * rnorm(p*3)
  w <- matrix(w, nrow = p)
  e <- matrix(rnorm(n * p), nrow = n)
  
  list(x = f %*% t(w) + e, f = f, w = w)
})

#wyjaśnienie PCA
#chcemy przekształcić dane przy pomocy przekształcenia liniowego aby dostać dane o ostatecznie niższym wymiarze
#mniejszy wymiar, który da nam możliwie dużo informacji
#zapisujemy kolumna po kolumnie
#nowe kolumny są kombinacjami wszystkich kolumn w macierzy danych
#maksymalizujemy kryterium żeby norma była maksymalna (wariancja nowej pierwszej zmiennej największa przy założeniu, że wektor przekształcający ma jednostkową normę)

PCA <- lapply(X, function(x) {
  x_stand <- scale(x[["x"]])
  lambdas <- eigen(t(x_stand) %*% x_stand)$values
  w <- eigen(t(x_stand) %*% x_stand)$vectors
  svd <- svd(x_stand)$v
  t1 <- x_stand %*% w
  t2 <- x_stand %*% svd
  list(t1 = t1, t2 = t2, lambdas = lambdas)
})

corr <- lapply(1:3, function(i) {
  cor_f1 <- cor(PCA[[i]][["t1"]][,1], X[[i]][["f"]][,1])
  cor_f2 <- cor(PCA[[i]][["t1"]][,2], X[[i]][["f"]][,2])
  cor_f3 <- cor(PCA[[i]][["t1"]][,3], X[[i]][["f"]][,3])
  list(cor_f1 = cor_f1,
       cor_f2 = cor_f2,
       cor_f3 = cor_f3)
})

sigma_hat <- lapply(1:3, function(i) {
  n <- n[i]
  lambdas <- PCA[[i]][["lambdas"]]
  lambda_1 <- lambdas[1]
  lambda_2 <- lambdas[2]
  lambda_3 <- lambdas[3]
  1/(n * (p - 3)) * (lambda_1 + lambda_2 + lambda_3)
})

sigma_cor <- replicate(100, lapply(1:3, function(i) {
  f <- matrix(rnorm(n[i] * 3), nrow = n[i])
  w <- rep((4 - 1:3), each = p) * rnorm(p*3)
  w <- matrix(w, nrow = p)
  e <- matrix(rnorm(n[i] * p), nrow = n[i])
  x <- f %*% t(w) + e
  x_stand <- scale(x)
  lambdas <- eigen(t(x_stand) %*% x_stand)$values
  w <- eigen(t(x_stand) %*% x_stand)$vectors
  t1 <- x_stand %*% w
  cor_f1 <- cor(t1[,1], f[,1])
  cor_f2 <- cor(t1[,2], f[,2])
  cor_f3 <- cor(t1[,3], f[,3])
  n <- n[i]
  lambda_1 <- lambdas[1]
  lambda_2 <- lambdas[2]
  lambda_3 <- lambdas[3]
  sigma_hat <- 1/(n * (p - 3)) * (lambda_1 + lambda_2 + lambda_3)
  return(list(cor1 = cor_f1, cor2 = cor_f2, cor3 = cor_f3, sigma_hat = sigma_hat))
}), simplify = FALSE)

cor_n1 <- matrix(NA, nrow = 3, ncol = 100)
cor_n2 <- matrix(NA, nrow = 3, ncol = 100)
cor_n3 <- matrix(NA, nrow = 3, ncol = 100)
sigma <- matrix(NA, nrow = 3, ncol = 100)

for (i in 1:100) {
  cor_n1[,i] <- c(sigma_cor[[i]][[1]][["cor1"]], sigma_cor[[i]][[1]][["cor2"]], sigma_cor[[i]][[1]][["cor3"]])
  cor_n2[,i] <- c(sigma_cor[[i]][[2]][["cor1"]], sigma_cor[[i]][[2]][["cor2"]], sigma_cor[[i]][[2]][["cor3"]])
  cor_n3[,i] <- c(sigma_cor[[i]][[3]][["cor1"]], sigma_cor[[i]][[3]][["cor2"]], sigma_cor[[i]][[3]][["cor3"]])
  sigma[,i] <- c(sigma_cor[[i]][[1]][["sigma_hat"]], sigma_cor[[i]][[2]][["sigma_hat"]], sigma_cor[[i]][[3]][["sigma_hat"]])
}

mean_cor_n1 <- rowMeans(cor_n1)
mean_cor_n2 <- rowMeans(cor_n2)
mean_cor_n3 <- rowMeans(cor_n3)

table_cor_mean <- data.frame(mean_cor_n1 = mean_cor_n1, mean_cor_n2 = mean_cor_n2, mean_cor_n3 = mean_cor_n3)
colnames(table_cor_mean) <- c("n = 150", "n = 500", "n = 2000")
rownames(table_cor_mean) <- c("kolumna 1", "kolumna 2", "kolumna 3")

library(knitr)
kable(table_cor_mean)

library(ggplot2)
library(tidyr)
library(dplyr)

df_cor_n1 <- data.frame(
  cor1 = abs(cor_n1[1,]),
  cor2 = abs(cor_n1[2,]),
  cor3 = abs(cor_n1[3,]),
  n_group = "n1"
)

df_cor_n2 <- data.frame(
  cor1 = abs(cor_n2[1,]),
  cor2 = abs(cor_n2[2,]),
  cor3 = abs(cor_n2[3,]),
  n_group = "n2"
)

df_cor_n3 <- data.frame(
  cor1 = abs(cor_n3[1,]),
  cor2 = abs(cor_n3[2,]),
  cor3 = abs(cor_n3[3,]),
  n_group = "n3"
)

df_all <- bind_rows(df_cor_n1, df_cor_n2, df_cor_n3)

df_long <- pivot_longer(df_all, cols = starts_with("cor"),
                        names_to = "cor_type", values_to = "cor_value")

ggplot(df_long, aes(x = interaction(n_group, cor_type), y = cor_value, fill = cor_type)) +
  geom_boxplot() +
  labs(x = "Grupa i typ korelacji", y = "Wartość korelacji",
       title = "Boxploty korelacji dla różnych grup i typów") +
  scale_x_discrete(labels = c("n1.cor1", "n1.cor2", "n1.cor3",
                              "n2.cor1", "n2.cor2", "n2.cor3",
                              "n3.cor1", "n3.cor2", "n3.cor3")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sigma_data_frame <- data.frame(
  sigma150 = sigma[1, ],
  sigma500 = sigma[2, ],
  sigma2000 = sigma[3, ]
)

colnames(sigma_data_frame) <- c("n = 150", "n = 500", "n = 2000")

sigma_long <- pivot_longer(
  sigma_data_frame,
  cols = everything(),
  names_to = "Sigma",
  values_to = "Value"
)

sigma_long$Sigma <- factor(sigma_long$Sigma, levels = names(sigma_data_frame))

ggplot(sigma_long, aes(x = Sigma, y = Value, fill = Sigma)) +
  geom_boxplot(alpha = 0.7, outlier.color = "red") +
  scale_fill_manual(values = c("lightblue", "orange", "lightgreen")) +
  labs(
    title = "Porównanie rozkładu sigma",
    x = "Rozmiar n",
    y = "Wartość sigma"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#zadanie 2.
#x obserwacje n na p
#3 składniki f, w i p
#e to błąd losowy - przez niego estymatory mają rozne wlasnosci
#dane ktore sa i nas interesuja to f
#sa one liniowo przez w i są zaszumione 

#normalny się przesuwa według właściwego wymiaru 5
#wykładniczy zawyża (powinniśmy się spodziewać przez asymetrię)
#cauchy jest jak trzeba (ciężkie ogony błędów nie psują)

#zadanie 3.
#varclust - robimy k-srednich na kolumnach i w ramach kazdego klastra szukamy prawdziwego wymiaru przy pomocy pesel
#3 a i b zrobilismy na zajeciach, 3c nie trzeba
#3e (tylko pesel i lasso)
#termin oddania 16.06 do 16

n <- 50
p <- 200

data <- replicate(100, {
  f1 <- matrix(rnorm(n * 3), nrow = n)
  w <- rep((4 - 1:3), each = p) * rnorm(p * 3)
  w1 <- matrix(w, nrow = p)
  e1 <- matrix(rnorm(n * p), nrow = n)
  x1 <- f1 %*% t(w1) + e1
  x1 <- scale(x1)
  
  f2 <- matrix(rnorm(n * 3), nrow = n)
  w <- rep((4 - 1:3), each = p) * rnorm(p*3)
  w2 <- matrix(w, nrow = p)
  e2 <- matrix(rnorm(n * p), nrow = n)
  x2 <- f2 %*% t(w2) + e2
  x2 <- scale(x2)
  
  f3 <- matrix(rnorm(n * 3), nrow = n)
  w <- rep((4 - 1:3), each = p) * rnorm(p*3)
  w3 <- matrix(w, nrow = p)
  e3 <- matrix(rnorm(n * p), nrow = n)
  x3 <- f3 %*% t(w3) + e3
  x3 <- scale(x3)
  
  f4 <- matrix(rnorm(n * 3), nrow = n)
  w <- rep((4 - 1:3), each = p) * rnorm(p*3)
  w4 <- matrix(w, nrow = p)
  e4 <- matrix(rnorm(n * p), nrow = n)
  x4 <- f4 %*% t(w4) + e4
  x4 <- scale(x4)
  
  beta <- c(0, 0.5, 0.5, 0, 0.5, 0.5, 0, 0, 0, 0, 0, 0)
  epsilon <- rnorm(n)
  f <- cbind(f1, f2, f3, f4)
  Y <- f %*% beta + epsilon
  
  return(list(X = cbind(x1, x2, x3, x4), f = cbind(f1, f2, f3, f4), Y = Y))
}, simplify = FALSE)

library(pesel)
library(glmnet)

sse <- sapply(data, function(d) {
  X <- d[["X"]]
  Y <- d[["Y"]]
  f <- d[["f"]]
  beta <- c(0, 0.5, 0.5, 0, 0.5, 0.5, 0, 0, 0, 0, 0, 0)
  
  obj_ridge <- cv.glmnet(X, Y, alpha = 0, intercept = FALSE, standardize = FALSE)
  hbeta_ridge <- coef(obj_ridge, s = 'lambda.min')[2:(4*p + 1)]
  sse_ridge <- sum((f %*% beta - X %*% hbeta_ridge)^2)
  
  obj_lasso <- cv.glmnet(X, Y, alpha = 1, intercept = FALSE, standardize = FALSE)
  hbeta_lasso <- coef(obj_lasso, s = 'lambda.min')[2:(4*p + 1)]
  sse_lasso <- sum((f %*% beta - X %*% hbeta_lasso)^2)
  
  w <- 1/abs(hbeta_ridge)
  X1 <- X%*%diag(1/w)
  obj_adlasso <- cv.glmnet(X1, Y, alpha = 1, intercept = FALSE, standardize = FALSE)
  hbeta_adlasso <- coef(obj_adlasso, s = 'lambda.min')[2:(4*p + 1)]
  hbeta_adlasso <- hbeta_adlasso/w
  sse_adlasso <- sum((f %*% beta - X %*% hbeta_adlasso)^2)
  
  p1 <- pesel(X, npc.max=15)$nPCs
  v <- eigen(t(X) %*% X)$vectors
  t <- X %*% v
  X2 <- t[, 1:p1]
  beta_pesel <- solve(t(X2) %*% X2) %*% t(X2) %*% Y
  sse_pesel <- sum((f %*% beta - X2 %*% beta_pesel)^2)
  
  list(sse_ridge = sse_ridge, sse_lasso = sse_lasso, sse_adlasso = sse_adlasso, sse_pesel = sse_pesel)
})

sse_df <- data.frame(
  sse_ridge = sapply(sse, function(x) x$sse_ridge),
  sse_lasso = sapply(sse, function(x) x$sse_lasso),
  sse_adlasso = sapply(sse, function(x) x$sse_adlasso),
  sse_pesel = sapply(sse, function(x) x$sse_pesel)
)

library(tidyr)

sse_long <- sse_df %>% 
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "SSE")

sse_long$Method <- factor(sse_long$Method,
                          levels = c("sse_ridge", "sse_lasso", "sse_adlasso", "sse_pesel"),
                          labels = c("Ridge", "Lasso", "Adaptive Lasso", "PESEL"))

ggplot(sse_long, aes(x = Method, y = SSE, fill = Method)) +
  geom_boxplot() +
  labs(title = "Porównanie metod regresji",
       x = "Metoda",
       y = "Suma kwadratów błędów (SSE)") +
  theme_minimal() +
  theme(legend.position = "none") 