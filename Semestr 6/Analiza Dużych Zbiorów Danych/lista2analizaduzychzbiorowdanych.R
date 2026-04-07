#pakiet pracma randortho moze wygenerowac losowa macierz ortonormalna
#z prawdopodobieństwem gamma mamy wartość z rozkładu normalnego 0, tał^2
#natomiast z prawdopodobieństwem 1 - gamma to będzie wartość 0
#wygenerujemy u jako rozkład jednostajny (0, 1), a potem zmienimy na 0 tam gdzie u >= gamma 
#i na zmienną z rozkładu normalnego poza tym
#nie lm piszemy samodzielnie test na podstawie podpunktu a (bo znamy wariancje)
#sprawdzić różnicę na przekątnej dla X'X oraz identyczności, a następnie maximum
#odległości od 0 dla tych poza przekątną

library(pracma)
X <- randortho(1000, "orthonormal")
#max(diag(t(X)%*%X) - rep(1, 1000))
#max(t(X)%*%X)
#min(t(X)%*%X)

gamma <- c(0.01, 0.05, 0.1)
tau <- c(1.5 * sqrt(2 * log(1000)), 3 * sqrt(2 * log(1000)))

beta_1 <- lapply(gamma, function(g){
  b <- rbinom(1000, size = 1, prob = gamma)
})

for (i in 1:length(beta)) {
  for (j in 1:length(beta_1[[i]])) {
    if (beta_1[[i]][j] == 1) {
      beta_1[[i]][j] <- rnorm(1, 0, tau[1])
    }
  }
}

beta_2 <- lapply(gamma, function(g){
  b <- rbinom(1000, size = 1, prob = gamma)
})

for (i in 1:3) {
  for (j in 1:length(beta_2[[i]])) {
    if (beta_2[[i]][j] == 1) {
      beta_2[[i]][j] <- rnorm(1, 0, tau[2])
    }
  }
}

Y_1 <- lapply(beta_1, function(beta) {
  X %*% beta + rnorm(1000)
})
Y_2 <- lapply(beta_2, function(beta){
  X %*% beta + rnorm(1000)
})

#zadanie 1

#a

beta1ls <- lapply(beta_1, function(beta) {
  beta + t(X) %*% rnorm(1000)
})

beta2ls <- lapply(beta_2, function(beta) {
  beta + t(X) %*% rnorm(1000)
})

#b

#do 0

beta1lsjs10 <- lapply(beta1ls, function(beta) {
  (1 - 998/sum(beta[[1]]**2)) * beta
})

beta2lsjs10 <- lapply(beta2ls, function(beta) {
  (1 - 998/sum(beta[[1]]**2)) * beta
})

#do średniej

beta1lsjsmean <- lapply(beta1ls, function(beta) {
  (1 - 997/999) * beta + 997/999 * mean(beta)
})

beta2lsjsmean <- lapply(beta2ls, function(beta) {
  (1 - 997/999) * beta + 997/999 * mean(beta)
})

#zadanie 2

sum(beta_1[[1]] != 0)

p_vals1 <- lapply(beta1ls, function(beta) {
  2 * (1 - pnorm(abs(beta)))
})

sum(p_vals1[[1]] <= 0.05)

p_vals2 <- lapply(beta2ls, function(beta) {
  2 * (1 - pnorm(abs(beta)))
})

#a

sum(p_vals1[[1]] <= 0.05/1000)

#b

sorted_pvals1 <- lapply(p_vals1, function(pvals) sort(pvals))
sorted_pvals2 <- lapply(p_vals2, function(pvals) sort(pvals))

bh1 <- lapply(sorted_pvals1, function(pvals) {
  ograniczenie <- ((1:1000/1000) * 0.05)
  pvals <= ograniczenie
})

sum(bh1[[1]])

bh2 <- lapply(sorted_pvals2, function(pvals) {
  ograniczenie <- ((1:1000/1000) * 0.05)
  pvals <= ograniczenie
})

#c



#prawdopodobieństwo h0 1-gamma i h1 gamma

