

fit <- function(lambda, y) {
  X <- G $ X
  colnames(X) <- G $ term.names
  Q <- matrix(0,ncol(X), ncol(X))
  idx <- G $ smooth[[1]] $ first.para:G $ smooth[[1]] $ last.para
  Q[idx, idx] <- lambda * G $ S[[1]]
  c(X %*% (solve((t(X) %*% X) + Q) %*% (t(X) %*% y)))
}

dfest <- function(lambda) {
  X <- G $ X
  colnames(X) <- G $ term.names
  Q <- matrix(0,ncol(X), ncol(X))
  idx <- G $ smooth[[1]] $ first.para:G $ smooth[[1]] $ last.para
  Q[idx, idx] <- lambda * G $ S[[1]]
  sum(diag(X %*% (solve((t(X) %*% X) + Q) %*% t(X))))
}



# now lets estimate the degrees of freedom based on simulation
simdf <- function(lambda, nsim = 10000) {
  muest <- fit(1e-15, y = y)
  sig2 <- var(y - muest)

  dfmat <- rep(NA, nsim)
  for (j in 1:nsim) {
    ysim <- rnorm(n, sd = sd) + muest
    dfmat[j] = sum((ysim - muest) * fit(lambda = lambda, y = ysim))/sig2
  }
  df = mean(dfmat)
  hist(dfmat, main = paste("sim:", round(df, 3), " est:", round(dfest(lambda), 2)))
}


# now lets estimate the degrees of freedom based on simulation
sim2df <- function(lambda, nsim = 10000) {
  muest <- mu
  sig2 <- var(y - muest)

  dfmat <- rep(NA, nsim)
  for (j in 1:nsim) {
    ysim <- rnorm(n, sd = sd) + muest
    dfmat[j] = sum((ysim - muest) * fit(lambda = lambda, y = ysim))/sd^2
  }
  df = mean(dfmat)
  hist(dfmat, main = paste("sim2:", round(df, 3), " est:", round(dfest(lambda), 2)))
}




n <- 30
x <- runif(n, 0, 1)
mu <- sin(x * 2 * pi) 
sd <- 1
e <- rnorm(n, sd = sd)
y <- mu + e

G <- mgcv::gam(y ~ s(x, k = 10), fit = FALSE)


par(mfrow = c(2,2))
simdf(exp(-10))
sim2df(exp(-10))

simdf(exp(1))
sim2df(exp(1))








