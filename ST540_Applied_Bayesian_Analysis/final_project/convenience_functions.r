library(poisbinom)
library(binGroup)


# Compute Clopper Pearson interval for the mean of success
# Probabilities in a Poisson Binomial distribution
CPInt <- function(x, p.vec, conf.level = 0.95) {
  n <- length(p.vec)
  a <- (1 - conf.level)/2
  if (x == 0) {
    ci <- c(0, 1 - a^(1/n))
  } else if (x == n){
    ci <- c(a^(1/n), 1)
  } else {
    ci <- c(qbeta(a, shape1 = x, shape2 = n-x+1),
            qbeta(1-a, shape1 = x+1, shape2 = n-x))
  }
  ci
}

# Compute Bayesian interval for the mean of success
# Probabilities in a Poisson Binomial distribution
BayesIntNO <- function(x, p.vec, conf.level = 0.95) {
  m <- length(p.vec)
  # order of the hit vector doesn't matter in this
  # case
  hit.vec <- c(rep(1, x), rep(0, m-x))
  # assuming that the probabilities are already ordered
  probs <- p.vec
  int <- vector(length = 2)
  names(int) <- c("LB", "UB")
  
  a <- hit.vec + .01
  b <- 1 - hit.vec + .01
  
  sum.samp <- vector(length = 1000)
  for (j in 1:m) {
    sum.samp <- sum.samp + rbeta(1000, a[j], b[j])
  }
  int[1] <- quantile(sum.samp, probs = .025)
  int[2] <- quantile(sum.samp, probs = .975)
  int/m
}

# Compute Bayesian interval for the mean of success
# Probabilities in a Poisson Binomial distribution
BayesIntBC <- function(x, p.vec, conf.level = 0.95) {
  m <- length(p.vec)
  int <- vector(length = 2)
  names(int) <- c("LB", "UB")
  
  a <- x + .5
  b <- m - x + .5
  
  int[1] <- qbeta(.025, a, b)
  int[2] <- qbeta(.975, a, b)
  
  int
}

# Compute Bayesian interval for the mean of success
# Probabilities in a Poisson Binomial distribution
BayesIntO <- function(x, p.vec, conf.level = 0.95) {
  m <- length(p.vec)
  
  # When the number of possible combinations gets above
  # 20 then take a random sample of size 10
  if (x == 0) {
    hit.mat <- matrix(rep(0, m), nrow = 1)
  } else if (x == m) {
    hit.mat <- matrix(rep(1, m), nrow = 1)
  } else {
    if (m <= 6) {
      col.orders <- combn(1:m, x)
    } else {
      col.orders <- unique(replicate(10, sample(1:m, x)), margin = 2)
    }
    hit.mat <- matrix(nrow = ncol(col.orders), ncol = m)
    for (i in 1:ncol(col.orders)) {
      hit.vec <- vector(length = m)
      hit.vec[col.orders[, i]] <- 1
      hit.vec[-col.orders[, i]] <- 0
      hit.mat[i, ] <- hit.vec 
    }
  }
  int.mat <- matrix(nrow = nrow(hit.mat), ncol = 2)
  colnames(int.mat) <- c("LB", "UB")
  for (i in 1:nrow(hit.mat)) {
    # assuming that the probabilities are already ordered
    probs <- p.vec
    u <- 1:m
    a <- hit.mat[i, ] + .01
    b <- 1 - hit.mat[i, ] + .01 + .0005*u
    
    sum.samp <- vector(length = 1000)
    for (j in 1:m) {
      sum.samp <- sum.samp + rbeta(1000, a[j], b[j])
    }
    int.mat[i, 1] <- quantile(sum.samp, probs = .025)
    int.mat[i, 2] <- quantile(sum.samp, probs = .975)
  }
  int.mat/m
}


# Modified from code by Frank Schaarschmidt
# in the BinWidth Function of the package BinGRoup
ExpWidth <- function (p.vec, conf.level = 0.95,
                      alternative = "two.sided", 
                      method = "CP") {
  n <- length(p.vec)
  L.Ind.bin <- function(x, p.vec, conf.level, alternative, method) {
    switch(method,
           CP = {
      int <- CPInt(x = x, p.vec = p.vec, conf.level = conf.level)
    },  BayesNO = {
      int <- BayesIntNO(x = x, p.vec = p.vec, conf.level = conf.level)
    },  BayesBC = {
      int <- BayesIntBC(x = x, p.vec = p.vec, conf.level = conf.level)
    })
    #Whx were thex subtracting p instead of 0 here?
    # if (alternative == "less") {
    #   CIlength <- int[2] - 0
    # }
    # if (alternative == "greater") {
    #   CIlength <- 1 - int[1]
    # }
    if (alternative == "two.sided") {
      CIlength <- int[2] - int[1]
    }
    list(int, CIlength)
  }
  # Thex are doing this to prevent underflow,
  # do I need to do this with the
  # PoiBin?
  # bin.prob <- function(x, n, p) {
  #   exp(lchoose(n, x) + x * log(p) + (n - x) * log(1 - p))
  # }
  xvec <- 0:n
  Lvec <- numeric(length = length(xvec))
  probvec <- numeric(length = length(xvec))
  intmat <- matrix(nrow = n+1, ncol = 2)
  for (i in 1:length(xvec)) {
    int.res <- L.Ind.bin(x = xvec[i], p.vec = p.vec,
                         conf.level = conf.level, 
                         alternative = alternative, method = method)
    Lvec[i] <- int.res[[2]]
    intmat[i, ] <- int.res[[1]]
    probvec[i] <- dpoisbinom(x = xvec[i], pp = p.vec)
  }
  # print(Lvec)
  # print(probvec)
  expCILength = sum(Lvec * probvec)
  out <- list(expCIWidth = expCILength, intervals = intmat, alternative = alternative, 
              p = p.vec, n = n)
  class(out) <- "IntResults"
  return(out)
}

CovProbs <- function (p.vec, conf.level = 0.95, method = "CP"){
  n <- length(p.vec)
  p.bar <- mean(p.vec)
  
  cov.prob <- 0
  distal.noncov.prob <- 0
  mesial.noncov.prob <- 0
  for (x in 0:n) {
    switch(method,
           CP = {
             int <- CPInt(x = x, p.vec = p.vec, conf.level = conf.level)
           },  BayesNO = {
             int <- BayesIntNO(x = x, p.vec = p.vec, conf.level = conf.level)
           },  BayesBC = {
             int <- BayesIntBC(x = x, p.vec = p.vec, conf.level = conf.level)
           })
    
    # Needed to consider confidence interval as closed only when lower or upper bound
    # was zero or one
    if (int[1] == 0 & int[1] == 1) {
      ind.cov <- T
    } else if (int[1] == 0) {
      ind.cov <- p.bar >= int[1] &&  p.bar < int[2]
    } else if (int[2] == 1) {
      ind.cov <- p.bar > int[1] &&  p.bar <= int[2]
    } else {
      ind.cov <- p.bar > int[1] &&  p.bar < int[2]
    }
    
    ind.cov <- p.bar >= int[1] &&  p.bar <= int[2]
    cov.prob <- cov.prob + ind.cov*dpoisbinom(x = x, pp = p.vec)
    
    ind.distal.noncov <- p.bar > int[2]
    distal.noncov.prob <- distal.noncov.prob + ind.distal.noncov*dpoisbinom(x = x, pp = p.vec)
    ind.mesial.noncov <- p.bar < int[1]
    mesial.noncov.prob <- mesial.noncov.prob + ind.mesial.noncov*dpoisbinom(x = x, pp = p.vec)
  }
  
  data.frame(CP = cov.prob, RNCP = distal.noncov.prob, LNCP = mesial.noncov.prob)
}
