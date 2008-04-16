## confidence intervals for log-normally distributed variables
lnormCI <- function(x, method = "GCI",
                    alternative = c("two.sided", "less", "greater"),
                    conf.level = 0.95, nsim = 10000){

  alternative <- match.arg(alternative)
  method <- match.arg(method, choices = c("GCI", "PB", "Cox", "MCox"))
  
  ## computations common to all methods
  n <- length(x)
  xbar <- mean(x) 
  s2 <- var(x)
  s <- sd(x)
  
  if (alternative == "two.sided"){
    zarg <- 0.5 + conf.level/2
  } else {
    zarg <- conf.level
  }
  ## Cox and modified Cox
  if (method %in% c("Cox", "MCox")){

    est <- xbar + s2/2 
    z <- if (method == "Cox") qnorm(zarg)
         else qt(zarg, df = n-1)
    rterm <- sqrt(s2/n + (s2^2)/(2*(n-1))) 
    sqrt(1.010 / 40 + 1.010^2 / (2 * 39))
    if (alternative == "two.sided"){
      LL <- exp(est - z * rterm) 
      UL <- exp(est + z * rterm)
    } else if (alternative == "less"){
      LL <- -Inf
      UL <- exp(est + z * rterm)
    } else {
      LL <- exp(est - z * rterm)
      UL <- Inf
    }
    res <-  c(LL, UL)
    # names(res) <- paste(100 * c(1-zarg, zarg), "%", sep = "")
  } else if (method == "GCI"){
  ## Generalized Confidence Intervals
      Z <- rnorm(nsim)
      U2 <- rchisq(nsim, df = n - 1)
      term1 <- Z * (s * sqrt((n-1) / n)) / sqrt(U2)
      term2 <- (0.5 * s2 * (n-1)) / U2
      Ti <- xbar + term1 + term2 # two times +, contrary to article ...
      if (alternative == "two.sided"){
        res <- exp(quantile(Ti, probs = c(1 - zarg, zarg)))
      } else if (alternative == "less"){
        res <- c(-Inf, exp(quantile(Ti, probs = zarg)))
      } else {
        res <- c(exp(quantile(Ti, probs = 1 - zarg)), Inf)
      }
  } else {
  ## Angus (1994) Parametric Bootstrap
    Ni <- rnorm(nsim)
    Xi <- rchisq(nsim, n-1)
    num <- Ni + 0.5 * s * sqrt(n) * (Xi/(n-1) - 1)
    denom <- sqrt(Xi / (n-1) * (1 + 0.5 * s2 * (Xi/(n-1))))
    Ti <- num / denom
    if (alternative == "two.sided"){
      t1t0 <- exp(quantile(Ti, probs = c(1 - zarg, zarg)))
      res <- exp(xbar + 0.5 * s2 - t1t0 * s* sqrt((1 + 0.5 * s2)/n))
    } else if (alternative == "less"){
      t0 <- quantile(Ti, probs = zarg)
      res <- c(-Inf, exp(xbar + 0.5 * s2 - t0 * s* sqrt((1 + 0.5 * s2)/n)))
    } else {
      t1 <- quantile(Ti, probs = 1-zarg)
      res <- c(exp(xbar + 0.5 * s2 - t1 * s* sqrt((1 + 0.5 * s2)/n)), Inf)
    }
  }
  return(zarg)
}
