## Id: BF.simul.R, last updated Oct 3, 2019
## Author: originally coded by Federico Crudu, with contributions of Felipe Osorio

simul.EE <- function(Nsize = 5000, nobs = 500, coef, sigma, alpha = 0.05, trace = TRUE, msg = NULL)
{ ## function to perform the simulation experiment (Section 3 of the manuscript)
  constr.A <- function(coef) {
    coef[2] - 1.0 / coef[3]
  }
  constr.B <- function(coef) {
    coef[2] * coef[3] - 1.0
  }
  deriv.A <- function(coef) {
    c(0.0, 1.0, 1.0 / coef[3]^2)
  }
  deriv.B <- function(coef) {
    c(0.0, coef[3], coef[2])
  }
  moments <- function(coef, ...) {
    res <- y - coef[1] - coef[2] * x1 - exp(coef[3] * x2)
    g2 <- res * x1
    g3 <- res * exp(coef[3] * x2) * x2
    cbind(res, g2, g3)
  }
  moments.deriv <- function(coef, ...) {
    res <- y - coef[1] - coef[2] * x1 - exp(coef[3] * x2)
    g11 <- -1
    g12 <- -x1
    g13 <- -exp(coef[3] * x2) * x2
    G1  <- c(sum(g11), sum(g12), sum(g13))
    g21 <- g11 * x1
    g22 <- g12 * x1
    g23 <- g13 * x1
    G2  <- c(sum(g21), sum(g22), sum(g23))
    g31 <- g11 * exp(coef[3] * x2) * x2
    g32 <- g12 * exp(coef[3] * x2) * x2
    g33 <- g13 * exp(coef[3] * x2) * x2 + res * exp(coef[3] * x2) * x2 * x2
    G3  <- c(sum(g31), sum(g32), sum(g33))
    cbind(G1, G2, G3)
  }
  c.objective <- function(coef) {
    res <- y - coef[1] - coef[2] * x1 - exp(coef[3] * x2)
    g2 <- res * x1
    g3 <- res * exp(coef[3] * x2) * x2
    g  <- cbind(res, g2, g3)
    crossprod(colSums(g))
  }

  if (is.null(msg))
    msg <- "Progress"

  ok <- matrix(FALSE, nrow = Nsize, ncol = 6) # results container
  now <- proc.time()

  if (trace) {
    cat(" ", paste(msg, ":", sep = ""), "\n")
    pb <- txtProgressBar(min = 0, max = Nsize, style = 3)
  }

  # Monte Carlo iterations
  for (i in 1:Nsize) {
    set.seed(123 + i)

    # model building
    coef.true <- coef
    x1  <- rnorm(nobs, 0, sigma)
    x2  <- rnorm(nobs, 0, sigma)
    eps <- rnorm(nobs, 0, sigma)
    y <- coef.true[1] + coef.true[2] * x1 + exp(coef.true[3] * x2) + eps

    # fitting unconstrained model
    X <- cbind(y,x1,x2)
    fm0 <- gmm(moments, X, coef.true)
    cf0 <- fm0$coefficients
    res0 <-fm0$residuals
    V0 <- fm0$vcov

    # fitting unconstrained models
    fm1 <- auglag(par = coef.true, fn = c.objective, heq = constr.A, control.outer = list(trace = FALSE))
    cf1 <- fm1$par
    fm2 <- auglag(par = coef.true, fn = c.objective, heq = constr.B, control.outer = list(trace = FALSE))
    cf2 <- fm2$par

    # distance metric statistic
    g0 <- moments(cf0)
    g1 <- moments(cf1)
    J.hat   <- crossprod(colSums(g0), solve(crossprod(g0), colSums(g0)))
    J.tilde <- crossprod(colSums(g1), solve(crossprod(g1), colSums(g1)))
    DM <- J.tilde - J.hat

    # Wald and BF statistics (hypothesis A)
    g1 <- moments(cf1)
    G1 <- moments.deriv(cf1)
    a1 <- constr.A(cf0)
    a1.dot <- deriv.A(cf0)
    Wald.A <- crossprod(a1, solve(crossprod(a1.dot, V0 %*% a1.dot), a1))
    BF.A <- crossprod(G1 %*% (cf1 - cf0), solve(crossprod(g1), colSums(g1)))

    # Wald, BF and LM statistics (hypothesis B)
    g2 <- moments(cf2)
    G2 <- moments.deriv(cf2)
    a2 <- constr.B(cf0)
    a2.dot <- deriv.B(cf0)
    Wald.B <- crossprod(a2, solve(crossprod(a2.dot, V0 %*% a2.dot), a2))
    BF.B <- crossprod(G2 %*% (cf2 - cf0), solve(crossprod(g2), colSums(g2)))
    LM <- crossprod(colSums(g2), solve(crossprod(g2), colSums(g2)))

    # saving results
    df <- length(a1)
    cutoff.A <- qchisq(1 - alpha, df)
    cutoff.B <- qchisq(1 - alpha, 1)

    ok[i,1] <- Wald.A > cutoff.A
    ok[i,2] <- Wald.B > cutoff.B
    ok[i,3] <- BF.A > cutoff.A
    ok[i,4] <- BF.B > cutoff.B
    ok[i,5] <- LM > cutoff.B
    ok[i,6] <- DM > cutoff.B

    # update progress bar
    if (trace)
      setTxtProgressBar(pb, i)
  }
  if (trace)
    close(pb)

  tnames <- c("Wald.A","Wald.B","BF.A","BF.B","LM","D")
  percentage <- apply(ok, 2, sum) / Nsize
  names(percentage) <- tnames
  cutoff <- c(cutoff.A, cutoff.B)
  names(cutoff) <- c("A","B")
  speed <- proc.time() - now
  out <- list(percentage = 100 * percentage, cutoffs = cutoff, speed = speed)
  out
}

summary.EE <- function(Nsize = 5000, coef, sigma, alpha = 0.05, trace = TRUE) {
  out <- matrix(0, nrow = 4, ncol = 6)
  now <- proc.time()

  out[1,] <- simul.EE(Nsize, nobs = 20,  coef, sigma, alpha, trace, msg = "Stage 1 of 4")$percentage
  out[2,] <- simul.EE(Nsize, nobs = 50,  coef, sigma, alpha, trace, msg = "Stage 2 of 4")$percentage
  out[3,] <- simul.EE(Nsize, nobs = 100, coef, sigma, alpha, trace, msg = "Stage 3 of 4")$percentage
  out[4,] <- simul.EE(Nsize, nobs = 500, coef, sigma, alpha, trace, msg = "Stage 4 of 4")$percentage

  rownames(out) <- c("20","50","100","500")
  colnames(out) <- c("Wald.A","Wald.B","BF.A","BF.B","LM","D")

  list(output = out, speed = proc.time() - now)
}

summary.power <- function(Nsize = 1000, nobs = 20, coef, delta, sigma, alpha = 0.05) {
  n <- length(delta)
  out <- matrix(0, nrow = n, ncol = 6)
  now <- proc.time()

  cat("  Progress:\n")
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n) {
    coef.alt <- c(coef[1], coef[2], coef[3] / delta[i])
    out[i,] <- simul.EE(Nsize, nobs, coef.alt, sigma, alpha, trace = FALSE, msg = NULL)$percentage
    setTxtProgressBar(pb, i)
  }
  close(pb)

  colnames(out) <- c("Wald.A","Wald.B","BF.A","BF.B","LM","D")

  list(power = out, delta = delta, speed = proc.time() - now)
}
