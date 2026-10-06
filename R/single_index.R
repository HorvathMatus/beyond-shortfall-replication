# ============================================================
# Single-Index expectile model for VaR and ES
# Jiang, Hu and Yu (2020), Journal of Financial Econometrics
#
# Port of the authors' code shipped with the previous paper
# ("No Shortfall/expectiele code/SP500.R" and SP500_find_alpha.R),
# with three deliberate changes, all flagged in comments below:
#   (C1) STEP 2 vectorised (the original O(n^2) double loop is far too slow
#        for P = 2099 out-of-sample days) — same estimator, same result.
#   (C2) the regressor matrix x is built the SAME way in-sample and
#        out-of-sample (the original negated the negative parts only in the
#        test block, so train and test features were inconsistent).
#   (C3) the expectile level matched to a quantile level is obtained from
#        the analytic expectile/quantile relation instead of a grid search.
# ============================================================
require(MAVE)
require(MASS)

# Epanechnikov kernel used throughout the original code
.si_kernel <- function(z) 15/16 * (1 - z^2) * (abs(z) <= 1)

# small solver: try solve(), fall back to pseudo-inverse
.si_solve <- function(A, b) {
  out <- tryCatch(solve(A, b), error = function(e) MASS::ginv(A) %*% b)
  if (any(!is.finite(out))) out <- MASS::ginv(A) %*% b
  out
}

# ---- regressor matrix: signed lag-1 and lag-2 returns ----------------------
# x_t = ( y_{t-1}1{y_{t-1}>0}, y_{t-1}1{y_{t-1}<0},
#         y_{t-2}1{y_{t-2}>0}, y_{t-2}1{y_{t-2}<0} )
# (C2) used for both training and prediction
si_design <- function(lag1, lag2) {
  cbind(lag1 * (lag1 > 0), lag1 * (lag1 < 0),
        lag2 * (lag2 > 0), lag2 * (lag2 < 0))
}

# ---- expectile level tau_e that matches quantile level alpha --------------
# (C3) For the tau-expectile e: tau/(1-tau) = E[(e-Y)+] / E[(Y-e)+].
# Setting e equal to the alpha-quantile q and using ES = E[Y | Y<=q] gives
#   tau_e = alpha (q - ES) / ( mu - q + 2 alpha (q - ES) )
si_expectile_level <- function(y, alpha) {
  q  <- as.numeric(quantile(y, alpha, type = 7))
  es <- mean(y[y <= q])
  mu <- mean(y)
  as.numeric(alpha * (q - es) / (mu - q + 2 * alpha * (q - es)))
}

# ---- fit the single-index expectile model ---------------------------------
# y : training responses (length n), X : training regressors (n x p)
# tau : expectile level
si_fit <- function(y, X, tau, max_outer = 4, max_inner = 10, tol = 0.01,
                   verbose = FALSE) {
  n <- length(y); p <- ncol(X)

  # --- initial direction from MAVE, then local linear a_i, b_i ---
  dr    <- MAVE::mave(y ~ X)
  b.mav <- coef(MAVE::mave.dim(dr))[, 1]
  b.mav <- b.mav * sign(b.mav[1])

  loc_init <- function(gam) {
    u  <- as.vector(X %*% gam)
    h  <- 1.06 * n^(-0.2) * sd(u)
    ag <- numeric(n); bg <- numeric(n)
    for (i in 1:n) {
      xi <- u - u[i]
      k  <- .si_kernel(xi / h)
      ab <- .si_solve(crossprod(cbind(k, k * xi), cbind(1, xi)),
                      crossprod(cbind(k, k * xi), y))
      ag[i] <- ab[1]; bg[i] <- ab[2]
    }
    list(u = u, h = h, ag = ag, bg = bg)
  }

  ini <- loc_init(b.mav)
  ag <- ini$ag; bg <- ini$bg
  gam <- b.mav

  # --- alternate: (1) update local a_i,b_i  (2) update direction gamma ---
  outer_it <- 1; eg <- 1
  while (outer_it < max_outer + 1 && eg > tol) {
    gam_old <- gam
    u <- as.vector(X %*% gam)
    h <- 1.06 * n^(-0.2) * sd(u)

    ## STEP 1 — asymmetric-least-squares local linear fit at each point
    for (i in 1:n) {
      xi <- u - u[i]
      k  <- .si_kernel(xi / h)
      z  <- cbind(1, xi)
      it1 <- 1; e1 <- 1
      while (it1 < max_inner && e1 > tol) {
        we  <- abs(tau - ((y - ag[i] - bg[i] * xi) <= 0))
        zkw <- cbind(we * k, we * k * xi)
        abw <- .si_solve(crossprod(zkw, z), crossprod(zkw, y))
        e1  <- sum((c(abw) - c(ag[i], bg[i]))^2)
        ag[i] <- abw[1]; bg[i] <- abw[2]
        it1 <- it1 + 1
      }
    }

    ## STEP 2 — update the index direction gamma
    ## (C1) vectorised version of the original i/j double loop:
    ##   G1 = sum_ij w_ij b_j^2 (x_i-x_j)(x_i-x_j)'
    ##   G2 = sum_ij w_ij b_j   (x_i-x_j)(y_i-a_j)
    ## expanded into matrix products via row/column sums.
    it2 <- 1; e2 <- 1
    while (it2 < max_inner && e2 > tol) {
      u  <- as.vector(X %*% gam)
      h  <- 1.06 * n^(-0.2) * sd(u)
      U  <- outer(u, u, "-")                      # U[i,j] = u_i - u_j
      Km <- .si_kernel(U / h)
      kj <- colSums(Km); kj[kj == 0] <- 1e-12
      Kij <- sweep(Km, 2, kj, "/")
      resid <- outer(y, ag, "-") - sweep(U, 2, bg, "*")
      W  <- abs(tau - (resid <= 0)) * Kij

      C  <- sweep(W, 2, bg^2, "*")
      XtCX <- crossprod(X, C %*% X)
      G1 <- crossprod(X, rowSums(C) * X) - XtCX - t(XtCX) +
            crossprod(X, colSums(C) * X)
      D  <- W * sweep(outer(y, ag, "-"), 2, bg, "*")
      G2 <- crossprod(X, rowSums(D)) - crossprod(X, colSums(D))

      gam_new <- MASS::ginv(G1) %*% G2
      nrm <- sqrt(sum(gam_new^2))
      if (!is.finite(nrm) || nrm == 0) break
      gam_new <- gam_new / nrm * sign(gam_new[1])
      e2  <- sum((gam - gam_new)^2)
      gam <- as.vector(gam_new)
      it2 <- it2 + 1
    }

    eg <- sum((gam - gam_old)^2)
    if (verbose) cat(sprintf("   outer %d: ||dgamma||^2 = %.5f\n", outer_it, eg))
    outer_it <- outer_it + 1
  }

  # --- final local fit at the estimated direction ---
  u <- as.vector(X %*% gam)
  h <- 1.06 * n^(-0.2) * sd(u)
  for (i in 1:n) {
    xi <- u - u[i]
    k  <- .si_kernel(xi / h)
    z  <- cbind(1, xi)
    it <- 1; e <- 1
    while (it < max_inner && e > tol) {
      we  <- abs(tau - ((y - ag[i] - bg[i] * xi) <= 0))
      zkw <- cbind(we * k, we * k * xi)
      abw <- .si_solve(crossprod(zkw, z), crossprod(zkw, y))
      e   <- sum((c(abw) - c(ag[i], bg[i]))^2)
      ag[i] <- abw[1]; bg[i] <- abw[2]
      it <- it + 1
    }
  }

  list(gamma = gam, h = h, y = y, X = X, u = u, ag = ag, bg = bg, tau = tau,
       fitted = ag, hit_rate = mean(y < ag))
}

# ---- one-step-ahead expectile forecasts at new regressor values -----------
si_predict <- function(fit, Xnew, max_inner = 10, tol = 0.01) {
  y <- fit$y; u <- fit$u; h <- fit$h; tau <- fit$tau
  ut <- as.vector(Xnew %*% fit$gamma)
  out <- numeric(length(ut))
  for (i in seq_along(ut)) {
    xi <- u - ut[i]
    k  <- .si_kernel(xi / h)
    z  <- cbind(1, xi); zk <- cbind(k, k * xi)
    ab <- .si_solve(crossprod(zk, z), crossprod(zk, y))
    a  <- ab[1]; b <- ab[2]
    it <- 1; e <- 1
    while (it < max_inner && e > tol) {
      we  <- abs(tau - ((y - a - b * xi) <= 0))
      zkw <- cbind(we * k, we * k * xi)
      abw <- .si_solve(crossprod(zkw, z), crossprod(zkw, y))
      e   <- sum((c(abw) - c(a, b))^2)
      a   <- abw[1]; b <- abw[2]
      it  <- it + 1
    }
    out[i] <- a
  }
  # guard: carry the last sane value forward if the local fit degenerates
  bad <- !is.finite(out)
  if (any(bad)) {
    fallback <- as.numeric(quantile(y, mean(y < fit$ag), type = 7))
    for (i in which(bad)) out[i] <- if (i > 1) out[i-1] else fallback
  }
  out
}

# ---- calibrate the expectile level to the target quantile level -----------
# The original SP500_find_alpha.R grid-searches tau so that the FITTED
# conditional expectile has in-sample exceedance rate equal to alpha
# (tauin <- mean(y < ag); pick tau minimising |tauin - alpha|).
# Same target here, but by bisection started at the analytic value, which
# needs ~10 fits instead of ~20 grid points. The exceedance rate is
# increasing in tau, and is a step function, so we keep the best candidate.
si_calibrate_tau <- function(y, X, alpha, tau0 = NULL, n_iter = 10,
                             verbose = FALSE, ...) {
  if (is.null(tau0)) tau0 <- si_expectile_level(y, alpha)
  lo <- tau0 / 6; hi <- min(tau0 * 6, 0.45)

  best <- list(err = Inf)
  eval_tau <- function(tt) {
    fit <- si_fit(y, X, tt, ...)
    err <- abs(fit$hit_rate - alpha)
    if (verbose) cat(sprintf("    tau=%.5f -> in-sample hit=%.4f (target %.3f)\n",
                             tt, fit$hit_rate, alpha))
    if (err < best$err) best <<- list(err = err, tau = tt, fit = fit)
    fit$hit_rate
  }

  # widen the bracket until it contains the target (small alpha needs very
  # small tau_e, so the initial multiplicative bracket can be too narrow)
  hr_lo <- eval_tau(lo)
  n_exp <- 0
  while (hr_lo > alpha && lo > 1e-7 && n_exp < 6) {
    lo <- lo / 8; hr_lo <- eval_tau(lo); n_exp <- n_exp + 1
  }
  hr_hi <- eval_tau(hi)
  n_exp <- 0
  while (hr_hi < alpha && hi < 0.49 && n_exp < 4) {
    hi <- min(hi * 2, 0.49); hr_hi <- eval_tau(hi); n_exp <- n_exp + 1
  }
  if (hr_lo > alpha || hr_hi < alpha) {          # still outside: keep best
    if (verbose) cat("    target outside bracket, using best endpoint\n")
    return(best)
  }
  for (k in seq_len(n_iter - 2)) {
    mid <- (lo + hi) / 2
    if (eval_tau(mid) < alpha) lo <- mid else hi <- mid
    if (best$err < 1e-4 || (hi - lo) < 1e-5) break
  }
  best
}

# ---- ES implied by an expectile forecast ----------------------------------
# Taylor (2008), eq. (12): with e the tau-expectile matched to quantile alpha,
#   ES_alpha = (1 + tau/((1-2tau) alpha)) e - (tau/((1-2tau) alpha)) mu
# NOTE the original code used (1 - 2*alpha) in place of (1 - 2*tau) and the
# out-of-sample mean for mu (look-ahead). Here: (1-2*tau) and the in-sample mu.
si_es_from_expectile <- function(e, tau, alpha, mu) {
  cf <- tau / ((1 - 2 * tau) * alpha)
  (1 + cf) * e - cf * mu
}
