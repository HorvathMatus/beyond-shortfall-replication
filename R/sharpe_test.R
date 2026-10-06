# ============================================================
# Test of equal Sharpe ratios, Ledoit & Wolf (2008), HAC version.
#
# H0: SR_1 = SR_2 for two return series observed on the same dates.
# The difference of the two (unannualised) Sharpe ratios is a smooth function
# of the first two moments of both series; its standard error follows from
# the delta method with a heteroskedasticity- and autocorrelation-consistent
# estimate of the long-run covariance of those moments (prewhitened
# quadratic-spectral kernel, Andrews and Monahan 1992, via sandwich::lrvar).
#
# Returns the difference of the annualised Sharpe ratios and the two-sided
# p-value.
# ============================================================
require(sandwich)

sharpe_diff_test <- function(r1, r2, ppy = 365) {
  ok <- is.finite(r1) & is.finite(r2)
  r1 <- r1[ok]; r2 <- r2[ok]
  m1 <- mean(r1); m2 <- mean(r2)
  g1 <- mean(r1^2); g2 <- mean(r2^2)
  s1 <- g1 - m1^2; s2 <- g2 - m2^2
  if (s1 <= 0 || s2 <= 0) return(c(dSR = NA, p = NA))
  delta <- m1 / sqrt(s1) - m2 / sqrt(s2)
  grad  <- c( g1 / s1^1.5, -g2 / s2^1.5, -m1 / (2 * s1^1.5), m2 / (2 * s2^1.5))
  Y <- cbind(r1 - m1, r2 - m2, r1^2 - g1, r2^2 - g2)
  V <- tryCatch(lrvar(Y, type = "Andrews", prewhite = TRUE),
                error = function(e) cov(Y) / nrow(Y))   # variance of the means
  se <- sqrt(as.numeric(t(grad) %*% V %*% grad))
  c(dSR = delta * sqrt(ppy), p = 2 * pnorm(-abs(delta / se)))
}

# ---- test of equal variances (Ledoit & Wolf 2011 style, HAC) ----------------
# H0: Var(r1) = Var(r2), tested on the log-variance ratio by the delta method
# with the same HAC long-run covariance as above. Returns the ratio of the
# standard deviations and the two-sided p-value.
var_ratio_test <- function(r1, r2) {
  ok <- is.finite(r1) & is.finite(r2)
  r1 <- r1[ok]; r2 <- r2[ok]
  m1 <- mean(r1); m2 <- mean(r2)
  g1 <- mean(r1^2); g2 <- mean(r2^2)
  s1 <- g1 - m1^2; s2 <- g2 - m2^2
  if (s1 <= 0 || s2 <= 0) return(c(sd_ratio = NA, p = NA))
  f    <- log(s1) - log(s2)
  grad <- c(-2 * m1 / s1, 2 * m2 / s2, 1 / s1, -1 / s2)
  Y <- cbind(r1 - m1, r2 - m2, r1^2 - g1, r2^2 - g2)
  V <- tryCatch(lrvar(Y, type = "Andrews", prewhite = TRUE),
                error = function(e) cov(Y) / nrow(Y))
  se <- sqrt(as.numeric(t(grad) %*% V %*% grad))
  c(sd_ratio = sqrt(s1 / s2), p = 2 * pnorm(-abs(f / se)))
}
