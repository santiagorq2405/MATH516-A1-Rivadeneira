# =============================================================================
# MATH 516 - Applied Statistics - A1: Bi-Log-Normal Mixture
# Utility functions for density, likelihood, EM, and bootstrap
# =============================================================================

# -- Bi-log-normal mixture density --------------------------------------------

dlnorm_mix <- function(x, pi1, mu1, sigma1, mu2, sigma2, log = FALSE) {
  d1 <- dlnorm(x, meanlog = mu1, sdlog = sigma1)
  d2 <- dlnorm(x, meanlog = mu2, sdlog = sigma2)
  dens <- pi1 * d1 + (1 - pi1) * d2

  if (log) return(log(pmax(dens, .Machine$double.xmin)))
  dens
}

# -- Bi-log-normal mixture CDF ------------------------------------------------

plnorm_mix <- function(q, pi1, mu1, sigma1, mu2, sigma2) {
  pi1 * plnorm(q, meanlog = mu1, sdlog = sigma1) +
    (1 - pi1) * plnorm(q, meanlog = mu2, sdlog = sigma2)
}

# -- Bin probabilities for the mixture ----------------------------------------

bin_probs <- function(startpoints, endpoints, pi1, mu1, sigma1, mu2, sigma2) {
  cdf_upper <- plnorm_mix(endpoints, pi1, mu1, sigma1, mu2, sigma2)
  cdf_lower <- plnorm_mix(startpoints, pi1, mu1, sigma1, mu2, sigma2)
  pmax(cdf_upper - cdf_lower, .Machine$double.xmin)
}

# -- Binned (multinomial) log-likelihood --------------------------------------

loglik_binned <- function(params, startpoints, endpoints, counts) {
  pi1    <- params[1]
  mu1    <- params[2]
  sigma1 <- params[3]
  mu2    <- params[4]
  sigma2 <- params[5]

  if (pi1 <= 0 || pi1 >= 1 || sigma1 <= 0 || sigma2 <= 0) return(-Inf)

  probs <- bin_probs(startpoints, endpoints, pi1, mu1, sigma1, mu2, sigma2)
  sum(counts * log(probs))
}

# Negative version for optim (minimization)
neg_loglik_binned <- function(params, startpoints, endpoints, counts) {
  -loglik_binned(params, startpoints, endpoints, counts)
}

# -- Jittered (continuous) log-likelihood -------------------------------------

loglik_jittered <- function(params, x) {
  pi1    <- params[1]
  mu1    <- params[2]
  sigma1 <- params[3]
  mu2    <- params[4]
  sigma2 <- params[5]

  if (pi1 <= 0 || pi1 >= 1 || sigma1 <= 0 || sigma2 <= 0) return(-Inf)

  sum(dlnorm_mix(x, pi1, mu1, sigma1, mu2, sigma2, log = TRUE))
}

# -- Jittering: generate pseudo-observations from binned data -----------------

jitter_binned <- function(startpoints, endpoints, counts, seed = 42) {
  set.seed(seed)
  x <- mapply(function(a, b, n) {
    if (n == 0) return(numeric(0))
    runif(n, min = a, max = b)
  }, startpoints, endpoints, counts, SIMPLIFY = FALSE)
  unlist(x)
}

# -- EM algorithm for mixture of two log-normals ------------------------------

em_lognormal_mix <- function(x, pi1_init = 0.3, mu1_init = -2, sigma1_init = 0.5,
                             mu2_init = -0.5, sigma2_init = 0.5,
                             max_iter = 500, tol = 1e-8, verbose = FALSE) {
  log_x <- log(x)
  n <- length(x)

  pi1    <- pi1_init
  mu1    <- mu1_init
  sigma1 <- sigma1_init
  mu2    <- mu2_init
  sigma2 <- sigma2_init

  ll_history <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {
    # E-step: compute responsibilities using log-sum-exp for stability
    log_d1 <- dnorm(log_x, mean = mu1, sd = sigma1, log = TRUE) + log(pi1)
    log_d2 <- dnorm(log_x, mean = mu2, sd = sigma2, log = TRUE) + log(1 - pi1)

    log_denom <- pmax(log_d1, log_d2) +
      log(exp(log_d1 - pmax(log_d1, log_d2)) +
            exp(log_d2 - pmax(log_d1, log_d2)))

    gamma1 <- exp(log_d1 - log_denom)

    # M-step
    n1 <- sum(gamma1)
    n2 <- n - n1

    pi1    <- n1 / n
    mu1    <- sum(gamma1 * log_x) / n1
    sigma1 <- sqrt(sum(gamma1 * (log_x - mu1)^2) / n1)
    mu2    <- sum((1 - gamma1) * log_x) / n2
    sigma2 <- sqrt(sum((1 - gamma1) * (log_x - mu2)^2) / n2)

    # Protect against degenerate solutions
    sigma1 <- max(sigma1, 1e-6)
    sigma2 <- max(sigma2, 1e-6)

    ll <- sum(dlnorm_mix(x, pi1, mu1, sigma1, mu2, sigma2, log = TRUE))
    ll_history[iter] <- ll

    if (verbose && iter %% 50 == 0) {
      cat(sprintf("  EM iter %d: loglik = %.4f\n", iter, ll))
    }

    if (iter > 1 && abs(ll_history[iter] - ll_history[iter - 1]) < tol) {
      ll_history <- ll_history[1:iter]
      if (verbose) cat(sprintf("  EM converged at iteration %d\n", iter))
      break
    }
  }

  # Enforce identifiability: component 1 has smaller mu
  if (mu1 > mu2) {
    tmp_mu <- mu1; mu1 <- mu2; mu2 <- tmp_mu
    tmp_s  <- sigma1; sigma1 <- sigma2; sigma2 <- tmp_s
    pi1 <- 1 - pi1
  }

  list(
    pi1 = pi1, mu1 = mu1, sigma1 = sigma1,
    mu2 = mu2, sigma2 = sigma2,
    loglik = ll_history[length(ll_history)],
    ll_history = ll_history,
    n_iter = length(ll_history)
  )
}

# -- Direct MLE on binned likelihood via optim --------------------------------

mle_binned <- function(startpoints, endpoints, counts, start_params,
                       method = "L-BFGS-B") {
  lower <- c(0.01, -5, 0.01, -5, 0.01)
  upper <- c(0.99,  3, 3.0,   3, 3.0)

  result <- optim(
    par     = start_params,
    fn      = neg_loglik_binned,
    startpoints = startpoints,
    endpoints   = endpoints,
    counts      = counts,
    method  = method,
    lower   = lower,
    upper   = upper,
    control = list(maxit = 5000, factr = 1e-10),
    hessian = TRUE
  )

  se <- tryCatch({
    sqrt(diag(solve(result$hessian)))
  }, error = function(e) rep(NA_real_, 5))

  list(
    pi1 = result$par[1], mu1 = result$par[2], sigma1 = result$par[3],
    mu2 = result$par[4], sigma2 = result$par[5],
    loglik = -result$value,
    convergence = result$convergence,
    se = se,
    optim_result = result
  )
}

# -- Parametric bootstrap goodness-of-fit -------------------------------------

chi_sq_stat <- function(observed, expected) {
  valid <- expected > 0
  sum((observed[valid] - expected[valid])^2 / expected[valid])
}

bootstrap_gof <- function(startpoints, endpoints, counts, fitted_params,
                          B = 1000, seed = 42, n_cores = 1) {
  set.seed(seed)
  N <- sum(counts)

  probs_hat <- bin_probs(startpoints, endpoints,
                         fitted_params[1], fitted_params[2], fitted_params[3],
                         fitted_params[4], fitted_params[5])
  probs_hat <- probs_hat / sum(probs_hat)

  expected_hat <- N * probs_hat
  T_obs <- chi_sq_stat(counts, expected_hat)

  run_one_boot <- function(b) {
    sim_counts <- as.numeric(rmultinom(1, size = N, prob = probs_hat))

    fit_b <- tryCatch({
      optim(
        par = fitted_params,
        fn  = neg_loglik_binned,
        startpoints = startpoints,
        endpoints   = endpoints,
        counts      = sim_counts,
        method  = "L-BFGS-B",
        lower   = c(0.01, -5, 0.01, -5, 0.01),
        upper   = c(0.99,  3, 3.0,   3, 3.0),
        control = list(maxit = 2000)
      )
    }, error = function(e) NULL)

    if (is.null(fit_b)) return(NA_real_)

    probs_b <- bin_probs(startpoints, endpoints,
                         fit_b$par[1], fit_b$par[2], fit_b$par[3],
                         fit_b$par[4], fit_b$par[5])
    probs_b <- probs_b / sum(probs_b)
    expected_b <- N * probs_b

    chi_sq_stat(sim_counts, expected_b)
  }

  if (n_cores > 1 && .Platform$OS.type != "windows") {
    T_boot <- parallel::mclapply(seq_len(B), run_one_boot,
                                  mc.cores = n_cores)
    T_boot <- unlist(T_boot)
  } else {
    T_boot <- sapply(seq_len(B), function(b) {
      if (b %% 100 == 0) cat(sprintf("  Bootstrap iteration %d / %d\n", b, B))
      run_one_boot(b)
    })
  }

  T_boot <- T_boot[!is.na(T_boot)]
  p_value <- mean(T_boot >= T_obs)

  list(
    T_obs   = T_obs,
    T_boot  = T_boot,
    p_value = p_value,
    B_valid = length(T_boot)
  )
}

# -- Bayesian MCMC: Random Walk Metropolis-Hastings ---------------------------
# Uses MLE Hessian to calibrate proposals for highly concentrated posteriors

to_unconstrained <- function(theta) {
  c(qlogis(theta[1]), theta[2], log(theta[3]), theta[4], log(theta[5]))
}

to_constrained <- function(phi) {
  c(plogis(phi[1]), phi[2], exp(phi[3]), phi[4], exp(phi[5]))
}

log_posterior <- function(phi, startpoints, endpoints, counts) {
  theta <- to_constrained(phi)
  pi1 <- theta[1]; mu1 <- theta[2]; sigma1 <- theta[3]
  mu2 <- theta[4]; sigma2 <- theta[5]

  ll <- loglik_binned(theta, startpoints, endpoints, counts)
  if (!is.finite(ll)) return(-Inf)

  lp <- dbeta(pi1, 2, 2, log = TRUE) +
    dnorm(mu1, 0, 2, log = TRUE) + dnorm(mu2, 0, 2, log = TRUE) +
    dexp(sigma1, 1, log = TRUE) + dexp(sigma2, 1, log = TRUE)

  lj <- log(pi1 * (1 - pi1)) + log(sigma1) + log(sigma2)

  ll + lp + lj
}

# Estimate proposal covariance from MLE Hessian via numerical differentiation
estimate_proposal_cov <- function(startpoints, endpoints, counts, mle_params) {
  phi_mle <- to_unconstrained(mle_params)
  eps <- 1e-4

  H <- matrix(0, 5, 5)
  f0 <- log_posterior(phi_mle, startpoints, endpoints, counts)

  for (i in 1:5) {
    for (j in i:5) {
      phi_pp <- phi_mle; phi_pp[i] <- phi_pp[i] + eps; phi_pp[j] <- phi_pp[j] + eps
      phi_pm <- phi_mle; phi_pm[i] <- phi_pm[i] + eps; phi_pm[j] <- phi_pm[j] - eps
      phi_mp <- phi_mle; phi_mp[i] <- phi_mp[i] - eps; phi_mp[j] <- phi_mp[j] + eps
      phi_mm <- phi_mle; phi_mm[i] <- phi_mm[i] - eps; phi_mm[j] <- phi_mm[j] - eps

      H[i, j] <- (log_posterior(phi_pp, startpoints, endpoints, counts) -
                     log_posterior(phi_pm, startpoints, endpoints, counts) -
                     log_posterior(phi_mp, startpoints, endpoints, counts) +
                     log_posterior(phi_mm, startpoints, endpoints, counts)) / (4 * eps^2)
      H[j, i] <- H[i, j]
    }
  }

  neg_H_inv <- tryCatch(solve(-H), error = function(e) diag(rep(1e-6, 5)))
  neg_H_inv <- (neg_H_inv + t(neg_H_inv)) / 2

  eig <- eigen(neg_H_inv, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-10)
  eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
}

mcmc_mh <- function(startpoints, endpoints, counts, start_params,
                    proposal_cov, n_iter = 5000, n_warmup = 2000,
                    seed = 42, chain_id = 1) {
  set.seed(seed + chain_id * 7)

  phi_current <- to_unconstrained(start_params)
  n_total <- n_iter + n_warmup
  samples <- matrix(NA, nrow = n_total, ncol = 5)
  log_posts <- numeric(n_total)
  accepted <- 0L
  window_accepted <- 0L

  lp_current <- log_posterior(phi_current, startpoints, endpoints, counts)

  # Cholesky of proposal covariance (scaled by 2.38^2/d for optimal rate)
  scale_factor <- 2.38^2 / 5
  chol_cov <- tryCatch(chol(scale_factor * proposal_cov),
                       error = function(e) diag(sqrt(diag(proposal_cov) * scale_factor)))

  adapt_interval <- 200
  adapt_scale <- 1.0

  for (i in seq_len(n_total)) {
    z <- rnorm(5)
    phi_proposal <- phi_current + adapt_scale * as.numeric(z %*% chol_cov)
    lp_proposal <- log_posterior(phi_proposal, startpoints, endpoints, counts)

    log_alpha <- lp_proposal - lp_current
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      phi_current <- phi_proposal
      lp_current <- lp_proposal
      accepted <- accepted + 1L
      window_accepted <- window_accepted + 1L
    }

    samples[i, ] <- to_constrained(phi_current)
    log_posts[i] <- lp_current

    # Adapt scale during warmup to target ~23% acceptance
    if (i <= n_warmup && i %% adapt_interval == 0) {
      window_rate <- window_accepted / adapt_interval
      if (window_rate < 0.15) {
        adapt_scale <- adapt_scale * 0.7
      } else if (window_rate > 0.35) {
        adapt_scale <- adapt_scale * 1.3
      }
      window_accepted <- 0L
    }
  }

  post_samples <- samples[(n_warmup + 1):n_total, ]
  colnames(post_samples) <- c("pi1", "mu1", "sigma1", "mu2", "sigma2")

  list(
    samples = post_samples,
    log_posterior = log_posts[(n_warmup + 1):n_total],
    acceptance_rate = accepted / n_total,
    chain_id = chain_id
  )
}

run_mcmc_chains <- function(startpoints, endpoints, counts, start_params,
                            n_chains = 4, n_iter = 5000, n_warmup = 2000,
                            seed = 42) {
  cat("  Computing proposal covariance from Hessian...\n")
  prop_cov <- estimate_proposal_cov(startpoints, endpoints, counts, start_params)

  chains <- lapply(seq_len(n_chains), function(ch) {
    perturbed <- start_params * exp(rnorm(5, 0, 0.02))
    perturbed[1] <- pmin(pmax(perturbed[1], 0.05), 0.95)
    cat(sprintf("  Chain %d running...\n", ch))
    mcmc_mh(startpoints, endpoints, counts, perturbed,
            proposal_cov = prop_cov,
            n_iter = n_iter, n_warmup = n_warmup, seed = seed, chain_id = ch)
  })

  all_samples <- do.call(rbind, lapply(chains, function(ch) ch$samples))
  acceptance_rates <- sapply(chains, function(ch) ch$acceptance_rate)

  rhat <- sapply(1:5, function(p) {
    chain_means <- sapply(chains, function(ch) mean(ch$samples[, p]))
    chain_vars  <- sapply(chains, function(ch) var(ch$samples[, p]))
    n <- n_iter; m <- n_chains
    W <- mean(chain_vars)
    B <- n * var(chain_means)
    V_hat <- (n - 1) / n * W + (1 + 1 / m) * B / n
    sqrt(V_hat / W)
  })
  names(rhat) <- c("pi1", "mu1", "sigma1", "mu2", "sigma2")

  cat(sprintf("  Acceptance rates: %s\n",
              paste(round(acceptance_rates, 3), collapse = ", ")))
  cat(sprintf("  Rhat values: %s\n",
              paste(round(rhat, 3), collapse = ", ")))

  list(
    chains = chains,
    all_samples = all_samples,
    acceptance_rates = acceptance_rates,
    rhat = rhat,
    summary = data.frame(
      parameter = c("pi1", "mu1", "sigma1", "mu2", "sigma2"),
      mean = colMeans(all_samples),
      sd = apply(all_samples, 2, sd),
      q2.5 = apply(all_samples, 2, quantile, 0.025),
      q97.5 = apply(all_samples, 2, quantile, 0.975),
      rhat = rhat
    )
  )
}

# -- Single log-normal MLE (baseline comparison) ------------------------------

mle_single_lnorm <- function(startpoints, endpoints, counts) {
  neg_ll <- function(params) {
    mu <- params[1]; sigma <- params[2]
    if (sigma <= 0) return(1e20)
    probs <- plnorm(endpoints, mu, sigma) - plnorm(startpoints, mu, sigma)
    probs <- pmax(probs, .Machine$double.xmin)
    -sum(counts * log(probs))
  }

  result <- optim(c(-0.5, 0.5), neg_ll, method = "L-BFGS-B",
                  lower = c(-5, 0.01), upper = c(3, 3),
                  control = list(maxit = 5000))

  list(mu = result$par[1], sigma = result$par[2], loglik = -result$value)
}
