# =============================================================================
# Author: Haodong Tian
# Description: Helper functions for MVMR-FA analysis including IV selection,
#              covariance estimation, and MR-SuSiE wrappers.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

library(BayesFM)
library(MendelianRandomization)
library(ggplot2)
library(coda)
library(tidyr)

# Main MVMR_BFA function
#' @title MVMR with Bayesian Factor Analysis
#' @description Implements Multivariate Mendelian Randomization using Bayesian Factor Analysis
#' to handle multiple correlated exposures
#'
#' @param exposure_beta Matrix or data frame of exposure variables
#' @param outcome_beta Vector of outcome betas
#' @param outcome_se Vector of outcome standard errors
#' @param Kmax Maximum number of factors to consider
#' @param model Character specifying either "fixed" or "random" effects model
#' @param Nid Integer specifying minimum number of variables per factor
#' @param burnin Integer specifying number of burnin iterations for MCMC
#' @param iter Integer specifying number of MCMC iterations
#' @param plot_diagnostics Logical, whether to generate diagnostic plots
#' @param generate_report Logical, whether to generate summary report
#' @return Object of class mvmr_bfa containing analysis results
MVMR_BFA <- function(exposure_beta,
                     outcome_beta,
                     outcome_se,
                     Kmax,
                     model        = "fixed",
                     Nid          = 2,
                     burnin       = 500,
                     iter         = 1000,
                     nu0          = Kmax + 1,
                     kappa        = 1 / Kmax,
                     kappa0       = 2,
                     xi0          = 1,
                     seed         = NULL,
                     plot_diagnostics = TRUE,
                     generate_report  = TRUE) {

  # Input validation
  validation <- check_mvmr_bfa_inputs(exposure_beta, outcome_beta, outcome_se)
  if (!validation$valid) {
    stop(paste(validation$messages, collapse = "\n"))
  }
  if (length(validation$messages) > 0) {
    warning(paste(validation$messages, collapse = "\n"))
  }

  if (!is.null(seed)) set.seed(seed)

  # Scale exposure data
  exposure_scaled <- scale(exposure_beta)

  # Perform Bayesian Factor Analysis
  bfa_fit <- befa(exposure_scaled,
                  Nid    = Nid,
                  Kmax   = Kmax,
                  nu0    = nu0,
                  kappa  = kappa,
                  kappa0 = kappa0,
                  xi0    = xi0,
                  burnin = burnin,
                  iter   = iter)

  # Post-process BFA results
  bfa_fit     <- post.column.switch(bfa_fit)
  bfa_fit     <- post.sign.switch(bfa_fit)
  bfa_summary <- summary(bfa_fit)

  # Get number of active factors from final MCMC iteration
  n_factors <- max(bfa_fit$dedic[nrow(bfa_fit$dedic), ])

  # Calculate factor scores using regression method
  p            <- ncol(exposure_scaled)
  alpha_post   <- bfa_fit$alpha
  alpha_means  <- matrix(0, nrow = p, ncol = n_factors)

  # Extract final factor allocation
  final_dedic <- bfa_fit$dedic[nrow(bfa_fit$dedic), ]

  # Extract factor loadings (posterior means)
  for (i in 1:p) {
    k <- final_dedic[i]
    if (k > 0) {  # variable loads on a factor
      alpha_means[i, k] <- mean(alpha_post[, i])
    }
  }

  # Get unique variances
  sigma_post <- colMeans(bfa_fit$sigma)
  sigma_mat  <- diag(sigma_post)

  # Get factor correlation matrix from posterior means
  R_means <- matrix(1, n_factors, n_factors)
  if (n_factors > 1) {
    R_post <- bfa_fit$R
    k <- 1
    for (i in 1:(n_factors - 1)) {
      for (j in (i + 1):n_factors) {
        R_means[i, j] <- R_means[j, i] <- mean(R_post[, k])
        k <- k + 1
      }
    }
  }

  # Calculate regression weights for factor scores
  Lambda <- alpha_means  # factor loadings matrix
  Psi    <- sigma_mat    # unique variances matrix
  R      <- R_means      # factor correlation matrix

  B <- R %*% t(Lambda) %*% solve(Lambda %*% R %*% t(Lambda) + Psi)

  # Calculate factor scores using regression method
  factor_scores <- exposure_scaled %*% t(B)
  colnames(factor_scores) <- paste0("Factor", 1:n_factors)

  # Prepare matrices for MVMR
  bx_matrix   <- as.matrix(factor_scores)
  bxse_matrix <- matrix(0, nrow = nrow(factor_scores), ncol = ncol(factor_scores))

  # Create MVMR input object
  mvmr_input <- mr_mvinput(bx   = bx_matrix,
                            bxse = bxse_matrix,
                            by   = outcome_beta,
                            byse = outcome_se)

  # Perform MVMR analysis
  mvmr_results <- mr_mvivw(model = model, mvmr_input)

  # Calculate diagnostics
  diagnostics <- calculate_diagnostics(bfa_fit, factor_scores)

  # Generate plots if requested
  plots <- NULL
  if (plot_diagnostics) {
    plots <- plot_mvmr_bfa_diagnostics(list(
      bfa_fit     = bfa_fit,
      bfa_summary = bfa_summary,
      factor_scores = factor_scores
    ), diagnostics)
  }

  # Generate report if requested
  report <- NULL
  if (generate_report) {
    report <- generate_mvmr_bfa_report(list(
      factor_scores = factor_scores,
      bfa_summary   = bfa_summary,
      mvmr_results  = list(model = mvmr_results),
      convergence   = list(mh_acceptance = mean(bfa_fit$MHacc))
    ), diagnostics)
  }

  # Return comprehensive results
  results <- list(
    factor_scores = factor_scores,
    bfa_summary   = bfa_summary,
    bfa_fit       = bfa_fit,
    mvmr_results  = list(model = mvmr_results, input = mvmr_input),
    diagnostics   = diagnostics,
    plots         = plots,
    report        = report,
    parameters    = list(Kmax = Kmax, model = model, Nid = Nid,
                         burnin = burnin, iter = iter,
                         nu0 = nu0, kappa = kappa, kappa0 = kappa0, xi0 = xi0)
  )

  class(results) <- "mvmr_bfa"
  return(results)
}

# Print method for mvmr_bfa objects
print.mvmr_bfa <- function(x, ...) {
  cat("\nMVMR-BFA Analysis Results\n")
  cat("=======================\n\n")

  cat("Factor Analysis Summary:\n")
  cat("- Number of factors:", ncol(x$factor_scores), "\n")
  cat("- Number of variables:", nrow(x$bfa_summary$alpha), "\n")
  cat("- MCMC iterations:", x$parameters$iter, "\n")
  cat("- Burn-in:", x$parameters$burnin, "\n\n")

  cat("MVMR Estimates:\n")
  estimates <- x$mvmr_results$model@Estimate
  se        <- x$mvmr_results$model@StdError

  # Two-sided p-values
  z_scores <- estimates / se
  p_values <- 2 * (1 - pnorm(abs(z_scores)))

  for (i in seq_along(estimates)) {
    cat(sprintf("Factor %d: %.4f (SE: %.4f, p-value: %.4e)\n",
                i, estimates[i], se[i], p_values[i]))
  }

  cat("\nDiagnostics Summary:\n")
  cat("- MH acceptance rate:", round(mean(x$bfa_fit$MHacc), 3), "\n")
  cat("- Mean factor dedication stability:",
      round(mean(x$diagnostics$dedication_stability), 3), "\n")
}


#' Check if inputs are valid for MVMR_BFA
#' @param exposure_beta Matrix or dataframe of exposure variables
#' @param outcome_beta Vector of outcome betas
#' @param outcome_se Vector of outcome standard errors
#' @return List with validation results and messages
check_mvmr_bfa_inputs <- function(exposure_beta, outcome_beta, outcome_se) {
  messages <- list()
  is_valid <- TRUE

  # Check data types
  if (!is.matrix(exposure_beta) && !is.data.frame(exposure_beta)) {
    messages <- c(messages, "exposure_beta must be a matrix or data frame")
    is_valid <- FALSE
  }

  # Check dimensions
  if (length(outcome_beta) != length(outcome_se)) {
    messages <- c(messages, "outcome_beta and outcome_se must have same length")
    is_valid <- FALSE
  }

  if (nrow(exposure_beta) != length(outcome_beta)) {
    messages <- c(messages, "Number of rows in exposure_beta must match length of outcome vectors")
    is_valid <- FALSE
  }

  # Check for missing values
  if (any(is.na(exposure_beta))) {
    messages <- c(messages, "Warning: Missing values found in exposure_beta")
  }

  if (any(is.na(outcome_beta)) || any(is.na(outcome_se))) {
    messages <- c(messages, "Warning: Missing values found in outcome data")
  }

  # Check for zero or negative standard errors
  if (any(outcome_se <= 0)) {
    messages <- c(messages, "Error: Standard errors must be positive")
    is_valid <- FALSE
  }

  return(list(valid = is_valid, messages = messages))
}


#' Calculate diagnostic statistics for MVMR_BFA
#' @param bfa_fit BFA fit object
#' @param factor_scores Matrix of factor scores
#' @return List of diagnostic statistics
calculate_diagnostics <- function(bfa_fit, factor_scores) {
  # Factor score metrics
  factor_correlations <- cor(factor_scores)
  factor_variances    <- apply(factor_scores, 2, var)

  # BFA convergence diagnostics
  factor_frequencies   <- table(bfa_fit$nfac) / length(bfa_fit$nfac)
  dedication_stability <- apply(bfa_fit$dedic, 2, function(x) {
    length(unique(x)) / length(x)
  })

  # Effective sample size for key parameters
  ess_alpha <- effectiveSize(mcmc(bfa_fit$alpha))
  ess_sigma <- effectiveSize(mcmc(bfa_fit$sigma))

  # Geweke convergence diagnostics
  geweke_alpha <- geweke.diag(mcmc(bfa_fit$alpha))
  geweke_sigma <- geweke.diag(mcmc(bfa_fit$sigma))

  return(list(
    factor_correlations  = factor_correlations,
    factor_variances     = factor_variances,
    factor_frequencies   = factor_frequencies,
    dedication_stability = dedication_stability,
    effective_sample_size = list(alpha = ess_alpha, sigma = ess_sigma),
    geweke_diagnostics    = list(alpha = geweke_alpha, sigma = geweke_sigma)
  ))
}


#' Plot diagnostics for MVMR_BFA results
#' @param mvmr_bfa_results Results from MVMR_BFA function
#' @param diagnostics Diagnostic statistics from calculate_diagnostics
#' @return List of ggplot objects
plot_mvmr_bfa_diagnostics <- function(mvmr_bfa_results, diagnostics) {
  library(ggplot2)
  library(tidyr)

  plots <- list()

  # Factor loadings heatmap from alpha_means matrix
  p         <- nrow(mvmr_bfa_results$bfa_summary$alpha)
  n_factors <- max(mvmr_bfa_results$bfa_fit$dedic[nrow(mvmr_bfa_results$bfa_fit$dedic), ])

  alpha_means <- matrix(0, nrow = p, ncol = n_factors)
  final_dedic <- mvmr_bfa_results$bfa_fit$dedic[nrow(mvmr_bfa_results$bfa_fit$dedic), ]
  alpha_post  <- mvmr_bfa_results$bfa_fit$alpha

  for (i in 1:p) {
    k <- final_dedic[i]
    if (k > 0) {
      alpha_means[i, k] <- mean(alpha_post[, i])
    }
  }

  clean_names <- gsub("alpha:", "", rownames(mvmr_bfa_results$bfa_summary$alpha))
  rownames(alpha_means) <- clean_names
  colnames(alpha_means) <- paste0("Factor", 1:n_factors)

  loadings_data          <- as.data.frame(alpha_means)
  loadings_data$Variable <- rownames(alpha_means)

  loadings_long <- tidyr::pivot_longer(
    loadings_data, cols = -Variable, names_to = "Factor", values_to = "Loading")

  plots$loadings <- ggplot(loadings_long, aes(x = Factor, y = Variable, fill = Loading)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_minimal() +
    labs(title = "Factor Loadings Heatmap") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Trace plot for number of active factors
  nfac_data <- data.frame(
    Iteration     = 1:length(mvmr_bfa_results$bfa_fit$nfac),
    NumberOfFactors = mvmr_bfa_results$bfa_fit$nfac)

  plots$nfac_trace <- ggplot(nfac_data, aes(x = Iteration, y = NumberOfFactors)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Trace Plot of Number of Factors", x = "Iteration", y = "Number of Factors")

  # Factor score distributions
  factor_scores_df   <- as.data.frame(mvmr_bfa_results$factor_scores)
  factor_scores_long <- tidyr::pivot_longer(
    factor_scores_df, cols = everything(), names_to = "Factor", values_to = "Score")

  plots$factor_dist <- ggplot(factor_scores_long, aes(x = Score, fill = Factor)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Factor) +
    theme_minimal() +
    labs(title = "Factor Score Distributions", x = "Score", y = "Density")

  # MCMC trace plots for representative loadings (first 5)
  alpha_mcmc <- as.data.frame(
    mvmr_bfa_results$bfa_fit$alpha[, 1:min(5, ncol(mvmr_bfa_results$bfa_fit$alpha))])
  colnames(alpha_mcmc) <- gsub("alpha:", "", colnames(alpha_mcmc))
  alpha_long <- tidyr::pivot_longer(
    alpha_mcmc, cols = everything(), names_to = "Loading", values_to = "Value")
  alpha_long$Iteration <- rep(1:nrow(alpha_mcmc), ncol(alpha_mcmc))

  plots$alpha_trace <- ggplot(alpha_long, aes(x = Iteration, y = Value, color = Loading)) +
    geom_line(alpha = 0.5) +
    theme_minimal() +
    labs(title = "MCMC Traces for Factor Loadings", x = "Iteration", y = "Loading Value")

  return(plots)
}


#' Generate summary report for MVMR_BFA analysis
#' @param mvmr_bfa_results Results from MVMR_BFA function
#' @param diagnostics Diagnostic statistics from calculate_diagnostics
#' @return Character string containing formatted report
generate_mvmr_bfa_report <- function(mvmr_bfa_results, diagnostics) {
  report <- ""

  report <- paste0(report, "MVMR-BFA Analysis Report\n", "=======================\n\n")

  report <- paste0(report,
    "Factor Analysis Results:\n", "-------------------------\n",
    "Number of factors identified: ", ncol(mvmr_bfa_results$factor_scores), "\n",
    "Number of variables: ", nrow(mvmr_bfa_results$bfa_summary$alpha), "\n\n")

  report <- paste0(report,
    "Convergence Diagnostics:\n", "------------------------\n",
    "MH acceptance rate: ",
    round(mean(mvmr_bfa_results$convergence$mh_acceptance), 3), "\n",
    "Factor stability: ",
    round(mean(diagnostics$dedication_stability), 3), "\n\n")

  report <- paste0(report, "Factor Loading Summary:\n", "---------------------\n")
  loading_summary <- mvmr_bfa_results$bfa_summary$alpha
  report <- paste0(report, "Mean absolute loadings by factor:\n")
  for (i in 1:ncol(loading_summary)) {
    report <- paste0(report, "Factor ", i, ": ",
                     round(mean(abs(loading_summary[, i])), 3), "\n")
  }
  report <- paste0(report, "\n")

  report <- paste0(report, "MVMR Results:\n", "-------------\n")
  mvmr_estimates <- mvmr_bfa_results$mvmr_results$model@Estimate
  mvmr_se        <- mvmr_bfa_results$mvmr_results$model@StdError

  for (i in seq_along(mvmr_estimates)) {
    report <- paste0(report,
      "Factor ", i, ": ",
      "Estimate = ", round(mvmr_estimates[i], 4),
      " (SE = ", round(mvmr_se[i], 4), ")\n")
  }

  return(report)
}


#' Perform parallel analysis to determine number of factors
#' @param data Matrix or data frame of variables
#' @param n.iter Number of iterations for parallel analysis
#' @param centile Percentile to use for comparison
#' @return Number of factors to retain
parallel_analysis <- function(data, n.iter = 1000, centile = 95) {
  n <- nrow(data)
  p <- ncol(data)

  # Eigenvalues of the observed correlation matrix
  orig.r      <- cor(data)
  orig.values <- eigen(orig.r)$values

  # Generate random data eigenvalues for comparison
  random.values <- matrix(0, nrow = n.iter, ncol = p)
  for (i in 1:n.iter) {
    random.data      <- matrix(rnorm(n * p), nrow = n)
    random.r         <- cor(random.data)
    random.values[i, ] <- eigen(random.r)$values
  }

  # Percentile threshold
  pa.values <- apply(random.values, 2, function(x) quantile(x, centile / 100))

  # Retain factors whose observed eigenvalue exceeds the random threshold
  n.factors <- sum(orig.values > pa.values)

  return(list(n.factors = n.factors, orig.values = orig.values, pa.values = pa.values))
}


#' Factor Analysis for MVMR using ML estimation
#' @param exposure_beta Matrix of exposure variables
#' @param outcome_beta Vector of outcome betas
#' @param outcome_se Vector of outcome standard errors
#' @param n_factor Optional, specific number of factors to extract
#' @param rotation Whether to apply varimax rotation
#' @param parallel_iter Number of iterations for parallel analysis
#' @param parallel_centile Percentile for parallel analysis
#' @param tol Convergence tolerance for ML iterations
#' @return MVMR results with factor analysis
MVMR_FA <- function(exposure_beta,
                    outcome_beta,
                    outcome_se,
                    n_factor         = NULL,
                    rotation         = TRUE,
                    parallel_iter    = 1000,
                    parallel_centile = 95,
                    min.eigval       = 0.7,
                    standardize      = TRUE,
                    maxit            = 100,
                    tol              = 1e-4) {  # tol: stop when change in uniquenesses < tol

  # Input validation
  if (!is.matrix(exposure_beta) && !is.data.frame(exposure_beta)) {
    stop("exposure_beta must be a matrix or data frame")
  }

  if (length(outcome_beta) != length(outcome_se)) {
    stop("outcome_beta and outcome_se must have same length")
  }

  if (nrow(exposure_beta) != length(outcome_beta)) {
    stop("Number of rows in exposure_beta must match length of outcome vectors")
  }

  if (standardize) {
    exposure_scaled <- scale(exposure_beta)
  } else {
    exposure_scaled <- as.matrix(exposure_beta)
  }

  R <- cor(exposure_scaled)
  p <- ncol(R)

  # Determine number of factors via parallel analysis if not specified
  if (is.null(n_factor)) {
    pa_results <- parallel_analysis(exposure_scaled,
                                    n.iter  = parallel_iter,
                                    centile = parallel_centile)
    n_factors  <- pa_results$n.factors

    # Drop factors with eigenvalues below the minimum threshold
    eigen_values <- pa_results$orig.values[1:n_factors]
    if (any(eigen_values < min.eigval)) {
      n_factors <- max(which(eigen_values >= min.eigval))
    }
  } else {
    n_factors  <- n_factor
    pa_results <- list(orig.values = eigen(R)$values, pa.values = rep(NA, ncol(R)))
  }

  # ML Factor Analysis
  # Initial communality estimates using squared multiple correlations
  diag_inv_R <- solve(R)
  h2    <- 1 - 1 / diag(diag_inv_R)
  U2    <- 1 - h2
  Udiag <- diag(U2)

  # Iterative ML estimation
  for (iter in 1:maxit) {
    old_U2 <- U2
    Radj   <- R - Udiag
    eig    <- eigen(Radj)

    loadings <- eig$vectors[, 1:n_factors] %*%
      diag(sqrt(pmax(eig$values[1:n_factors], 0)))

    h2    <- rowSums(loadings^2)
    U2    <- 1 - h2
    Udiag <- diag(U2)

    if (max(abs(U2 - old_U2)) < tol) break  # stop when uniquenesses have converged
  }

  # Varimax rotation if requested
  if (rotation) {
    h2       <- rowSums(loadings^2)
    loadings <- loadings / sqrt(h2)
    rotated  <- loadings

    for (i in 1:50) {  # maximum 50 varimax iterations
      d <- 0
      for (j in 1:(n_factors - 1)) {
        for (k in (j + 1):n_factors) {
          u <- rotated[, j]; v <- rotated[, k]
          a <- sum(u^2 - v^2); b <- 2 * sum(u * v)
          if (abs(b) > 1e-10) {
            theta       <- atan2(b, a) / 4
            c           <- cos(theta); s <- sin(theta)
            temp        <- rotated[, j]
            rotated[, j] <- c * u - s * v
            rotated[, k] <- s * temp + c * v
            d <- d + abs(theta)
          }
        }
      }
      if (d < 1e-6) break  # varimax convergence criterion
    }

    loadings <- rotated * sqrt(h2)
  }

  # Residual correlation and model fit (chi-square)
  residual  <- R - (loadings %*% t(loadings) + Udiag)
  chi_square <- (nrow(R) - 1 - (2 * p + 5) / 6 - (2 * n_factors) / 3) * sum(residual^2)
  df        <- ((p - n_factors)^2 - p - n_factors) / 2
  p_value   <- 1 - pchisq(chi_square, df)

  ml_results <- list(
    uniquenesses  = U2, communalities = h2, residual = residual,
    converged     = iter < maxit, iterations = iter,
    chi_square    = chi_square, df = df, p_value = p_value)

  # Factor scores via regression method
  weights <- solve(R) %*% loadings
  scores  <- exposure_scaled %*% weights

  bx_matrix   <- as.matrix(scores)
  bxse_matrix <- matrix(0, nrow = nrow(scores), ncol = ncol(scores))

  mvmr_input   <- mr_mvinput(bx = bx_matrix, bxse = bxse_matrix,
                              by = outcome_beta, byse = outcome_se)
  mvmr_results <- mr_mvivw(mvmr_input)

  results <- list(
    factor_analysis = list(
      loadings          = loadings, factor_scores = scores,
      n_factors         = n_factors, parallel_analysis = pa_results,
      ml_results        = ml_results, rotated = rotation,
      variable_names    = names(exposure_beta)),
    mvmr_results = list(model = mvmr_results, input = mvmr_input),
    parameters   = list(
      rotation         = rotation, parallel_iter = parallel_iter,
      parallel_centile = parallel_centile, min.eigval = min.eigval,
      standardize      = standardize, maxit = maxit, tol = tol,
      specified_n_factor = !is.null(n_factor))
  )

  class(results) <- "mvmr_fa"
  return(results)
}


# Print method for mvmr_fa objects
print.mvmr_fa <- function(x, ..., digits = 3) {
  cat("\nMVMR-FA Analysis Results (Maximum Likelihood)\n")
  cat("=======================================\n\n")

  cat("Factor Analysis Summary:\n")
  cat("- Number of factors:", x$factor_analysis$n_factors, "\n")
  cat("- Number of variables:", nrow(x$factor_analysis$loadings), "\n")
  cat("- Rotation applied:", x$parameters$rotation, "\n")
  cat("- ML iterations:", x$factor_analysis$ml_results$iterations, "\n\n")

  cat("Model Fit:\n")
  cat("Chi-square test of model fit:\n")
  cat("- H0: The factor model explains the correlations\n")
  cat("- H1: The factor model does not explain the correlations\n")
  cat("- Chi-square:", round(x$factor_analysis$ml_results$chi_square, 3), "\n")
  cat("- Degrees of freedom:", x$factor_analysis$ml_results$df, "\n")
  cat("- P-value:", format.pval(x$factor_analysis$ml_results$p_value), "\n")
  cat("Note: P-value > 0.05 suggests good model fit\n\n")

  cat("Factor Loadings Matrix:\n")
  loadings <- x$factor_analysis$loadings
  if (!is.null(x$factor_analysis$variable_names)) {
    rownames(loadings) <- x$factor_analysis$variable_names
  }
  colnames(loadings) <- paste0("Factor", 1:ncol(loadings))
  print(round(loadings, digits))
  cat("\n")

  cat("Communalities:\n")
  communalities       <- x$factor_analysis$ml_results$communalities
  names(communalities) <- x$factor_analysis$variable_names
  print(round(communalities, digits))
  cat("\n")

  cat("MVMR Estimates:\n")
  estimates <- x$mvmr_results$model@Estimate
  se        <- x$mvmr_results$model@StdError
  z_scores  <- estimates / se
  p_values  <- 2 * (1 - pnorm(abs(z_scores)))
  for (i in seq_along(estimates)) {
    cat(sprintf("Factor %d: %.4f (SE: %.4f, p-value: %.4e)\n",
                i, estimates[i], se[i], p_values[i]))
  }
}


# Plot method for mvmr_fa objects
plot.mvmr_fa <- function(x, exposure_names = NULL, ...) {
  plots <- list()

  # Scree plot with parallel analysis threshold
  scree_data <- data.frame(
    Component  = 1:length(x$factor_analysis$parallel_analysis$orig.values),
    Eigenvalue = x$factor_analysis$parallel_analysis$orig.values,
    PA_Threshold = x$factor_analysis$parallel_analysis$pa.values)

  plots$scree <- ggplot(scree_data, aes(x = Component)) +
    geom_line(aes(y = Eigenvalue,    color = "Observed"),          size = 1) +
    geom_line(aes(y = PA_Threshold,  color = "Parallel Analysis"), size = 1, linetype = "dashed") +
    geom_vline(xintercept = x$factor_analysis$n_factors, linetype = "dotted") +
    theme_minimal() +
    labs(title = "Scree Plot with Parallel Analysis", y = "Eigenvalue") +
    scale_color_manual(values = c("Observed" = "blue", "Parallel Analysis" = "red")) +
    theme(legend.title = element_blank())

  # Factor loadings heatmap
  loadings_data          <- as.data.frame(x$factor_analysis$loadings)
  loadings_data$Variable <- x$factor_analysis$variable_names
  loadings_long <- tidyr::pivot_longer(
    loadings_data, cols = -Variable, names_to = "Factor", values_to = "Loading")

  plots$loadings <- ggplot(loadings_long, aes(x = Factor, y = Variable, fill = Loading)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_minimal() +
    labs(title = "Factor Loadings Heatmap") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8))

  # Factor scores distribution
  scores_data <- as.data.frame(x$factor_analysis$factor_scores)
  scores_long <- tidyr::pivot_longer(
    scores_data, cols = everything(), names_to = "Factor", values_to = "Score")

  plots$scores <- ggplot(scores_long, aes(x = Score, fill = Factor)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Factor) +
    theme_minimal() +
    labs(title = "Factor Scores Distribution")

  return(plots)
}
