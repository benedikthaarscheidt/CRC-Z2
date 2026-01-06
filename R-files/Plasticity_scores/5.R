#this files contains the following indices/functions: 
#Coefficient sum (calculate_plasticity),
#Relative Plasticity Index (calculate_RPI) - tested,
#Plasticity Quotient (calculate_PQ) - tested,
#Phenotypic Range (PR) (calculate_PR) - tested,
#Norm of reaction width (calculate_NRW) - tested,
#Environment-Specific Plasticity (ESP) (calculate_ESP) - tested,
#Calculate Plasticity Differential (PD) (calculate_PD) - tested,
#Calculate Fitness Plasticity Index (FPI) (calculate_FPI) - tested,
#Calculate Transplant Plasticity Score (TPS)(calculate_TPS) - tested,
#Calculate Finlay Wilkinson regression across a environmental gradient(TPS)(calculate_finlay_wilkinson) - tested,

#this has 3 functions


########################################

# Function to check and install missing packages
check_and_install_packages = function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
}



############################################


#' Calculate the Plasticity Score from a Polynomial Reaction Norm
#'
#' This function fits polynomial regression models of varying degrees to describe 
#' the relationship between a phenotype and an environmental gradient.
#' The best-fitting polynomial is selected using both AIC and BIC criteria, with 
#' a preference for lower-degree models when performance differences are small.
#'
#' The plasticity score is calculated as the sum of the absolute values of the polynomial coefficients.
#'
#' @param trait_values A numeric vector representing phenotype values.
#' @param env_values (Optional) A numeric vector representing environmental values. 
#' If NULL, equidistant values will be generated.
#' @param max_degree The maximum polynomial degree to consider (default is 3).
#' @param criterion Selection criterion for model complexity. Options are "AIC" (default) or "BIC".
#' @return A list containing:
#'   - `best_degree`: The best-fitting polynomial degree.
#'   - `plasticity_score`: The sum of absolute values of the polynomial coefficients.
#'   - `coefficients`: The estimated polynomial coefficients for the best model.
#'
#' @examples
#' trait_values = c(100, 110, 120, 105, 115, 125)
#' env_values = c(1, 2, 3, 1, 2, 3)
#'
#' # Calculate Plasticity Score with given environmental values
#' plasticity_result = calculate_plasticity(trait_values, env_values)
#' print(plasticity_result)
#'
#' # Calculate Plasticity Score assuming equidistant environments
#' plasticity_result_no_env = calculate_plasticity(trait_values)
#' print(plasticity_result_no_env)
#'@references
#' GENETICS AND EVOLUTION OF PHENOTYPIC PLASTICITY - Samuel M. Scheiner
#' @export
calculate_plasticity = function(trait_values, env_values = NULL, max_degree = 3, criterion = "BIC") {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector")
  }
  
  n_values = length(trait_values)
  
  # If no environmental values are provided, create equidistant values
  if (is.null(env_values)) {
    env_values = seq_len(n_values)
  }
  
  # Ensure environmental values are numeric
  if (!is.numeric(env_values)) {
    stop("env_values must be numeric")
  }
  
  # Check for length mismatch
  if (length(trait_values) != length(env_values)) {
    stop("trait_values and env_values must have the same length")
  }
  
  # Store best model information
  best_degree = 1
  best_criterion_value = Inf
  best_model = NULL
  
  # Try polynomial degrees from 1 to max_degree
  for (degree in 1:max_degree) {
    formula = as.formula(paste("trait_values ~ poly(env_values, ", degree, ", raw = TRUE)", sep = ""))
    model = lm(formula)
    
    # Use AIC or BIC for model selection
    model_criterion = ifelse(criterion == "AIC", AIC(model), BIC(model))
    
    # Penalize higher-degree models if difference is small (<2)
    if ((model_criterion) < best_criterion_value) {
      best_criterion_value = model_criterion
      best_degree = degree
      best_model = model
    }
  }
  
  
  coefficients = coef(best_model)[-1]
  
  # Compute plasticity score as sum of absolute values of coefficients
  plasticity_score = sum(abs(coefficients))
  
  return(list(
    best_degree = best_degree,
    plasticity_score = plasticity_score,
    coefficients = coefficients
  ))
}

##############################################



#' Calculate Cross-Environment Covariance and Correlation for Linear Reaction Norms
#'
#' This function computes the covariance (and optionally correlation) of trait values across different environments.
#' If no environment vector is provided, it assumes equidistant environments.
#'
#' @param trait_values A numeric vector representing trait values measured across different environments.
#' @param env_values (Optional) A numeric vector specifying the environment associated with each trait measurement.
#'                   If NULL, environments are assumed to be equidistant.
#' @param return_correlation Logical; if TRUE, returns both covariance and correlation (default: FALSE).
#'
#' @return A named list containing:
#'   - `covariance`: The covariance between trait values across environments.
#'   - `correlation` (if `return_correlation = TRUE`): The Pearson correlation coefficient.
#'
#' @examples
#' trait_values = c(10, 12, 15, 17, 20, 22, 25, 27, 30, 32)
#'
#' # Calculate covariance with equidistant environments
#' result = cross_env_cov(trait_values)
#' print(result)
#'
#' # Provide an explicit environment vector
#' env_values = seq_along(trait_values)  # Example of environmental values
#' result_with_env = cross_env_cov(trait_values, env_values, return_correlation = TRUE)
#' print(result_with_env)
#'
#' @export
cross_env_cov = function(trait_values, env_values = NULL, return_correlation = FALSE) {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector")
  }
  
  n = length(trait_values)
  
  # If no env_values are provided, assume equidistant environments
  if (is.null(env_values)) {
    env_values = seq_len(n)  # Equidistant environments
  }
  
  # Ensure correct length
  if (length(env_values) != n) {
    stop("trait_values and env_values must have the same length")
  }
  
  # Calculate covariance
  covariance = cov(trait_values, env_values)
  
  # Calculate correlation if requested
  if (return_correlation) {
    correlation = cor(trait_values, env_values)
    return(list(covariance = covariance, correlation = correlation))
  } else {
    return(list(covariance = covariance))
  }
}
##############################################


#' Finlay–Wilkinson Fit with Centered Phenotype-Derived Environment Index
#'
#' Fits genotype-specific Finlay–Wilkinson regressions using the classic
#' phenotype-derived environment index with mandatory centering.
#' Each input row must be an aggregated genotype × environment mean.
#'
#' The model is
#' \deqn{Y_{ij} = M + G_i + \beta_i E_j + e_{ij},}
#' where \eqn{M} is the grand mean, \eqn{G_i = \bar Y_{i\cdot} - M} is the
#' genotype main effect, and \eqn{E_j = \bar Y_{\cdot j} - M} is the centered
#' environment index (environment mean across genotypes minus the grand mean).
#' For each genotype \eqn{i}, \eqn{\beta_i} is estimated by OLS regression
#' of \eqn{Y_{ij}} on \eqn{E_j}.
#'
#' @param y Numeric vector of aggregated trait means \eqn{Y_{ij}} (one per row).
#' @param genotype Genotype identifier (factor or coercible to factor), same length as \code{y}.
#' @param environment Environment identifier (factor or coercible to factor), same length as \code{y}.
#' @param plot Logical; if \code{TRUE}, plots points \eqn{(E_j, Y_{ij})} and fitted FW lines per genotype.
#'
#' @return A \code{data.frame} with one row per genotype containing:
#' \itemize{
#'   \item \code{genotype}: Genotype ID.
#'   \item \code{beta}: Finlay–Wilkinson slope \eqn{\beta_i}.
#'   \item \code{intercept}: OLS intercept; equals \eqn{M + G_i} because \eqn{E_j} is centered.
#'   \item \code{r2}: Coefficient of determination for the per-genotype regression.
#'   \item \code{rmse}: Root-mean-square error (deviation from regression).
#'   \item \code{n_env}: Number of environments used for that genotype.
#'   \item \code{G}: Genotype main effect \eqn{G_i = \bar Y_{i\cdot} - M}.
#'   \item \code{M}: Grand mean \eqn{M}.
#' }
#'
#' @details
#' \itemize{
#'   \item Input must already be aggregated to genotype × environment means.
#'   \item The environment index is \eqn{E_j = \bar Y_{\cdot j} - M}; no leave-one-out adjustment is applied.
#'   \item Centering of \eqn{E_j} does not affect slopes; it makes intercepts interpretable at the average environment.
#'   \item Genotypes with fewer than two usable environments yield \code{NA} for \code{beta}/\code{intercept}.
#' }
#'
#' @examples
#' y  <- c(4,6,5,7,  6,8,7,9)                 # aggregated means
#' g  <- rep(c("G1","G2"), each=4)
#' e  <- rep(c("E1","E2","E3","E4"), times=2)
#' fw_fit(y, g, e, plot=TRUE)
#'
#' @seealso \code{\link{lm}}
#' @export
calculate_finlay_wilkinson <- function(Y, genotype_ids=NULL, env_values=NULL, plot=FALSE) {
  Y <- as.matrix(Y)
  if (is.null(genotype_ids)) {
    gnames <- rownames(Y)
    if (is.null(gnames)) gnames <- paste0("G", seq_len(nrow(Y)))
  } else {
    if (length(genotype_ids) != nrow(Y)) stop("genotype_ids length must match nrow(Y)")
    gnames <- as.character(genotype_ids)
  }
  M <- mean(Y, na.rm=TRUE)
  if (is.null(env_values)) {
    X <- colMeans(Y, na.rm=TRUE) - M
    xlab <- "Environment index E_j (centered)"
  } else {
    if (length(env_values) != ncol(Y)) stop("env_values length must match ncol(Y)")
    X <- as.numeric(env_values) - mean(as.numeric(env_values), na.rm=TRUE)
    xlab <- "Covariate (centered)"
  }
  G_eff <- rowMeans(Y, na.rm=TRUE) - M
  res_list <- lapply(seq_len(nrow(Y)), function(i) {
    y <- as.numeric(Y[i, ])
    ok <- is.finite(y) & is.finite(X)
    x <- X[ok]; yy <- y[ok]
    if (length(yy) < 2 || var(x, na.rm=TRUE) == 0) {
      return(data.frame(genotype=gnames[i], beta=NA_real_, intercept=NA_real_, r2=NA_real_, rmse=NA_real_, n_env=length(yy), G=G_eff[i], M=M))
    }
    m <- lm(yy ~ x)
    pr <- predict(m)
    data.frame(genotype=gnames[i], beta=coef(m)[2], intercept=coef(m)[1], r2=summary(m)$r.squared, rmse=sqrt(mean((yy - pr)^2)), n_env=length(yy), G=G_eff[i], M=M)
  })
  res <- do.call(rbind, res_list); rownames(res) <- res$genotype
  if (plot) {
    plot(NA, xlim=range(X, na.rm=TRUE), ylim=range(Y, na.rm=TRUE), xlab=xlab, ylab="Trait Y_ij", main="Finlay–Wilkinson")
    cols <- setNames(seq_len(nrow(Y)), gnames)
    for (i in seq_len(nrow(Y))) {
      y <- as.numeric(Y[i, ])
      ok <- is.finite(y) & is.finite(X)
      points(X[ok], y[ok], pch=19, col=cols[gnames[i]])
      if (!is.na(res[i, "beta"])) abline(a=res[i, "intercept"], b=res[i, "beta"], col=cols[gnames[i]], lwd=2)
    }
    legend("topleft", legend=gnames, col=seq_len(nrow(Y)), pch=19, lwd=2, bty="n", title="Genotype")
  }
  res[, c("genotype","beta","intercept","r2","rmse","n_env","G","M")]
}

