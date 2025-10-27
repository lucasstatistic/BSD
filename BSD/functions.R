#--------------------------------------------------------------------
# Bivariate Simplex/Beta Distribution Copula (BSBC) Negative Log-Likelihood
#
# Calculates the negative log-likelihood of the Bivariate Simplex or Beta
# distribution with a specified copula structure. This function is typically
# used as the objective function for maximum likelihood estimation (MLE).
#
densitBSB <- function(par, y, dist = c("Simplex", "Beta"),
                      copula = c("FGM","Clayton","Frank","Gaussian", "Gumbel")) {
  
  # Distribution validation with friendly message
  valid_dists <- c("Simplex", "Beta")
  if (!dist %in% valid_dists) {
    message("Invalid distribution: '", dist, "'")
    message("Valid distributions: ", paste(valid_dists, collapse = ", "))
    return(1e10) # Return high value for optimization failure
  }
  
  # Copula validation with friendly message
  valid_copulas <- c("FGM","Clayton","Frank","Gaussian","Gumbel")
  if (!copula %in% valid_copulas) {
    message("Invalid copula: '", copula, "'")
    message("Valid copulas: ", paste(valid_copulas, collapse = ", "))
    return(1e10) # Return high value for optimization failure
  }
  
  # Function bodies for dSIMPLEX, pSIMPLEX, dBE, pBE are assumed to be available
  # (e.g., loaded from gamlss.dist or provided internally)
  
  suppressWarnings({
    mu1    <- plogis(par[1])
    mu2    <- plogis(par[2])
    
    sigma1 <- exp(par[3])
    sigma2 <- exp(par[4])
    
    # Check for invalid means (must be in (0, 1))
    if (any(is.na(mu1) | is.na(mu2) | mu1 <= 0 | mu1 >= 1 | mu2 <= 0 | mu2 >= 1)) {
      return(1e6)
    }
    
    if (copula %in% c("FGM","Gaussian")) {
      lambda <- tanh(par[5])    # restricts to (-1,1)
    } else {
      lambda <- par[5]
    }
    
    y1 <- y[,1]
    y2 <- y[,2]
    
    # --- Marginal Densities and CDFs ---
    if(dist == "Simplex") {
      # Check if Simplex functions are available
      if (!exists("dSIMPLEX") || !exists("pSIMPLEX")) {
        message("Simplex distribution functions (dSIMPLEX, pSIMPLEX) not found. Please load gamlss.dist package")
        return(1e10)
      }
      d1 <- dSIMPLEX(y1, mu1, sigma1, log = TRUE)
      d2 <- dSIMPLEX(y2, mu2, sigma2, log = TRUE)
      u  <- pSIMPLEX(y1, mu1, sigma1)
      v  <- pSIMPLEX(y2, mu2, sigma2)
      
    } else if(dist == "Beta") {
      # Check if Beta functions are available
      if (!exists("dBE") || !exists("pBE")) {
        message("Beta distribution functions (dBE, pBE) not found. Please load gamlss.dist package")
        return(1e10)
      }
      d1 <- dBE(y1, mu1, sigma1, log = TRUE)
      d2 <- dBE(y2, mu2, sigma2, log = TRUE)
      u  <- pBE(y1, mu1, sigma1)
      v  <- pBE(y2, mu2, sigma2)
    }
    
    # Clamp CDFs to avoid numerical issues at boundaries
    u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
    v <- pmin(pmax(v, 1e-6), 1 - 1e-6)
    
    # --- Copula Density (Log-scale) ---
    if (copula == "FGM") {
      dens_cop <- log1p(lambda * (1 - 2*u) * (1 - 2*v))
      
    } else if (copula == "Clayton") {
      dens_cop <- log((1+lambda) * u^(-1-lambda) * v^(-1-lambda) *
                        (u^(-lambda) + v^(-lambda) - 1)^(-2 - 1/lambda))
      
    } else if (copula == "Frank") {
      th <- lambda
      # Note: The original code's expression for Frank's derivative is complex.
      # Using the standard expression: c(u,v;theta) = -th*exp(-th*(u+v))*(1-exp(-th))/((1-exp(-th*u)-exp(-th*v)+exp(-th*(u+v)))^2)
      a <- th*exp(-th*(u + v))*(-1+exp(-th))
      b <- (1 - exp(-th*u) - exp(-th*v) + exp(-th*(u+v)))^2
      dens_cop <- log(a/b)
      
    } else if (copula == "Gaussian") {
      z1 <- qnorm(u)
      z2 <- qnorm(v)
      rho <- lambda
      # Log-density of Gaussian Copula
      dens_cop <- ((-0.5*log(1-rho^2) -
                      (rho^2*(z1^2+z2^2) - 2*rho*z1*z2) / (2*(1-rho^2))))
      
    } else if (copula == "Gumbel") {
      th <- lambda
      lu <- -log(u)
      lv <- -log(v)
      A  <- (lu^th + lv^th)^(1/th)
      
      logC <- -A
      log_term1 <- (th - 1)*(log(lu) + log(lv))
      log_term2 <- (2/th - 2)*log(lu^th + lv^th)
      log_term3 <- -log(u) - log(v)
      log_term4 <- log1p((th - 1)*lu^th*lv^th/(lu^th + lv^th)) # log(1 + (term))
      
      dens_cop <- logC + log_term1 + log_term2 + log_term3 + log_term4
    }
    
    # Log-likelihood = log(f1) + log(f2) + log(c(u,v))
    ll <- d1 + d2 + dens_cop
    # Negative sum of log-likelihood (for minimization)
    return(-sum(ll))
  })
}
#
#--------------------------------------------------------------------
#
# Maximum Likelihood Estimation for BSBC Models
#
# Estimates the parameters of the Bivariate Simplex/Beta Distribution Copula
# (BSBC) model using Maximum Likelihood via the \code{nlminb} optimizer.
#
estimBSB <- function(y, dist = c("Simplex","Beta"),
                     copula = c("FGM","Clayton","Frank","Gaussian", "Gumbel","all")) {
  
  # Friendly distribution validation
  valid_dists <- c("Simplex", "Beta")
  if (!dist %in% valid_dists) {
    message("Invalid distribution: '", dist, "'")
    message("Valid distributions: ", paste(valid_dists, collapse = ", "))
    return(invisible(NULL))
  }
  
  # Friendly copula validation
  valid_copulas <- c("FGM","Clayton","Frank","Gaussian","Gumbel","all")
  if (!copula %in% valid_copulas) {
    message("Invalid copula: '", copula, "'")
    message("Valid copulas: ", paste(valid_copulas, collapse = ", "))
    return(invisible(NULL))
  }
  
  # Check data structure
  if (!is.data.frame(y) && !is.matrix(y)) {
    message("Data 'y' must be a data.frame or matrix with two columns")
    return(invisible(NULL))
  }
  
  if (ncol(y) != 2) {
    message("Data 'y' must have exactly two columns")
    return(invisible(NULL))
  }
  
  # Load numDeriv for hessian calculation
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    message("Package 'numDeriv' is required for standard error calculation")
    return(invisible(NULL))
  }
  
  N <- nrow(y)
  
  # --- Initial Values (based on transformed means and approximate dispersion) ---
  init <- c(
    qlogis(pmin(pmax(mean(y[,1], na.rm = TRUE),0.05),0.95)), # mu1
    qlogis(pmin(pmax(mean(y[,2], na.rm = TRUE),0.05),0.95)), # mu2
    log(sd(y[,1], na.rm = TRUE)/mean(y[,1], na.rm = TRUE) + 0.1), # sigma1 (log-transformed)
    log(sd(y[,2], na.rm = TRUE)/mean(y[,2], na.rm = TRUE) + 0.1), # sigma2 (log-transformed)
    0.05                                      # lambda (copula parameter)
  )
  
  fit_fun <- function(cop) {
    cat("\n--- Copula", cop, "---\n")
    
    # Check for missing values
    if (any(is.na(y))) {
      message("Data contains missing values. These observations will be removed.")
      y_complete <- y[complete.cases(y), ]
    } else {
      y_complete <- y
    }
    
    fit <- tryCatch({
      nlminb(
        start = init,
        objective = densitBSB,
        y = y_complete,
        dist = dist,
        copula = cop,
        control = list(eval.max = 100000, iter.max = 100000)
      )
    }, error = function(e) {
      message("Optimization failed for copula ", cop, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(fit)) {
      return(NULL)
    }
    
    # --- Parameter Transformation (back to natural scale) ---
    est <- c(
      mu1    = plogis(fit$par[1]),
      mu2    = plogis(fit$par[2]),
      sigma1 = exp(fit$par[3]),
      sigma2 = exp(fit$par[4]),
      lambda = if (cop %in% c("FGM","Gaussian")) tanh(fit$par[5]) else fit$par[5]
    )
    
    # --- Hessian and Standard Error Calculation ---
    H <- try(numDeriv::hessian(densitBSB, fit$par,
                               y = y_complete, dist = dist, copula = cop), silent = TRUE)
    
    if (inherits(H,"try-error") || any(is.na(H)) || det(H) <= 0) {
      se <- rep(NA, length(est))
    } else {
      # Jacobian matrix (for delta method)
      J <- diag(c(
        est["mu1"]*(1-est["mu1"]), # d(mu)/d(eta)
        est["mu2"]*(1-est["mu2"]), # d(mu)/d(eta)
        est["sigma1"],             # d(sigma)/d(log(sigma))
        est["sigma2"],             # d(sigma)/d(log(sigma))
        if (cop %in% c("FGM","Gaussian")) (1 - est["lambda"]^2) # d(lambda)/d(atanh(lambda))
        # Note: Jacobian for Gumbel/Clayton/Frank depends on the true link,
        # but here we use the original code's simplification for other copulas.
        else 1
      ))
      
      # Covariance matrix on natural scale: J * H^(-1) * J^T
      cov_transf <- J %*% solve(H) %*% t(J)
      se <- sqrt(diag(cov_transf))
    }
    
    # --- Results Table ---
    results <- data.frame(
      Estimate = round(est,4),
      Std_Error = round(se,4),
      CI95_Lower  = round(est - 1.96 * se, 4),
      CI95_Upper  = round(est + 1.96 * se, 4)
    )
    
    cat(sprintf("Estimates for %s (Marginal: %s):\n", cop, dist))
    print(results)
    
    # --- Correlation Measures ---
    tau_val <- cor(y_complete[,1], y_complete[,2], method="kendall")
    rho_val <- cor(y_complete[,1], y_complete[,2], method="pearson")
    
    cat("\n")
    cat("Kendall tau:", round(tau_val,3), "\n")
    cat("Pearson rho:", round(rho_val,3), "\n\n")
    
    # --- Information Criteria ---
    logLik_val <- -fit$objective
    k <- length(fit$par)
    n <- nrow(y_complete)
    AIC_val <- 2*k - 2*logLik_val
    BIC_val <- log(n)*k - 2*logLik_val
    
    # Remove the fitted parameters (encoded) for cleaner output structure
    fit_mod <- within(fit, rm(par))
    class(fit_mod) <- "BSBC" # Mark as a BSBC fit object
    
    list(copula = cop,
         dist = dist,
         parameters = est,
         fit = fit_mod,
         table = results,
         AIC = AIC_val,
         BIC = BIC_val)
  }
  
  # --- Fit all copulas or just the specified one ---
  if (copula == "all") {
    copulas_to_fit <- c("FGM","Clayton","Frank","Gaussian","Gumbel")
    fits <- list()
    
    for (cop in copulas_to_fit) {
      fit_result <- fit_fun(cop)
      if (!is.null(fit_result)) {
        fits <- c(fits, list(fit_result))
      }
    }
    
    if (length(fits) == 0) {
      message("No copulas could be successfully fitted")
      return(invisible(NULL))
    }
    
    AICs <- sapply(fits, function(x) x$AIC)
    names(AICs) <- sapply(fits, function(x) x$copula)
    
    BICs <- sapply(fits, function(x) x$BIC)
    names(BICs) <- sapply(fits, function(x) x$copula)
    
    best_AIC <- names(which.min(AICs))
    best_BIC <- names(which.min(BICs))
    
    cat("\n=== Final Comparison ===\n")
    cat("Best Copula (min AIC):", best_AIC, "(", round(AICs[best_AIC], 4), ")\n")
    cat("Best Copula (min BIC):", best_BIC, "(", round(BICs[best_BIC], 4), ")\n")
    class(fits) <- "BSBC_list" # Mark as a list of BSBC fits
    return(fits)
  } else {
    fit <- fit_fun(copula)
    if (is.null(fit)) {
      return(invisible(NULL))
    }
    class(fit) <- "BSBC" # Mark as a single BSBC fit
    return(fit)
  }
}
#
#--------------------------------------------------------------------
#
# Summary Method for BSBC Model Fits
#
# Provides a summary of a single fitted BSBC model or a comparison table for
# multiple fitted BSBC models.
#
summaryBSBC <- function(fit_list, digits = 4) {
  
  if (inherits(fit_list, "BSBC")) {
    fit <- fit_list
    cat("=== BIVARIATE MODEL SUMMARY ===\n")
    cat("Marginal Distribution:", fit$dist, "\n")
    cat("Copula:", fit$copula, "\n")
    cat("Convergence:", ifelse(fit$fit$convergence==0,"Successful","Issues (Check Message)"), "\n")
    if(fit$fit$convergence!=0) cat("Message:", fit$fit$message,"\n")
    cat("Log-Likelihood:", round(-fit$fit$objective,digits), "\n")
    cat("AIC:", round(fit$AIC,digits), "\n")
    cat("BIC:", round(fit$BIC,digits), "\n\n")
    cat("PARAMETER ESTIMATES:\n")
    print(fit$table)
    
    # Note: Requires 'y' data to be available in the environment for correlation calculation
    # Since 'y' is not passed into the summary function, we cannot accurately calculate
    # the sample correlations here. It is better to rely on the correlations printed
    # in the estimBSB function or require 'y' as an argument.
    # For now, we will assume the data 'y' is available in the global environment
    # or rely on the user to check the output of estimBSB.
    # cat("\nSample Kendall's tau:", round(cor(y[,1], y[,2], method="kendall"), 3), "\n")
    # cat("Sample Pearson's rho:", round(cor(y[,1], y[,2], method="pearson"), 3), "\n")
    
  } else if (inherits(fit_list, "BSBC_list")) {
    cat("=== MODEL COMPARISON ===\n")
    comp_df <- do.call(rbind, lapply(fit_list, function(x) {
      data.frame(Copula     = x$copula,
                 Dist       = x$dist,
                 Converged  = (x$fit$convergence==0),
                 LogLik     = round(-x$fit$objective,digits),
                 AIC        = round(x$AIC,digits),
                 BIC        = round(x$BIC,digits),
                 stringsAsFactors = FALSE)
    }))
    print(comp_df[order(comp_df$AIC),]) # Print sorted by AIC
    invisible(comp_df) # Invisibly return the comparison data frame
  } else {
    stop("Invalid input format for summary.BSBC. Must be a single BSBC fit object or a list of BSBC fits.")
  }
}
#
#--------------------------------------------------------------------
#
# Simulate Data from BSBC Models
#
# Generates random samples from a Bivariate Simplex or Beta distribution
# with a specified copula structure.
#
rBSBC <- function(n, copula = c("FGM","Clayton","Frank","Gumbel","Gaussian"),
                  theta = 0.5,
                  dist = c("Simplex","Beta"),
                  par1 = c(0.5, 0.3),
                  par2 = c(0.5, 0.3)) {
  
  # Friendly copula validation
  valid_copulas <- c("FGM","Clayton","Frank","Gumbel","Gaussian")
  if (!copula %in% valid_copulas) {
    message("Invalid copula: '", copula, "'")
    message("Valid copulas: ", paste(valid_copulas, collapse = ", "))
    return(invisible(NULL))
  }
  
  # Friendly distribution validation
  valid_dists <- c("Simplex", "Beta")
  if (!dist %in% valid_dists) {
    message("Invalid distribution: '", dist, "'")
    message("Valid distributions: ", paste(valid_dists, collapse = ", "))
    return(invisible(NULL))
  }
  
  # --- Theta range checks ---
  is_valid <- TRUE; msg <- ""
  if (copula == "FGM" && (theta < -1 || theta > 1)) {
    is_valid <- FALSE
    msg <- paste0("For FGM copula, theta must be in [-1, 1]. Value provided: ", theta)
  } else if (copula == "Clayton" && theta <= 0) {
    is_valid <- FALSE
    msg <- paste0("For Clayton copula, theta must be > 0. Value provided: ", theta)
  } else if (copula == "Frank" && theta == 0) {
    is_valid <- FALSE
    msg <- paste0("For Frank copula, theta cannot be 0. Value provided: ", theta)
  } else if (copula == "Gumbel" && theta < 1) {
    is_valid <- FALSE
    msg <- paste0("For Gumbel copula, theta must be ≥ 1. Value provided: ", theta)
  } else if (copula == "Gaussian" && (theta <= -1 || theta >= 1)) {
    is_valid <- FALSE
    msg <- paste0("For Gaussian copula, theta must be in (-1, 1). Value provided: ", theta)
  }
  
  if (!is_valid) {
    message("Invalid theta value: ", msg)
    return(invisible(NULL))
  }
  
  # --- Package check ---
  if (!requireNamespace("copula", quietly = TRUE)) {
    message("Package 'copula' is required for simulation")
    return(invisible(NULL))
  }
  if (!requireNamespace("gamlss.dist", quietly = TRUE)) {
    message("Package 'gamlss.dist' is required for Simplex/Beta quantiles")
    return(invisible(NULL))
  }
  
  # --- Define Copula Object ---
  cop_obj <- switch(copula,
                    "FGM"     = copula::fgmCopula(theta, dim = 2),
                    "Clayton" = copula::claytonCopula(theta, dim = 2),
                    "Frank"   = copula::frankCopula(theta, dim = 2),
                    "Gumbel"  = copula::gumbelCopula(theta, dim = 2),
                    "Gaussian"= copula::normalCopula(theta, dim = 2)
  )
  
  # --- Generate Uniform Samples from Copula ---
  U <- copula::rCopula(n, copula = cop_obj)
  
  # --- Transform by Marginal Quantile Functions ---
  if (dist == "Beta") {
    # Check if Beta functions are available
    if (!exists("qBE")) {
      message("Beta quantile function (qBE) not found. Please load gamlss.dist package")
      return(invisible(NULL))
    }
    Y1 <- gamlss.dist::qBE(U[,1], mu = par1[1], sigma = par1[2])
    Y2 <- gamlss.dist::qBE(U[,2], mu = par2[1], sigma = par2[2])
  } else if (dist == "Simplex") {
    # Check if Simplex functions are available
    if (!exists("qSIMPLEX")) {
      message("Simplex quantile function (qSIMPLEX) not found. Please load gamlss.dist package")
      return(invisible(NULL))
    }
    Y1 <- gamlss.dist::qSIMPLEX(U[,1], mu = par1[1], sigma = par1[2])
    Y2 <- gamlss.dist::qSIMPLEX(U[,2], mu = par2[1], sigma = par2[2])
  }
  
  return(data.frame(y1 = Y1, y2 = Y2))
}
#
#------------------------------------------------------------------
#
# Bivariate Simplex Joint Expectation (FGM Copula)
#
# Calculates the joint expectation E[X*Y] for the Bivariate Simplex
# distribution with the FGM copula, using either Bessel function integration
# or Struve function approximations.
#
eBsimplex <- function(theta,
                      method = c("struve","bessel"),
                      tol_struve = 1e-12,
                      maxk = 2000,
                      int_rel_tol = 1e-8,
                      int_abs_tol = 0,
                      int_subdivisions = 500) {
  
  method <- match.arg(method)
  
  mu_x <- theta[1]; mu_y <- theta[2]
  sigma_x <- theta[3]; sigma_y <- theta[4]
  lambda <- theta[5]
  
  # --- Parameter checks ---
  if (!(0 <= mu_x && mu_x <= 1)) stop("mu_x out of range [0,1]")
  if (!(0 <= mu_y && mu_y <= 1)) stop("mu_y out of range [0,1]")
  if (!(sigma_x > 0)) stop("sigma_x must be > 0")
  if (!(sigma_y > 0)) stop("sigma_y must be > 0")
  if (!(-1 <= lambda && lambda <= 1)) stop("lambda out of range [-1,1]")
  
  # --- Pre-calculations ---
  # ... (ex, ey, cx, cy, ax, ay as in your original code) ...
  ex <- 1/mu_x - 1
  ey <- 1/mu_y - 1
  
  cx <- 1 / (sigma_x * sqrt(2*pi))
  cy <- 1 / (sigma_y * sqrt(2*pi))
  
  ax <- (ex + 1)^2 / (sigma_x^2 * ex)
  ay <- (ey + 1)^2 / (sigma_y^2 * ey)
  
  K <- function(x, v) base::besselK(x, v)
  
  # --- Bessel Method ---
  if (method == "bessel") {
    # ... (Bessel integration logic from your original code) ...
    # Assumes safe_integrate is defined locally or available (using stats::integrate)
    safe_integrate <- function(f) {
      val <- tryCatch(
        integrate(f, lower = 1, upper = Inf,
                  rel.tol = int_rel_tol, abs.tol = int_abs_tol,
                  subdivisions = int_subdivisions)$value,
        error = function(e) NA_real_
      )
      return(val)
    }
    
    J0x <- safe_integrate(function(x) K(ax * (x^2 + 1) / x, -1))
    J1x <- safe_integrate(function(x) (1/x) * K(ax * (x^2 + 1) / x, 0))
    J0y <- safe_integrate(function(y) K(ay * (y^2 + 1) / y, -1))
    J1y <- safe_integrate(function(y) (1/y) * K(ay * (y^2 + 1) / y, 0))
    
    if (any(is.na(c(J0x, J1x, J0y, J1y)))) {
      warning("Bessel integrals failed; result returned as NA")
      return(NA_real_)
    }
    
    E <- mu_x*mu_y + lambda * (2*cx^2*exp(2*ax)*(1/ex * J1x + J0x) - mu_x) *
      (2*cy^2*exp(2*ay)*(1/ey * J1y + J0y) - mu_y)
    
    # Numerical safeguard
    E <- min(max(E, 0), 1)
    
    return(E)
  }
  
  # --- Struve Method (Approximation) ---
  if (method == "struve") {
    # The modified_struve helper function must be defined within the package scope.
    # ... (Struve approximation logic from your original code, including RandomFieldsUtils check) ...
    
    # Helper: modified_struve via series (for internal use if RFUtils not available)
    modified_struve <- function(z, nu, tol = 1e-12, max.iter = 2000) {
      if (length(z) > 1) {
        return(sapply(z, modified_struve, nu = nu, tol = tol, max.iter = max.iter))
      }
      if (z == 0) return(0)
      
      k <- 0L
      term <- ((-1)^0 * (z/2)^(2*0 + nu + 1)) / (gamma(1.5) * gamma(nu + 1.5))
      s <- term
      while (abs(term) > tol && k < max.iter) {
        k <- k + 1L
        term <- ((-1)^k * (z/2)^(2*k + nu + 1)) / (gamma(k + 1.5) * gamma(k + nu + 1.5))
        s <- s + term
      }
      return(s)
    }
    
    
    if (requireNamespace("RandomFieldsUtils", quietly = TRUE)) {
      L <- function(x, v) RandomFieldsUtils::struveH(x, v)
    } else {
      L <- function(x, v) modified_struve(x, v, tol = tol_struve, max.iter = maxk)
    }
    
    # Calculate Struve and Bessel values at 2*ax and 2*ay
    K0_2ax <- K(2*ax, 0); K1_2ax <- K(2*ax, 1)
    K0_2ay <- K(2*ay, 0); K1_2ay <- K(2*ay, 1)
    Lm1_2ax <- L(2*ax, -1); L0_2ax <- L(2*ax, 0)
    Lm1_2ay <- L(2*ay, -1); L0_2ay <- L(2*ay, 0)
    
    # Revised E expression
    Ax <- (K0_2ax * Lm1_2ax + K1_2ax * L0_2ax)
    Ay <- (K0_2ay * Lm1_2ay + K1_2ay * L0_2ay)
    
    E <- mu_x*mu_y + lambda * ((cx^2 * pi * Ax) - mu_x) * ((cy^2 * pi * Ay) - mu_y)
    
    # Numerical safeguard
    E <- min(max(E, 0), 1)
    
    return(E)
  }
}
#
#---------------------------------------------------------------
#
# Helper: Clamps values between (0,1)
clamp01 <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)

# --- Helper function for joint density (rescaled for plotting) ---
# Note: This function computes the joint density f(y1, y2) = f1(y1) * f2(y2) * c(u, v)
# It is rescaled (log(f) - max(log(f))) for numerical stability and visualization.
fxy <- function(y1, y2, param,
                dist = c("Simplex","Beta"),
                copula_type = c("FGM","Clayton","Frank","Gumbel","Gaussian")) {
  
  dist <- match.arg(dist)
  copula_type <- match.arg(copula_type)
  
  mu1 <- param[1]; mu2 <- param[2]
  sigma1 <- param[3]; sigma2 <- param[4]
  lambda <- param[5]
  
  if (length(y1) != length(y2)) stop("y1 and y2 must have the same length (use expand.grid).")
  
  # --- Marginal CDFs and Densities ---
  # Requires gamlss.dist::pSIMPLEX, dSIMPLEX, pBE, dBE to be available
  if (dist == "Simplex") {
    u <- gamlss.dist::pSIMPLEX(y1, mu1, sigma1)
    v <- gamlss.dist::pSIMPLEX(y2, mu2, sigma2)
    f1 <- gamlss.dist::dSIMPLEX(y1, mu1, sigma1)
    f2 <- gamlss.dist::dSIMPLEX(y2, mu2, sigma2)
  } else { # Beta
    u <- gamlss.dist::pBE(y1, mu1, sigma1)
    v <- gamlss.dist::pBE(y2, mu2, sigma2)
    f1 <- gamlss.dist::dBE(y1, mu1, sigma1)
    f2 <- gamlss.dist::dBE(y2, mu2, sigma2)
  }
  
  # --- Numerical protection ---
  u <- clamp01(u)
  v <- clamp01(v)
  logf1 <- log(pmax(f1, 1e-300))
  logf2 <- log(pmax(f2, 1e-300))
  
  # --- Copula Density (Log-scale) ---
  n <- length(u)
  dens_cop_log <- rep(NA_real_, n)
  eps <- 1e-12
  
  if (copula_type == "FGM") {
    # Lambda for Gaussian/FGM is restricted to (-1, 1) in estimBSB, but param[5] is raw.
    # We must ensure lambda is within valid bounds for plotting, here we rely on param[5] being the estimated value.
    if (copula_type == "Gaussian" || copula_type == "FGM") {
      lambda <- pmin(pmax(lambda, -1 + eps), 1 - eps)
    }
    
    dens_cop_log <- log1p(lambda * (1 - 2*u) * (1 - 2*v))
    
  } else if (copula_type == "Clayton") {
    val <- u^(-lambda) + v^(-lambda) - 1
    val <- pmax(val, eps)
    dens_cop_log <- log1p(lambda) + (-1 - lambda)*log(u) + (-1 - lambda)*log(v) +
      (-2 - 1/lambda) * log(val)
    
  } else if (copula_type == "Frank") {
    th <- lambda
    if (abs(th) < 1e-10) {
      dens_cop_log <- rep(0, n)
    } else {
      num <- th * exp(-th*(u + v)) * (1 - exp(-th))
      den <- (1 - exp(-th*u) - exp(-th*v) + exp(-th*(u+v)))^2
      den <- pmax(den, eps)
      dens_cop_log <- log(pmax(abs(num), eps)) - log(den)
    }
  } else if (copula_type == "Gumbel") {
    th <- lambda
    lu <- pmax(-log(u), eps)
    lv <- pmax(-log(v), eps)
    A <- (lu^th + lv^th)^(1/th)
    logC <- -A
    log_term1 <- (th - 1) * (log(lu) + log(lv))
    log_term2 <- (2/th - 2) * log(lu^th + lv^th)
    log_term3 <- -log(u) - log(v)
    frac <- (th - 1) * lu^th * lv^th / pmax(lu^th + lv^th, eps)
    frac <- pmax(frac, -0.999999)
    log_term4 <- log1p(frac)
    dens_cop_log <- logC + log_term1 + log_term2 + log_term3 + log_term4
  } else if (copula_type == "Gaussian") {
    # Use the same logic as in densitBSB for consistency
    rho <- lambda
    z1 <- qnorm(u)
    z2 <- qnorm(v)
    dens_cop_log <- ((-0.5*log(1-rho^2) -
                        (rho^2*(z1^2+z2^2) - 2*rho*z1*z2) / (2*(1-rho^2))))
  }
  
  
  # --- Joint Density ---
  log_z <- logf1 + logf2 + dens_cop_log
  
  # Rescale for visualization (prevents overflow and underflow)
  maxlog <- max(log_z, na.rm = TRUE)
  z_scaled <- exp(log_z - maxlog)
  attr(z_scaled, "scale_factor") <- exp(maxlog)
  z_scaled[is.na(z_scaled) | is.infinite(z_scaled)] <- min(z_scaled[!is.na(z_scaled) & !is.infinite(z_scaled)], na.rm = TRUE)
  return(z_scaled)
}
#-------------------------
# Plot BSD Joint Density
plotBSD <- function(y1, y2, param,
                    dist = c("Simplex","Beta"),
                    copula_type = c("FGM","Clayton","Frank","Gumbel","Gaussian"),
                    plot_type = c("surface3d","contour"),
                    title = NULL,
                    n_grid = 500) {
  
  plot_type <- match.arg(plot_type)
  dist <- match.arg(dist)
  copula_type <- match.arg(copula_type)
  
  # --- Grid creation ---
  y1_seq <- seq(min(y1), max(y1), length.out = n_grid)
  y2_seq <- seq(min(y2), max(y2), length.out = n_grid)
  
  # --- Grid completo ---
  grid <- expand.grid(y1 = y1_seq, y2 = y2_seq)
  grid$z <- fxy(grid$y1, grid$y2, param, dist, copula_type) # CHAMA fxy AQUI
  
  # --- Limpeza de NA/Inf ---
  min_z <- min(grid$z[!is.na(grid$z) & !is.infinite(grid$z)], na.rm = TRUE)
  if (is.infinite(min_z) | length(min_z) == 0) min_z <- 0
  grid$z[is.na(grid$z) | is.infinite(grid$z)] <- min_z
  
  # --- Matriz zmat ---
  zmat <- matrix(grid$z, nrow = n_grid, ncol = n_grid, byrow = TRUE)
  
  # --- Paleta azul → amarelo ---
  col_fun <- grDevices::colorRampPalette(c("lightblue", "lemonchiffon"))
  cols <- col_fun(100)
  
  z_range <- range(zmat, na.rm = TRUE)
  map_to_col <- function(z) {
    idx <- round((z - z_range[1]) / diff(z_range) * (length(cols)-1)) + 1
    cols[pmin(pmax(idx,1), length(cols))]
  }
  
  if (plot_type == "surface3d") {
    # --- Superfície 3D ---
    graphics::persp(x = y1_seq, y = y2_seq, z = zmat,
                    theta = 50, phi = 20, expand = 0.6,
                    col = "lemonchiffon",
                    border = "lightblue",
                    ticktype = "detailed",
                    xlab = "y1", ylab = "y2", zlab = "dCopula",
                    main = ifelse(is.null(title), "3D Surface", title),
                    box = TRUE,
                    cex.main = 1,
                    cex.lab = 1,
                    cex.axis = 0.7)
  } else if (plot_type == "contour") {
    # --- Contorno 2D ---
    graphics::filled.contour(x = y1_seq, y = y2_seq, z = zmat,
                             color.palette = function(n) col_fun(n),
                             zlim = z_range,
                             plot.title = graphics::title(main = ifelse(is.null(title), "Contour Plot", title),
                                                          xlab = "y1", ylab = "y2"),
                             plot.axes = {
                               graphics::axis(1); graphics::axis(2)
                               graphics::contour(x = y1_seq, y = y2_seq, z = zmat, add = TRUE,
                                                 lwd = 1.2, col = "gray40")
                             },
                             key.title = graphics::title("dCopula"))
  }
}
#----------------------------------------------------------------