# load necessary libraries
library(MASS)
library(glmnet)
library(ncvreg)
library(SIS) 

# set seed for reproducibility
set.seed(123456789)

# simulation parameters 
n <- 600
p <- 2000
num_simulations <- 200
s_vals <- c(3, 6, 12, 24)

# calculate nsis based on the paper (n/log(n))
nsis_val <- floor(n / log(n)) 
cat("Using nsis =", nsis_val, "for SIS methods.\n")

# create a data frame to store results
results <- data.frame(
  s = s_vals,
  M_lambda_max = NA, RSD_lambda_max = NA,  
  M_MMLE_rank = NA, RSD_MMLE_rank = NA,
  SIS_MLR = NA, RSD_MLR = NA,
  SIS_MMLE = NA, RSD_MMLE = NA,
  LASSO = NA, RSD_LASSO = NA,
  SCAD = NA, RSD_SCAD = NA
)

# data generation function for setting 3 (S3)
generate_data <- function(n, p, s) {
  # generate underlying standard normal variables for all p potential positions first
  Z <- matrix(rnorm(n * p), nrow = n)

  # create the final X matrix
  X <- matrix(0, nrow = n, ncol = p)

  # first p-50 variables are just independent standard normals
  X[, 1:(p-50)] <- Z[, 1:(p-50)]

  # generate the last 50 variables based on the first s variables from Z
  for (k in 1:50) {
    col_idx <- p - 50 + k

    # sum over the first s columns of Z
    sum_part <- 0
    if (s > 0) {
      # avoid loop if s=0
      for (j in 1:s) {
        # if s is large, check if Z has enough columns
        if(j <= ncol(Z)) {
           sum_part <- sum_part + Z[,j] * ((-1)^(j+1)) / 5
        }
      }
    }

    # add the error term using independent N(0,1) noise
    error_term <- sqrt(max(0, 25 - s)) / 5 * rnorm(n) # max(0,...) handles potential s>25 case
    X[, col_idx] <- sum_part + error_term
  }

  # standardize all variables after generation
  X <- scale(X)
  # handle cases where a column might have 0 variance after construction
  X[, apply(X, 2, var) == 0] <- 0

  # define sets of active coefficient vector based on s
  if (s == 3) {
    beta_active <- c(3, 4, 3)
  } else if (s == 6) {
    beta_active <- c(3, -3, 3, -3, 3, -3)
  } else if (s == 12) {
    beta_active <- rep(c(3, 4), 6)
  } else if (s == 24) {
    beta_active <- rep(c(3, 4), 12)
  } else {
    # use this beta coefficients pattern in case of error
    beta_active <- rep(3, s) 
  }

  # ensure beta_active has length s
  if (length(beta_active) != s) {
      stop("Length of beta_active does not match s")
  }

  # generate beta (true coefficients are in positions p-50+1 to p-50+s)
  beta <- c(rep(0, p-50), beta_active, rep(0, 50-s))
  true_vars <- which(beta != 0)

  # generate response using logistic regression
  eta <- X %*% beta
  # prevent extreme probabilities due to large eta values
  eta <- pmax(-15, pmin(15, eta))
  prob <- exp(eta) / (1 + exp(eta))
  Y <- rbinom(n, 1, prob)

  # return list including X, Y, and true_vars
  list(X = X, Y = Y, true_vars = true_vars)
}


# evaluate SIS method using MLR screening
evaluate_sis_mlr <- function(X, Y, true_vars, nsis) {
  p <- ncol(X)
  n <- nrow(X)
  
  if (nsis <= 0 || nsis > p) {
      stop("nsis must be between 1 and p")
  }
  
  # initialize list of selected variables
  selected <- NULL
  
  # marginal Likelihood Ratio screening
  lr_stats <- numeric(p)
  null_model <- tryCatch({
      glm(Y ~ 1, family = binomial)
  }, error = function(e) NULL)

  if (is.null(null_model)) {
      warning("Null model failed to converge. Returning p.")
      # return if null model fails
      return(p) 
  }
  null_loglik <- tryCatch(logLik(null_model), error = function(e) -Inf)

  # calculate likelihood ratio for each variable
  for (j in 1:p) {
    full_model <- tryCatch({
        glm(Y ~ X[,j], family = binomial)
    }, warning = function(w){ NULL }, error = function(e){ NULL })
    
    if (!is.null(full_model)) {
       full_loglik <- tryCatch(logLik(full_model), error = function(e) -Inf)
       # check if logLik values are valid before calculating LR stat
       if (is.finite(full_loglik) && is.finite(null_loglik)) {
          lr_stats[j] <- max(0, 2 * (full_loglik - null_loglik)) # ensure non-negative
       } else {
          # if logLik calculation failed, assign 0
          lr_stats[j] <- 0
       }
    } else {
      # assign 0 if model failed
       lr_stats[j] <- 0
    }
  }
  
  # select variables based on likelihood ratio ranking
  ranked_indices <- order(lr_stats, decreasing = TRUE)
  selected <- ranked_indices[1:nsis]
    
  # find the minimum model size (k) within the top nsis variables
  for (k in 1:nsis) {
    # ensure selected indices are valid before checking
    current_selection <- selected[1:k]
    if (length(current_selection) > 0 && !any(is.na(current_selection))) {
        # check if all true variables are in the current top k selection
        if (all(true_vars %in% current_selection)) {
          return(k)
        }
    } else {
      # if selection itself is invalid with NA indices, treat as failure and return p
      warning("NA indices found in MLR selection.")
      return(p) 
    }
  }
  
  # if the loop completes without finding all true variables within top nsis, return p
  return(p)
}

# evaluate SIS method using mmle screening 
evaluate_sis_mmle <- function(X, Y, true_vars, nsis){
  p <- ncol(X)
  n <- nrow(X)

  if (nsis <= 0 || nsis > p) {
      stop("nsis must be between 1 and p")
  }

  # use SIS package to get the indices of the selected variables
  sis_result <- tryCatch({
    SIS::SIS(X, Y, family = "binomial", nsis = nsis, iter = FALSE, standardize = FALSE)
  }, error = function(e){
    warning(paste("SIS package function failed:", e$message))
    return(NULL)
  })

  if(is.null(sis_result) || is.null(sis_result$sis.ix0)){
      warning("SIS package did not return valid indices (sis.ix0). Returning p.")
      return(p)
  }
  
  # get the indices selected by the initial screening
  selected_by_pkg <- sis_result$sis.ix0
  
  # check if the length matches nsis (it should)
  if(length(selected_by_pkg) != nsis){
      warning(paste("SIS package returned", length(selected_by_pkg), "indices, expected", nsis))
      if(length(selected_by_pkg) < nsis) {
          selected_by_pkg <- c(selected_by_pkg, rep(NA, nsis - length(selected_by_pkg)))
      } else {
          selected_by_pkg <- selected_by_pkg[1:nsis]
      }
  }

  # find the minimum model size (k) within the top nsis variables
  for (k in 1:nsis) {
    current_selection <- selected_by_pkg[1:k]
    # ensure indices are valid before checking subset
    if (length(current_selection) > 0 && !any(is.na(current_selection))) {
        # check if all true variables are in the current top k selection
        if (all(true_vars %in% current_selection)) {
          return(k)
        }
    } else {
       # if selection itself is invalid (contains NA), treat as failure and return p
       warning("NA indices found in SIS package selection.")
       return(p)
    }
  }

  # if the loop completes without finding all true variables within top nsis, return p
  return(p)
}


# evaluate ranking by mmle
evaluate_mmle_rank <- function(X, Y, true_vars) {
  p <- ncol(X)
  n <- nrow(X)
  
  # calculate marginal regression coefficients
  beta_marginal <- numeric(p)
  for (j in 1:p) {
    marginal_model <- tryCatch({
        glm(Y ~ X[,j], family = binomial)
    }, warning = function(w){ NULL }, error = function(e){ NULL })
    
    # check if model converged and has coefficients
    if (!is.null(marginal_model) && length(coef(marginal_model)) > 1) {
      beta_marginal[j] <- coef(marginal_model)[2] 
    } else {
      # if model failed or is intercept-only, assign 0
      beta_marginal[j] <- 0
    }
  }
  
  # rank all variables based on magnitude of marginal coefficients
  beta_marginal[is.na(beta_marginal)] <- 0 
  ranked_indices <- order(abs(beta_marginal), decreasing = TRUE)
  
  # find the minimum model size (k) in the full ranking
  for (k in 1:p) {
     # ensure ranked indices are valid
    current_selection <- ranked_indices[1:k]
     if (length(current_selection) > 0 && !any(is.na(current_selection))) {
        if (all(true_vars %in% current_selection)) {
           return(k)
        }
     } else {
        # if ranking is invalid, treat as failure and return p
        warning("NA indices found in MMLE full ranking.")
        return(p)
     }
  }
  return(p) 
}


# evaluate LASSO with 10-fold cross-validation
evaluate_glmnet <- function(X, Y, true_vars) {
  p <- ncol(X)
  # use 10-fold cross-validation to select lambda
  cv_fit <- tryCatch({
    cv.glmnet(X, Y, family = "binomial", alpha = 1, 
                      standardize = FALSE, nfolds = 10)
  }, error = function(e) NULL)
  
  if (is.null(cv_fit)) {
    warning("cv.glmnet failed. Returning p.")
    return(p)
  }

  # get coefficients at lambda.1se 
  coeffs <- tryCatch({
      coef(cv_fit, s = "lambda.1se")
  }, error = function(e) NULL)

  if (is.null(coeffs)) {
      warning("Getting coefficients from cv.glmnet failed. Returning p.")
      return(p)
  }
  
  selected_coeffs <- which(coeffs[-1] != 0)
  selected_indices <- as.numeric(selected_coeffs) 
  
  # check if all true variables are included
  if (length(selected_indices) > 0 && all(true_vars %in% selected_indices)) {
    return(length(selected_indices))
  }
  # if not all true variables are selected, or if selection failed, return p
  return(p)  
}

# evaluate SCAD with 10-fold cross-validation
evaluate_ncvreg <- function(X, Y, true_vars) {
  p <- ncol(X)
  # use 10-fold cross-validation to select lambda
   cv_fit <- tryCatch({
    cv.ncvreg(X, Y, family = "binomial", penalty = "SCAD", 
                       gamma = 3.7, nfolds = 10)
   }, error = function(e) NULL)

   if (is.null(cv_fit)) {
      warning("cv.ncvreg failed. Returning p.")
      return(p)
   }
  
  # get coefficients at lambda.min (using lambda.min as is common for ncvreg)
   coeffs <- tryCatch({
     coef(cv_fit, lambda = cv_fit$lambda.min)
   }, error = function(e) NULL)

   if (is.null(coeffs)) {
     warning("Getting coefficients from ncvreg failed. Returning p.")
     return(p)
   }
  
  selected_coeffs <- which(coeffs[-1] != 0) # exclude intercept
  # ensure selected indices are numeric
  selected_indices <- as.numeric(selected_coeffs)

  # check if all true variables are included
  if (length(selected_indices) > 0 && all(true_vars %in% selected_indices)) {
    # return the size of the selected set
    return(length(selected_indices))
  }
  # if not all true variables are selected, or if selection failed, return p
  return(p) 
}


# simulation loop
cat("\nStarting simulations (", num_simulations, " runs per 's')...\n")
for (i in seq_along(s_vals)) {
  s <- s_vals[i]
  cat("Running simulations for s =", s, "\n")
  
  # initialize size vectors for this s value
  max_eigenvalues_sim <- numeric(num_simulations)
  sizes_mmle_rank <- numeric(num_simulations)
  sizes_mlr <- numeric(num_simulations)
  sizes_mmle <- numeric(num_simulations) # Will use SIS package results
  sizes_lasso <- numeric(num_simulations)
  sizes_scad <- numeric(num_simulations)
  
  # inner loop for simulations
  for (sim in 1:num_simulations) {
    # progress indicator
    if(sim %% 20 == 0) cat("  Simulation", sim, "/", num_simulations, "...\n") 
    
    # generate data
    dat <- generate_data(n, p, s)
    X <- dat$X # predictor matrix
    Y <- dat$Y # response vector
    true_vars <- dat$true_vars # indices of true variables
    
    # calculate maximum eigenvalue
    max_eigenvalue <- tryCatch({
        # calculate sample covariance matrix
        Sigma_hat <- try(cov(X), silent = TRUE)
        if(inherits(Sigma_hat, "try-error")) {
            stop("Covariance calculation failed.") # Stop if cov fails
        }
        # calculate eigenvalues
        eigen_vals <- eigen(Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
        # return maximum eigenvalues
        max(eigen_vals)
      }, error = function(e) {
          warning(paste("Eigenvalue calculation failed for sim", sim, "s=", s, ":", e$message))
          return(NA)
      })
    max_eigenvalues_sim[sim] <- max_eigenvalue

    # check for issues in generated data (Y)
    if(length(unique(Y)) < 2) {
        warning(paste("Skipping simulation", sim, "for s=", s, "due to constant Y response."))
        # assign p (failure) to all MMS methods for this sim
        sizes_mmle_rank[sim] <- p
        sizes_mlr[sim] <- p
        sizes_mmle[sim] <- p
        sizes_lasso[sim] <- p
        sizes_scad[sim] <- p
        # skip to next simulation
        next 
    }

    # evaluate all MMS methods
    sizes_mmle_rank[sim] <- evaluate_mmle_rank(X, Y, true_vars)
    sizes_mlr[sim] <- evaluate_sis_mlr(X, Y, true_vars, nsis = nsis_val) 
    sizes_mmle[sim] <- evaluate_sis_mmle(X, Y, true_vars, nsis = nsis_val) 
    sizes_lasso[sim] <- evaluate_glmnet(X, Y, true_vars)
    sizes_scad[sim] <- evaluate_ncvreg(X, Y, true_vars)
  }
  
  # store results
  results[i, "M_lambda_max"] <- round(median(max_eigenvalues_sim, na.rm = TRUE), 2)
  results[i, "RSD_lambda_max"] <- round(IQR(max_eigenvalues_sim, na.rm = TRUE) / 1.349, 2)
  results[i, "M_MMLE_rank"] <- round(median(sizes_mmle_rank, na.rm = TRUE), 2)
  results[i, "RSD_MMLE_rank"] <- round(IQR(sizes_mmle_rank, na.rm = TRUE) / 1.349, 2)
  results[i, "SIS_MLR"] <- round(median(sizes_mlr, na.rm = TRUE), 2)
  results[i, "RSD_MLR"] <- round(IQR(sizes_mlr, na.rm = TRUE) / 1.349, 2)
  results[i, "SIS_MMLE"] <- round(median(sizes_mmle, na.rm = TRUE), 2)
  results[i, "RSD_MMLE"] <- round(IQR(sizes_mmle, na.rm = TRUE) / 1.349, 2)
  results[i, "LASSO"] <- round(median(sizes_lasso, na.rm = TRUE), 2)
  results[i, "RSD_LASSO"] <- round(IQR(sizes_lasso, na.rm = TRUE) / 1.349, 2)
  results[i, "SCAD"] <- round(median(sizes_scad, na.rm = TRUE), 2)
  results[i, "RSD_SCAD"] <- round(IQR(sizes_scad, na.rm = TRUE) / 1.349, 2)
}
cat("Simulations complete.\n\n")

# print descriptive headers for the output table
cat("TABLE 5 (Reproduced - Setting S3)\n")
cat("The Median MMS / M-lambda_max and the associated RSD (in parentheses)\n")
cat("Setting 3 (S3), p =", p, ", n =", n, ", num_simulations =", num_simulations, "\n\n")

# print headers
header <- sprintf("%-3s | %-16s | %-17s | %-15s | %-15s | %-15s | %-15s\n",
                  "s", "M-lambda_max(RSD)", "MMLE-Rank(RSD)", "SIS-MLR(RSD)", "SIS-MMLE(RSD)", "LASSO(RSD)", "SCAD(RSD)")
cat(header)
cat(rep("-", nchar(header)), "\n", sep = "")

# print results
for (i in 1:nrow(results)) {
  row <- results[i, ]
  cat(sprintf("%-3d | %-16s | %-17s | %-15s | %-15s | %-15s | %-15s\n",
              row$s,
              paste0(row$M_lambda_max, " (", row$RSD_lambda_max, ")"),
              paste0(row$M_MMLE_rank, " (", row$RSD_MMLE_rank, ")"), 
              paste0(row$SIS_MLR, " (", row$RSD_MLR, ")"),
              paste0(row$SIS_MMLE, " (", row$RSD_MMLE, ")"), # Result from SIS package + MMS calc
              paste0(row$LASSO, " (", row$RSD_LASSO, ")"),
              paste0(row$SCAD, " (", row$RSD_SCAD, ")") 
  )) 
}