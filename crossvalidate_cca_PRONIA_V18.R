library(tidyverse)
library(CCA)
library(caret)
library(dplyr)
library(clue)

cca_dti_rop_chr_cv <- read.csv('~/Documents/Harvard/projects/DTI_Pronia/in_between_data/cca_dti_rop_chr_md5_16.csv')
#cca_vars_rop_chr_cv <- read.csv('~/Documents/Harvard/projects/DTI_Pronia/in_between_data/cca_vars_rop_chr_for_cv.csv')
cca_vars_rop_chr_cv <- read.csv('~/Documents/Harvard/projects/DTI_Pronia/in_between_data/cca_vars_rop_chr_for_cv_NA_V18.csv')

# Function to create permutations of row indices within a data frame
# The function generates a shuffled version of the row numbers (without changing the data itself)
# For cross-validation we permute the data dynamically without altering the original dataset
create_permutation <- function(df, num_perm) {
  # Number of rows in the original data frame
  n <- nrow(df)
  # Pre-allocate a matrix to store the original and shuffled row numbers
  shuffled_matrix <- matrix(nrow = n, ncol = num_perm)
  # Set seed for reproducibility
  set.seed(123)
  # Fill the first column with original row numbers
  shuffled_matrix[, 1] <- seq_len(n)
  # Fill the remaining columns with shuffled row numbers
  for (i in seq_len(num_perm - 1)) {
    shuffled_matrix[, i + 1] <- sample(shuffled_matrix[, 1])
  }
  # Convert to frame for faster manipulation
  result <- as.data.frame(shuffled_matrix)
  # Rename the columns
  names(result) <- c("id", paste0("shuffled_", 1:(num_perm - 1)))
  return(result)
}

# Function to center a data matrix (X) and apply the same centering transformation to a test set (X_test_cv_perm)
# The centering is performed by subtracting the column means of the training data (X) and adjusting for columns that have no variation.
center_cv <- function(X, X_test_cv_perm = NULL, apply_mean = NULL) {
  # Check if apply_mean exists (indicating whether we should apply centering to the test set as well)
  if (!is.null(apply_mean)){
    # Identify constant columns in the training data X (those with no variation after differencing)
    icte <- apply(diff(X, 1, 1)^2, 2, sum) == 0
    # Calculate column means of the training data X and save them
    X_saved_means <- colMeans(X)
    # Center the training data X by subtracting the column means from each column
    X <- sweep(X, 2, colMeans(X), "-")
    # Remove the constant columns from the training data
    X <- X[, !icte]
    # Similarly, identify constant columns in the test data X_test_cv_perm
    icte_test <- apply(diff(X_test_cv_perm, 1, 1)^2, 2, sum) == 0
    # Center the test data by subtracting the saved column means of the X_saved_means derived from training sample from X_test_cv_perm
    X_test_cv_perm <- sweep(X_test_cv_perm, 2, X_saved_means, "-")
    # Remove the constant columns from the test data
    X_test_cv_perm <- X_test_cv_perm[, !icte_test]
    # Return the centered and adjusted test data
    return(X_test_cv_perm)
  }else{
    # If apply_mean does not exist, only center the training data X without adjusting the test set
    # Identify constant columns in the training data X
    icte <- apply(diff(X, 1, 1)^2, 2, sum) == 0
    X <- sweep(X, 2, colMeans(X), "-")
    X_saved_means <- colMeans(X)
    X <- X[, !icte]
    return(X)
  }
}

# Function to perform imputation on training and test datasets within a cross-validation framework
# The function imputes missing values based on group-level statistics (mean for most variables and median for a specified variable)
impute_cv <- function(vars_train_df, vars_test_df, group_var, median_var, which_imputation){
  # Extract the names of the columns in the training data frame for later selection
  vars_to_select <- vars_train_df %>%
    as.data.frame(.)%>%
    names
  # Imputation for the training data
  if (which_imputation == 'train'){
    # Convert all variables to numeric and apply group-based imputation
    imputed_df <- vars_train_df %>%
      as.data.frame(.)%>%
      mutate(across(everything(), ~ as.numeric(.)))%>%
      group_by(get(group_var))%>%

      # Impute missing values for all columns except median_var using group-wise mean
      mutate(across(-(c(median_var)), 
                    ~ case_when(is.na(.) ~ mean(., na.rm = TRUE),
                                TRUE ~ .)))%>%
      
      # Impute missing values for median_var using group-wise median
      mutate(across(median_var, ~ case_when(is.na(.) ~ median(., na.rm = TRUE),
                                            TRUE ~ .)))%>%
      ungroup()%>%
      
      # Ensure columns are in the original order
      dplyr::select(vars_to_select)%>%
      as.matrix(.)
    
    # Imputation for the test data
  }else if(which_imputation == 'test'){
    # Calculate mean imputation values for all columns except median_var based on training data
    imputed_values_mean <- vars_train_df %>%
      as.data.frame(.)%>%
      mutate(across(everything(), ~as.numeric(.)))%>%
      group_by(get(group_var))%>%
      summarise(across(-c(median_var), ~mean(., na.rm = TRUE)))%>%
      dplyr::select(-`get(group_var)`)
    
    # Calculate median imputation values for median_var based on training data
    imputed_values_median <- vars_train_df %>%
      as.data.frame(.)%>%
      mutate(across(everything(), ~as.numeric(.)))%>%
      group_by(get(group_var))%>%
      summarise(across(median_var, ~median(., na.rm = TRUE)))%>%
      dplyr::rename(!!group_var := `get(group_var)`)
    
    # Apply the imputed values to the test data
    imputed_df <- vars_test_df %>%
      as.data.frame() %>%
      mutate(across(everything(), ~ as.numeric(.))) %>%
      # Left join to get the imputed mean values from the training data
      left_join(imputed_values_mean, by = group_var, suffix = c("", "_mean")) %>%
      # Left join to get the imputed median values from the training data
      left_join(imputed_values_median, by = group_var, suffix = c("", "_median")) %>%
      group_by(get(group_var)) %>%
      
      # Impute missing values in the test data using the pre-calculated mean values
      mutate(across(-(c(median_var, group_var, paste0(median_var, "_median"), ends_with('_mean'))), ~ case_when(
        is.na(.) ~ get(paste0(cur_column(), "_mean")),  # Use the pre-calculated mean
        TRUE ~ .
      ))) %>%
      
      # Impute missing values for median_var using the pre-calculated median values
      mutate(across(median_var, ~ case_when(
        is.na(.) ~ get(paste0(cur_column(), "_median")),  # Use the pre-calculated median
        TRUE ~ .
      ))) %>%
      ungroup() %>%
      
      # Ensure columns are in the original order
      dplyr::select(vars_to_select)%>%
      as.matrix(.)
  }
  
  # Return the imputed dataset either training or test data
  return(imputed_df)
}

align_cca_pairs <- function(X_fold, Y_fold, X_ref, Y_ref, verbose = FALSE) {
  K <- ncol(X_ref)
  
  # Compute similarity (absolute correlation) matrices
  sim_X <- abs(cor(X_fold, X_ref))
  sim_Y <- abs(cor(Y_fold, Y_ref))
  
  # Joint similarity matrix: average similarity for each component pair
  joint_sim <- (sim_X + sim_Y) / 2
  
  # Solve assignment problem (find best one-to-one mapping)
  assignment <- solve_LSAP(joint_sim, maximum = TRUE)
  
  aligned_X <- matrix(NA, nrow = nrow(X_fold), ncol = K)
  aligned_Y <- matrix(NA, nrow = nrow(Y_fold), ncol = K)
  
  for (k in 1:K) {
    j <- assignment[k]
    
    # Compute best sign for aligned pair
    sign_X <- sign(cor(X_fold[, j], X_ref[, k]))
    sign_Y <- sign(cor(Y_fold[, j], Y_ref[, k]))
    
    # Choose consistent sign (based on joint correlation direction)
    sign_flip <- ifelse(sign_X + sign_Y >= 0, 1, -1)
    
    aligned_X[, k] <- sign_flip * X_fold[, j]
    aligned_Y[, k] <- sign_flip * Y_fold[, j]
    
    if (verbose) {
      cat(sprintf("Aligning fold comp %d → ref comp %d | sign = %s | cor_X = %.3f | cor_Y = %.3f\n",
                  j, k, ifelse(sign_flip == 1, "+", "-"),
                  cor(X_fold[, j], X_ref[, k]),
                  cor(Y_fold[, j], Y_ref[, k])))
    }
  }
  
  return(list(X_aligned = aligned_X, Y_aligned = aligned_Y, 
              # Return or store which is the assignment so you could use it in 
              # permutation testing
              alignment_indices = assignment  # e.g., assignment[1] = 5 means fold comp1 → ref comp5
  ))
}




# Canonical Correlation Analysis (CCA) function
# This function computes the canonical correlation between two matrices, Y and X.
# It uses QR decomposition and Singular Value Decomposition (SVD) to find the canonical variables and correlations.
cca <- function(Y, X, R, S) {
  # N is the number of observations (rows) in the response matrix Y
  N <- nrow(Y)
  # Compute QR decomposition of Y (response matrix)
  qr_decomp <- qr(Y, LAPACK = TRUE)
  Qy <- qr.Q(qr_decomp)
  Ry <- qr.R(qr_decomp)
  iY <- qr_decomp$pivot
  # Compute QR decomposition of X (predictor matrix)
  qr_decomp <- qr(X, LAPACK = TRUE)
  Qx <- qr.Q(qr_decomp)
  Rx <- qr.R(qr_decomp)
  iX <- qr_decomp$pivot
  # Calculate the rank of matrices Y and X, and set K to the minimum of the two ranks
  K <- min(qr(Y)$rank, qr(X)$rank)
  # Perform Singular Value Decomposition (SVD) on the cross-product of the orthogonal matrices Qy and Qx
  svd_result <- svd(t(Qy) %*% Qx)
  D <- svd_result$d
  L <- svd_result$u
  M <- svd_result$v
  # Calculate the canonical correlations, ensuring values are between 0 and 1
  cc <- pmin(pmax(D[1:K], 0), 1)
  # Compute the canonical coefficients (A and B) for the response and predictor matrices
  A <- solve(Ry) %*% L[, 1:K] * sqrt(N - R)
  B <- solve(Rx) %*% M[, 1:K] * sqrt(N - S)
  # Rearrange A and B to match the row order of the original matrices Y and X
  A[iY, ] <- A
  B[iX, ] <- B
  return(list(A = A, B = B, cc = cc))
}

# Permutation-based Canonical Correlation Analysis (CCA) with Cross-Validation
# This function applies a permutation test within cross-validation to assess canonical correlation
# between the training and test data matrices for two sets of variables Y and X.

# ----------------------------------------------------------- @ edit for align
permcca_winkler_cv <- function(Y, X, Y_test_cv, X_test_cv, nP, Pset, Pset_test_cv,
                               ref_X_loadings_perm, ref_Y_loadings_perm){

#permcca_winkler_cv <- function(Y, X, Y_test_cv, X_test_cv, nP, Pset, Pset_test_cv){
  # Initialize the result matrix Z to NULL
  Z = NULL
  # Get the number of rows in the response (Y) and predictor (X) matrices
  Ny = nrow(Y)
  Nx = nrow(X)
  # Check if the number of rows in Y and X are the same
  if(Ny != Nx){
    print('Y and X do not have the same number of rows.')
  }
  # Set the total number of observations (N) for consistency
  N = Ny 
  rm(Ny, Nx)
  
  # Create identity matrices for the training data (Y and X)
  I = diag(N)
  
  # Use the same procedure for the test data matrices
  Ny_test = nrow(Y_test_cv)
  Nx_test = nrow(X_test_cv)
  # Check if the number of rows in the test data Y and X are the same
  if(Ny_test != Nx_test){
    print('Y and X do not have the same number of rows.')
  }
  # Set the total number of test observations (N_test) for consistency
  N_test = Ny_test 
  rm(Ny_test, Nx_test)
  # Create identity matrices for the test data
  I_test = diag(N_test)
  
  # Define residualization parameters (no residualization here, so we set to identity)
  Qz = I
  Qz_test = I_test
  
  # Center the test and training data
  # Center the test Y and X matrices based on the training data's mean values
  Y_test_cv = center_cv(Y, Y_test_cv, 'apply')
  Y  = center_cv(Y)
  
  # Apply the transformation matrices Qz to both Y and Y_test_cv
  Y <- t(Qz) %*% Y
  Y_test_cv <- t(Qz_test) %*% Y_test_cv
  
  # Recalculate the number of rows after centering
  Ny = nrow(Y)
  Ny_test = nrow(Y_test_cv)
  # Set residualization parameter R = 0, indicating no correction for nuisance variables
  R = 0
  
  # Initialize variables for weight matrices
  W = Z
  Qw = Qz
  Qw_test = Qz_test
  
  # Center the test and training X matrices
  X_test_cv = center_cv(X, X_test_cv, 'apply')
  X = center_cv(X)
  
  # Apply the transformation matrices Qz to both X and X_test_cv
  X <- t(Qz) %*% X
  X_test_cv <- t(Qz_test) %*% X_test_cv
  
  # Recalculate the number of rows after centering for X
  Nx = nrow(X)
  Nx_test = nrow(X_test_cv)
  # Set residualization parameter S = 0, indicating no correction for nuisance variables
  S  = 0
  
  # Perform initial Canonical Correlation Analysis (CCA) on the training data
  Y = Qz %*% Y
  X = Qw %*% X
  R = R
  S = S
  
  # Run CCA between the centered Y and X matrices
  ccresults = cca(Qz %*% Y, Qw %*% X, R, S)
  
  # Extract the canonical coefficients and correlations
  A = ccresults$A
  B = ccresults$B
  r = ccresults$cc
  
  # Get the number of canonical components (K)
  K <- length(r)
  # Compute the canonical variates for both Y and X
  U <- Y %*% cbind(A, Null(A))
  V <- X %*% cbind(B, Null(B))
  
  # Test Phase: Apply CCA to the training data and extract canonical coefficients
  # Perform canonical correlation on the training data using the canonical functions
  cca_result_withinperm <- cancor(Y, X)
  
  # Extract canonical coefficients (loadings) for both Y and X from the training CCA
  Y_loadings_perm <- cca_result_withinperm$xcoef
  X_loadings_perm <- cca_result_withinperm$ycoef
  
  # ------------------------------------------------------------ @ edit fo align
  Y_loadings_perm_raw <- cca_result_withinperm$xcoef
  X_loadings_perm_raw <- cca_result_withinperm$ycoef
  # Truncate all to common component dimension (safe)
  K <- min(ncol(X_loadings_perm_raw), ncol(Y_loadings_perm_raw), 
           ncol(ref_X_loadings_perm), ncol(ref_Y_loadings_perm))
  X_loadings_perm_raw <- X_loadings_perm_raw[, 1:K, drop = FALSE]
  Y_loadings_perm_raw <- Y_loadings_perm_raw[, 1:K, drop = FALSE]
  ref_X <- ref_X_loadings_perm[, 1:K, drop = FALSE]
  ref_Y <- ref_Y_loadings_perm[, 1:K, drop = FALSE]
  # Align the weights to reference
  aligned <- align_cca_pairs(X_loadings_perm_raw, Y_loadings_perm_raw, ref_X, ref_Y)
  X_loadings_perm <- aligned$X_aligned
  Y_loadings_perm <- aligned$Y_aligned
  alignment_indices <- aligned$alignment_indices
  reordering_map <- match(1:length(alignment_indices), alignment_indices)
  
  # Check if the canonical correlations from both methods match
  # If the results do not match, print an error message
  if (any(round(cca_result_withinperm$cor, digits = 5) != round(ccresults$cc, digits = 5))){
    print('There is a problem with the 2 different cc functions')
  }
  
  # Project the test data onto the canonical variates using the loadings from the training data
  Y_test_canonical_perm <- Y_test_cv %*% Y_loadings_perm
  X_test_canonical_perm <- X_test_cv %*% X_loadings_perm
  
  # Compute canonical variates for the test data
  U_test <- Y_test_cv %*% cbind(A, Null(A))
  V_test <- X_test_cv %*% cbind(B, Null(B))
  # Calculate canonical correlations for the holdout fold (test data)
  canonical_correlations_test <- diag(cor(Y_test_canonical_perm, X_test_canonical_perm))
  
  # --------------------------------------------------------------------------
  # Project the training data onto the canonical variates to test whether it produces
  # the same outcomes as the test data projections
  Y_train_canonical_perm <- Y %*% Y_loadings_perm
  X_train_canonical_perm <- X %*% X_loadings_perm
  
  # Calculate canonical correlations for the training set
  canonical_correlations_train <- diag(cor(Y_train_canonical_perm, X_train_canonical_perm))
  # -------------------------------------------------------------------------- #

  # Initialise counters and matrices for tracking changes across permutations
  cnt = matrix(0, nrow = 1, ncol = K);
  lW  = matrix(0, nrow = 1,ncol = K);
  lW_test  = matrix(0, nrow = 1,ncol = K);
  lW_test_weights_change = matrix(0, nrow = 1,ncol = K);
  cnt_test = matrix(0, nrow = 1, ncol = K);
  cnt_test_weights_change = matrix(0, nrow = 1, ncol = K)
  
  # Loop through each permutation
  for (p in 1:nP) {
    if (is.logical(Pset)) {
      # First permutation is no permutation
      if (p == 1) {
        idxY <- 1:Ny
        idxX <- 1:Nx
      } else {
        # Randomly shuffle the rows for permutation
        idxY <- sample(Ny)
        idxX <- sample(Nx)
      }
      # When Pset is provided, use it directly
    } else {
      idxY <- as.matrix(Pset[,p, drop = FALSE])
      idxX <- 1:Nx
    }
    if (is.logical(Pset_test_cv)) {
      # First permutation is no permutation
      if (p == 1) {
        idxY_test <- 1:Ny_test
        idxX_test <- 1:Nx_test
      } else {
        # Randomly shuffle the rows for the test data permutation
        idxY_test <- sample(Ny_test)
        idxX_test <- sample(Nx_test)
      }
    } else {
      # When Pset_test_cv is provided, use it directly
      idxY_test <- as.matrix(Pset_test_cv[,p, drop = FALSE])
      idxX_test <- 1:Nx_test
    }
    
    # Calculate the canonical correlation within the current cross-validation fold
    rperm_k1 <- cancor(Y[idxY,1:ncol(Y)], X[idxX, 1:ncol(X)])
    # Extract canonical coefficients (loadings) from the permuted training data
    Y_loadings_k1 <- rperm_k1$xcoef
    X_loadings_k1 <- rperm_k1$ycoef
    
    # ---------------------------------------------------------- @ edit for align
    # Extract raw permuted loadings
    Y_loadings_k1_raw <- rperm_k1$xcoef
    X_loadings_k1_raw <- rperm_k1$ycoef
    # Truncate safely to match ref dimension
    K <- min(ncol(Y_loadings_k1_raw), ncol(X_loadings_k1_raw), length(reordering_map))
    Y_loadings_k1_raw <- Y_loadings_k1_raw[, 1:K, drop = FALSE]
    X_loadings_k1_raw <- X_loadings_k1_raw[, 1:K, drop = FALSE]
    # Reorder permuted components using non-permuted alignment map
    Y_loadings_k1 <- Y_loadings_k1_raw[, reordering_map[1:K], drop = FALSE]
    X_loadings_k1 <- X_loadings_k1_raw[, reordering_map[1:K], drop = FALSE]
    
    # Test whether the weights have changed for the training data
    rperm_train_weights_change <- diag(cor(Y[idxY, 1:ncol(Y)] %*% Y_loadings_k1, 
                                           X[idxX, 1:ncol(X)] %*% X_loadings_k1))
    # Test if the test sample remains similar (compared to permuted weights)
    rperm_test_sample_changes <- diag(cor(Y_test_cv[idxY_test, 1:ncol(Y_test_cv)] %*% Y_loadings_perm,
                                          X_test_cv[idxX_test, 1:ncol(X_test_cv)]  %*% X_loadings_perm))
    rperm_test_weights_change <- diag(cor(Y_test_cv %*% Y_loadings_k1,
                                          X_test_cv %*% X_loadings_k1))

    # --------------------------------------------------------- @ edit for align
    if (p == 1) {
      cor_reordered <- cca_result_withinperm$cor[reordering_map[1:K]]
      if (any(round(cor_reordered, 5) != round(rperm_train_weights_change, 5))) {
        warning("Mismatch between reordered canonical correlations and projections. Check reordering.")
      }
    }
    
#    if (p == 1){
#      # Check if the canonical correlations between the permuted training data and the original data are identical
#      if (any(round(cca_result_withinperm$cor, digits = 5) != round(rperm_train_weights_change, digits = 5))){
#        print('There is a problem with the cc function between outside and inside the permutation')
#      }
#    }
    for (k in 1:K) {
      # This part is needed only for permutation testing of the training fold
      # A separate test is applied to the testing fold, so this block is specific to the training set
      rperm_all <- cca(Qz %*% U[idxY, k:ncol(U)], Qw %*% V[idxX, k:ncol(V)], R, S)
      rperm <- rperm_all$cc
      # Calculate Wilk's test statistic (used to test the significance of the canonical correlations)
      lWtmp <- -rev(cumsum(rev(log(1 - rperm^2))))
      lW[k] <- lWtmp[1]
    }
    
    # Store the results of the permutation test for the test fold
    lW_test <- rperm_test_sample_changes
    lW_test_weights_change <- rperm_test_weights_change

    if (p == 1) {
      # Store initial values for later comparison in the permutation process
      lW1 <- lW
      lW1_test <- lW_test
      lW1_test_weights_change <- lW_test_weights_change
    }
    # Update the counters for each permutation where the observed statistic is greater than or equal to the initial statistic
    cnt <- cnt + (lW >= lW1)
    cnt_test <- cnt_test + (lW_test >= lW1_test)
    cnt_test_weights_change <- cnt_test_weights_change  + (lW_test_weights_change >= lW1_test_weights_change)
  }
  # Calculate the permutation uncorrected p-values
  punc  = cnt/nP
  punc_test = cnt_test/nP
  punc_test_weights_change = cnt_test_weights_change/nP
  
  # Compute the family-wise error rate (FWER) as the cumulative maximum of the uncorrected p-values
  pfwer = cummax(punc) # pfwer
  
  # Store the canonical correlations for the test and training datasets
  cc = r               # canonical correlations
  cc_test = canonical_correlations_test 
  
  # Store the FWER for the test dataset
  pfwer_test = punc_test 
  pfwer_test_weights_change = punc_test_weights_change
  
  # Store the canonical weights for the left and right datasets
  cwl = A              # canonical weights (left)
  cwr = B              # canonical weights (right)
  
  # Store the canonical variables for the left and right datasets
  cvl = Qz%*%Y%*%A         # canonical variables (left)
  cvr = Qw%*%X%*%B         # canonical variables (right)
  
  # Return a list of results including p-values, canonical correlations, weights, and variables
  return(list(pfwer = pfwer, 
              cc = r,
              cc_test = cc_test,
              pfwer_test = pfwer_test,
              pfwer_test_weights_change = pfwer_test_weights_change,
              cwl = cwl,
              cwr = cwr,
              cvl = cvl,
              cvr = cvr))
}

# new CCA function:
run_cv_aligned <- function(X_mat,
                           Y_mat,
                           group_labels,
                           ref_X_load,
                           ref_Y_load,
                           n_folds = 10,
                           n_repeats = 5,
                           n_perm_fold = 5000) {
  
  K <- min(ncol(ref_X_load), ncol(ref_Y_load))
  cc_test_all <- matrix(NA, nrow = n_folds * n_repeats, ncol = K)
  fold_index  <- 0
  
  for (rep_idx in 1:n_repeats) {
    folds <- createFolds(sample(group_labels), k = n_folds, list = TRUE)
    
    for (fold_idx in 1:n_folds) {
      fold_index <- fold_index + 1
      
      train_id <- unlist(folds[-fold_idx])
      test_id  <- unlist(folds[fold_idx])
      
      # ---------- your imputation ----------
      X_train <- impute_cv(X_mat[train_id , ], X_mat[test_id , ],
                           'studygroup_num', 'uhr_1stdegree_re', 'train')
      X_test  <- impute_cv(X_mat[train_id , ], X_mat[test_id , ],
                           'studygroup_num', 'uhr_1stdegree_re', 'test')
      Y_train <- Y_mat[train_id , ]
      Y_test  <- Y_mat[test_id  , ]
      
      # ---------- CCA on training ----------
      cca_fold <- cancor(X_train, Y_train)
      
      # ---------- align to reference ----------
      aligned <- align_cca_pairs(cca_fold$xcoef[, 1:K, drop = FALSE],
                                 cca_fold$ycoef[, 1:K, drop = FALSE],
                                 ref_X_load[, 1:K, drop = FALSE],
                                 ref_Y_load[, 1:K, drop = FALSE])
      
      X_load_fold <- aligned$X_aligned
      Y_load_fold <- aligned$Y_aligned
      
      # ---------- test-set canonical variates ----------
      X_test_can <- X_test %*% X_load_fold
      Y_test_can <- Y_test %*% Y_load_fold
      cc_test_all[fold_index, ] <- diag(cor(X_test_can, Y_test_can))
      
      # ---------- OPTIONAL: keep any extra diagnostics ----------
      #   If you still need pfwer_test, weight-change checks etc.,
      #   call permcca_winkler_cv() here exactly once *without*
      #   its internal permutation (nP = 1) – or skip it entirely
      #   because global permutation will be outside.
    }
  }
  return(cc_test_all)
}


# Subset the DWI (diffusion-weighted imaging) measure by arranging and processing the dataset
df2shuffle <- cca_vars_rop_chr_cv %>% 
  column_to_rownames('PSN')%>%
  arrange(studygroup_num)

# Create an ordered dataframe for studygroup labels and PSN
df2order <- cca_vars_rop_chr_cv %>% 
  arrange(studygroup_num)%>%
  dplyr::select(PSN, studygroup_num)

# Merge the DTI dataset with the group labels and arrange the data accordingly
dti2shuffle <- cca_dti_rop_chr_cv %>%
  left_join(., df2order, by = 'PSN')%>%
  arrange(studygroup_num)%>%
  dplyr::select(-studygroup_num)%>%
  column_to_rownames('PSN')
  
# Calculate the number of samples in each group (Group 1 and Group 2)
group1_size <- df2shuffle %>% filter(studygroup_num %in% 0)%>%nrow()  # Group 1 size
group2_size <- df2shuffle %>% filter(studygroup_num %in% 1)%>%nrow()  # Group 2 size

# Set the number of folds and repetitions for cross-validation
# n_folds <- 10 
# n_repeats <- 5
n_folds <- 10 
n_repeats <- 5 

# Create a vector of group labels based on the group sizes
group_labels <- c(rep(0, group1_size), rep(1, group2_size))

# Initialize an empty list to store cross-validation results
cv_splits <- list()

# Set seed for reproducibility to ensure the same splits each time
set.seed(123)  # Set a seed for reproducibility

# Prepare the data for CCA (Canonical Correlation Analysis) 
X_orig <- df2shuffle %>%
  as.matrix(.) # Original cohort predictors
X_imputed <- df2shuffle %>%
  mutate(across(everything(), ~ as.numeric(.)))%>%
  group_by(studygroup_num)%>%
  mutate(across(-(c(uhr_1stdegree_re)), 
                ~ case_when(is.na(.) ~ mean(., na.rm = TRUE),
                            TRUE ~ .)),
         uhr_1stdegree_re = 
           case_when(is.na(uhr_1stdegree_re) ~ 
                       median(uhr_1stdegree_re, na.rm = TRUE),
                     TRUE ~ uhr_1stdegree_re))%>%
  ungroup()%>%
  as.matrix(.) # Original cohort predictors
Y_orig <- dti2shuffle  %>%
  as.matrix(.)# Original cohort outcomes

# Cross-validation setup
cv_results <- list()
cv_results_train <- list()
cv_permutation_significance_train <- list()
cv_permutation_significance_test_weights_change <- list()
cv_permutation_significance_test <- list()

# ------------------------------------------------------@inserted for alignment: 
# run overall CCA to align the weights:
# # @inserted for alignment: Calculate reference weights from overall CCA (before CV loop)
overall_cca <- cancor(X_imputed, Y_orig)
#ref_X_loadings <- overall_cca$xcoef
#ref_Y_loadings <- overall_cca$ycoef

# Reference weights from full sample
overall_cca <- cancor(X_imputed, Y_orig)
ref_X <- overall_cca$xcoef
ref_Y <- overall_cca$ycoef

cc_obs_mat  <- run_cv_aligned(X_imputed, Y_orig, group_labels,
                              ref_X, ref_Y,
                              n_folds, n_repeats)
cc_obs_mean <- colMeans(cc_obs_mat, na.rm = TRUE)

# shuffle the rows. now do this outside the CV. 
n_perm <- 1000            # <- set as needed
K      <- length(cc_obs_mean)
null_dist <- matrix(NA, nrow = n_perm, ncol = K)

for (p in 1:n_perm) {
  cat("Permutation", p, "/", n_perm, "\n")
  
  # Shuffle *rows* of Y (or X) once
  Y_perm <- Y_orig[sample(nrow(Y_orig)), ]
  
  # Run the same CV with aligned weights, but on permuted Y
  cc_perm_mat <- run_cv_aligned(X_imputed, Y_perm, group_labels,
                                ref_X, ref_Y,
                                n_folds, n_repeats)
  
  null_dist[p, ] <- colMeans(cc_perm_mat, na.rm = TRUE)
}

emp_p <- sapply(1:K, function(k)
  (sum(null_dist[, k] >= cc_obs_mean[k]) + 1) / (n_perm + 1))
print(emp_p)



# Loop over each shuffle
for (repeat_idx in 1:n_repeats) {
  # Shuffle the data and create folds
  shuffled_data <- data.frame(id = 1:nrow(df2shuffle), group = sample(group_labels))
  folds <- createFolds(shuffled_data$group, k = n_folds, list = TRUE, returnTrain = FALSE)
  
  # Store fold-wise canonical correlations
  fold_canonical_correlations <- vector("list", n_folds)
  fold_canonical_correlations_train <- vector("list", n_folds)
  pfwer_rc_result <- vector("list", n_folds)
  pfwer_rc_result_test <- vector("list", n_folds)
  pfwer_rc_result_test_weights_change <- vector("list", n_folds)
  
  # Loop over each fold
  for (fold_idx in 1:n_folds) {
    # Split data into training and holdout (test)
    train_indices <- unlist(folds[-fold_idx])
    test_indices <- unlist(folds[fold_idx])
    
    # Training data
    X_train_before_impute <- X_orig[train_indices, ]
    Y_train <- Y_orig[train_indices, ]
    
    # Holdout data
    X_test_before_impute <- X_orig[test_indices, ]
    Y_test <- Y_orig[test_indices, ]
    
    X_train <- impute_cv(X_train_before_impute, X_test_before_impute,
                         'studygroup_num', 'uhr_1stdegree_re', 'train') 
    X_test <- impute_cv(X_train_before_impute, X_test_before_impute,
                        'studygroup_num', 'uhr_1stdegree_re', 'test') 
    
    # Apply CCA to the training data
    cca_result <- cancor(X_train, Y_train)
    
    perms_training <- create_permutation(X_train %>% as.data.frame(.), 5000)
    perms_test <- create_permutation(X_test %>% as.data.frame(.), 5000)
    
    # Extract canonical coefficients
    X_loadings <- cca_result$xcoef
    Y_loadings <- cca_result$ycoef
    
    # --------------------------------------------------@inserted for alignment:
    # Enforce equal number of components
    K <- min(ncol(X_loadings), ncol(Y_loadings), ncol(ref_X_loadings), ncol(ref_Y_loadings))
    # Truncate to shared component space
    X_loadings <- X_loadings[, 1:K, drop = FALSE]
    Y_loadings <- Y_loadings[, 1:K, drop = FALSE]
    ref_X_loadings <- ref_X_loadings[, 1:K, drop = FALSE]
    ref_Y_loadings <- ref_Y_loadings[, 1:K, drop = FALSE]
    # Align fold weights to overall reference weights
    aligned <- align_cca_pairs(X_loadings, Y_loadings, ref_X_loadings, ref_Y_loadings, verbose = TRUE)
    X_loadings_aligned <- aligned$X_aligned
    Y_loadings_aligned <- aligned$Y_aligned
    # test the different outputs:
    X_canonical_nonaligned <- X_train %*% X_loadings
    Y_canonical_nonaligned <- Y_train %*% Y_loadings
    canonical_correlations_nonaligned <- diag(cor(X_canonical_nonaligned, Y_canonical_nonaligned))
    X_canonical_aligned <- X_train %*% X_loadings_aligned
    Y_canonical_aligned <- Y_train %*% Y_loadings_aligned
    canonical_correlations_aligned <- diag(cor(X_canonical_aligned, Y_canonical_aligned))
    
    # Project the test data onto the canonical variates
    X_test_canonical <- X_test %*% X_loadings
    Y_test_canonical <- Y_test %*% Y_loadings
    # Calculate canonical correlations for the holdout fold
    canonical_correlations_test <- diag(cor(X_test_canonical, Y_test_canonical))
    
    # --------------------------------------------------@inserted for alignment:
    # Project the test data onto the canonical variates
    X_test_canonical_aligned <- X_test %*% X_loadings_aligned
    Y_test_canonical_aligned <- Y_test %*% Y_loadings_aligned
    canonical_correlations_test <- diag(cor(X_test_canonical_aligned, Y_test_canonical_aligned))
    
    # Store the canonical correlations for this fold
    fold_canonical_correlations[[fold_idx]] <- canonical_correlations_test

    # Project the train data onto the canonical variates
    X_train_canonical <- X_train %*% X_loadings
    Y_train_canonical <- Y_train %*% Y_loadings
    # Calculate canonical correlations for the training fold
    canonical_correlations_train <- diag(cor(X_train_canonical, Y_train_canonical))
   
    # --------------------------------------------------@inserted for alignment:
    # Project the train data onto the canonical variates
    X_train_canonical_aligned <- X_train %*% X_loadings_aligned
    Y_train_canonical_aligned <- Y_train %*% Y_loadings_aligned
    # Calculate canonical correlations for the training fold
    canonical_correlations_train <- diag(cor(X_train_canonical_aligned, Y_train_canonical_aligned))
    
    # Store the canonical correlations for this fold
    fold_canonical_correlations_train[[fold_idx]] <- canonical_correlations_train    
    
    print(paste0('Here we will run again the permutation test, we are in fold: ',
                 fold_idx, '/', n_folds, ' and repetition: ', repeat_idx, '/', n_repeats, sep = ''))
    perm_result_train_test <- permcca_winkler_cv(Y = as.matrix(X_train),
                                                 X = as.matrix(Y_train),
                                                 Y_test_cv = as.matrix(X_test),
                                                 X_test_cv = as.matrix(Y_test),
                                                 nP = 5000,
                                                 Pset = perms_training,
                                                 Pset_test_cv = perms_test,
                                                 ref_X_loadings_perm = ref_Y_loadings, 
                                                 ref_Y_loadings_perm = ref_X_loadings)
    pfwer_rc_result[[fold_idx]] <- perm_result_train_test$pfwer
    pfwer_rc_result_test[[fold_idx]] <- perm_result_train_test$pfwer_test
    pfwer_rc_result_test_weights_change[[fold_idx]] <- perm_result_train_test$pfwer_test_weights_change
    
    number_to_keep <- min(ncol(X_train), ncol(Y_train))
  }
  
  # Store results for this shuffle
  cv_results[[repeat_idx]] <- fold_canonical_correlations
  cv_results_train[[repeat_idx]] <- fold_canonical_correlations_train
  cv_permutation_significance_test[[repeat_idx]] <- pfwer_rc_result_test
  cv_permutation_significance_test_weights_change[[repeat_idx]] <- pfwer_rc_result_test_weights_change
  cv_permutation_significance_train[[repeat_idx]] <- pfwer_rc_result
}

output <- list(cv_results_test = cv_results, 
               cv_results_train = cv_results_train,
               cv_permutation_train = cv_permutation_significance_train,
               cv_permutation_test = cv_permutation_significance_test,
               cv_permutation_test_weights_change = cv_permutation_significance_test_weights_change)

# saveRDS(output, 'simulated_data/CV_output_aligned.Rds')

command = 'calculated'
if (command == 'loaded'){
  loaded_output <- readRDS('simulated_data/CV_output.Rds')
  }else if(command == 'calculated'){
  loaded_output <- output
}
cv_results = loaded_output$cv_results_test 
cv_results_train = loaded_output$cv_results_train
cv_permutation_significance_train = loaded_output$cv_permutation_train
cv_permutation_significance_test = loaded_output$cv_permutation_test
cv_permutation_significance_test_weights_change = loaded_output$cv_permutation_test_weights_change

# Initialize a list to store the canonical correlations for each position (1 to 15)
canonical_correlations_list <- vector("list", 15)
canonical_correlations_list_permutation <- vector("list", 15)
canonical_correlations_list_permutation_test <- vector("list", 15)
canonical_correlations_list_permutation_test_weights_change <- vector("list", 15)
canonical_correlations_list_train <- vector("list", 15)
# calculate everything also for train
for (repeat_idx in 1:n_repeats) {
  # Loop through each fold within the current shuffle
  for (fold_idx in 1:n_folds) {
    # Extract all 15 canonical correlations for this fold
    for (i in 1:15) {
      canonical_correlations_list_train[[i]] <- c(canonical_correlations_list_train[[i]], cv_results_train[[repeat_idx]][[fold_idx]][i])
      canonical_correlations_list[[i]] <- c(canonical_correlations_list[[i]], cv_results[[repeat_idx]][[fold_idx]][i])
      canonical_correlations_list_permutation[[i]] <- c(canonical_correlations_list_permutation[[i]],
                                                        cv_permutation_significance_train[[repeat_idx]][[fold_idx]][i])
      canonical_correlations_list_permutation_test[[i]] <- c(canonical_correlations_list_permutation_test[[i]],
                                                        cv_permutation_significance_test[[repeat_idx]][[fold_idx]][i])
      canonical_correlations_list_permutation_test_weights_change[[i]] <- c(canonical_correlations_list_permutation_test_weights_change[[i]],
                                                             cv_permutation_significance_test_weights_change[[repeat_idx]][[fold_idx]][i])
    }
  }
}

# Calculate the mean for each of the 15 canonical correlations across all folds and shuffles
mean_canonical_correlations <- sapply(canonical_correlations_list, mean)
median_canonical_correlations <- sapply(canonical_correlations_list, median)
sd_canonical_correlations <- sapply(canonical_correlations_list, sd)
std_error_test <- sd_canonical_correlations / sqrt(n_repeats * n_folds)
# Calculate the 95% confidence interval
conf_interval_lower <- mean_canonical_correlations - 1.96 * std_error_test
conf_interval_higher <- mean_canonical_correlations + 1.96 * std_error_test
# Calculate the mean for each of the 15 canonical correlations across all folds and shuffles
mean_canonical_correlations_train <- sapply(canonical_correlations_list_train, mean)
median_canonical_correlations_train <- sapply(canonical_correlations_list_train, median)
sd_canonical_correlations_train <- sapply(canonical_correlations_list_train, sd)
# Calculate the standard error
std_error_train <- sd_canonical_correlations_train / sqrt(n_repeats * n_folds)
# Calculate the 95% confidence interval
conf_interval_lower_train <- mean_canonical_correlations_train - 1.96 * std_error_train
conf_interval_higher_train <- mean_canonical_correlations_train + 1.96 * std_error_train
# Permutation: Calculate the mean for each of the 15 canonical correlations across all folds and shuffles
mean_canonical_correlations_perm <- sapply(canonical_correlations_list_permutation, mean)
median_canonical_correlations_perm <- sapply(canonical_correlations_list_permutation, median)
sd_canonical_correlations_perm <- sapply(canonical_correlations_list_permutation, sd)
std_error_perm <- sd_canonical_correlations_perm / sqrt(n_repeats * n_folds)
# Calculate the 95% confidence interval
conf_interval_lower_perm <- mean_canonical_correlations_perm - 1.96 * std_error_perm
conf_interval_higher_perm <- mean_canonical_correlations_perm + 1.96 * std_error_perm
# Permutation: Calculate the mean for each of the 15 canonical correlations across all folds and shuffles
mean_canonical_correlations_perm_test <- sapply(canonical_correlations_list_permutation_test, mean)
median_canonical_correlations_perm_test <- sapply(canonical_correlations_list_permutation_test, median)
sd_canonical_correlations_perm_test <- sapply(canonical_correlations_list_permutation_test, sd)
std_error_perm_test <- sd_canonical_correlations_perm_test / sqrt(n_repeats * n_folds)
# Calculate the 95% confidence interval
conf_interval_lower_perm_test <- mean_canonical_correlations_perm_test - 1.96 * std_error_perm_test
conf_interval_higher_perm_test <- mean_canonical_correlations_perm_test + 1.96 * std_error_perm_test
# Permutation: Calculate the mean for each of the 15 canonical correlations across all folds and shuffles
mean_canonical_correlations_perm_test_weights_change <- sapply(canonical_correlations_list_permutation_test_weights_change, mean)
median_canonical_correlations_perm_test_weights_change <- sapply(canonical_correlations_list_permutation_test_weights_change, median)
sd_canonical_correlations_perm_test_weights_change <- sapply(canonical_correlations_list_permutation_test_weights_change, sd)
std_error_perm_test_weights_change <- sd_canonical_correlations_perm_test_weights_change / sqrt(n_repeats * n_folds)
# Calculate the 95% confidence interval
conf_interval_lower_perm_test_weights_change <- mean_canonical_correlations_perm_test_weights_change - 1.96 * std_error_perm_test_weights_change
conf_interval_higher_perm_test_weights_change <- mean_canonical_correlations_perm_test_weights_change + 1.96 * std_error_perm_test_weights_change

# Print the result in a readable format
for (i in 1:2) {
  cat("Test folds: Mean of Canonical Correlation", i, ":", mean_canonical_correlations[i], '(', sd_canonical_correlations[i], ')', "\n")
  cat("Test folds: Median of Canonical Correlation", i, ":", median_canonical_correlations[i], "\n")
  cat("Test folds: Confidence interval", i, ":", conf_interval_lower[i], '-', conf_interval_higher[i], "\n")
  cat("Training folds: Mean of Canonical Correlation", i, ":", mean_canonical_correlations_train[i], '(', sd_canonical_correlations_train[i], ')', "\n")
  cat("Training folds: Median of Canonical Correlation", i, ":", median_canonical_correlations_train[i], "\n")
  cat("Training folds: Confidence interval", i, ":", conf_interval_lower_train[i], '-', conf_interval_higher_train[i], "\n")
  cat("PERMUTATION: Training folds: Mean of Canonical Correlation", i, ":", mean_canonical_correlations_perm[i], '(', sd_canonical_correlations_perm[i], ')', "\n")
  cat("PERMUTATION: Training folds: Median of Canonical Correlation", i, ":", median_canonical_correlations_perm[i], "\n")
  cat("PERMUTATION: Training folds: Confidence interval", i, ":", conf_interval_lower_perm[i], '-', conf_interval_higher_perm[i], "\n")
  cat("PERMUTATION: Test folds: Mean of Canonical Correlation", i, ":", mean_canonical_correlations_perm_test[i], '(', sd_canonical_correlations_perm_test[i], ')', "\n")
  cat("PERMUTATION: Test folds: Median of Canonical Correlation", i, ":", median_canonical_correlations_perm_test[i], "\n")
  cat("PERMUTATION: Test folds: Confidence interval", i, ":", conf_interval_lower_perm_test[i], '-', conf_interval_higher_perm_test[i], "\n")
  cat("PERMUTATION weights change: Test folds: Mean of Canonical Correlation", i, ":", mean_canonical_correlations_perm_test_weights_change[i], '(', sd_canonical_correlations_perm_test_weights_change[i], ')', "\n")
  cat("PERMUTATION weights change: Test folds: Median of Canonical Correlation", i, ":", median_canonical_correlations_perm_test_weights_change[i], "\n")
  cat("PERMUTATION weights change: Test folds: Confidence interval", i, ":", conf_interval_lower_perm_test_weights_change[i], '-', conf_interval_higher_perm_test_weights_change[i], "\n")
}