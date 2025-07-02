library(tidyverse)
library(CCA)
library(caret)
library(dplyr)
library(clue)

cca_dti_rop_chr_cv <- read.csv('~/Documents/Harvard/projects/DTI_Pronia/in_between_data/cca_dti_rop_chr_for_cv_V18.csv')
cca_vars_rop_chr_cv <- read.csv('~/Documents/Harvard/projects/DTI_Pronia/in_between_data/cca_vars_rop_chr_for_cv_NA_V18.csv')

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

center <- function(X) {
  icte <- apply(diff(X, 1, 1)^2, 2, sum) == 0
  X <- sweep(X, 2, colMeans(X), "-")
  X <- X[, !icte]
  return(X)
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

# i will save U_train, V_train and U_test, V_test
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
    #sign_flip <- ifelse(sign_X + sign_Y >= 0, 1, -1)
    flip_score <- cor(X_fold[, j], X_ref[, k]) + cor(Y_fold[, j], Y_ref[, k])
    sign_flip <- sign(flip_score)
    
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
# permcca_winkler function
permcca_winkler <- function(Y, X, nP, Pset){
  # Permutation inference for canonical correlation analysis (CCA)
  
  # Inputs:
  # - Y        : Left set of variables, size N by P.
  # - X        : Right set of variables, size N by Q.
  # - nP       : An integer representing the number
  #              of permutations.
  #              Default is 1000 permutations.
  # - Pset     : Predefined set of permutations (e.g., that respect
  #                                              %              exchangeability blocks). For information on how to
  #              generate these, see:
  #              https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM
  #              If a selection matrix is provided (for the Theil method),
  #              Pset will have to have fewer rows than the original N, i.e.,
  #              it will have as many rows as the effective number of
  #              subjects determined by the selection matrix.
  #
  # Outputs:
  # - p   : p-values, FWER corrected via closure.
  # - r   : Canonical correlations.
  # - A   : Canonical coefficients, left side.
  # - B   : Canonical coefficients, right side.
  # - U   : Canonical variables, left side.
  # - V   : Canonical variables, right side.
  #
  # ___________________________________________
  # AM Winkler, O Renaud, SM Smith, TE Nichols
  # NIH - Univ. of Geneva - Univ. of Oxford
  # Mar/2020 (first version)
  # Apr/2021 (this version)
  Z = NULL
  Ny = nrow(Y)
  Nx = nrow(X)
  if(Ny != Nx){
    print('Y and X do not have the same number of rows.')
  }
  N = Ny 
  rm(Ny, Nx)
  I = diag(N)
  
  # since we do not residualize
  Qz = I
  # when numbers are very small ~ e-15 the sign is different in R
  Y  = center(Y)
  Y <- t(Qz) %*% Y
  
  Ny = nrow(Y)
  # since I do not write this with the possibility to correct for nuisance
  # variates I specify that R = 0
  R = 0
  
  W = Z
  Qw = Qz
  
  X = center(X)
  X <- t(Qz) %*% X
  
  Nx = nrow(X)
  # since I do not write this with the possibility to correct for nuisance
  # variates I specify that S = 0
  S  = 0
  
  # Initial CCA
  Y = Qz %*% Y
  X = Qw %*% X
  R = R
  S = S
  ccresults = cca(Qz %*% Y, Qw %*% X, R, S) 
  A = ccresults$A
  B = ccresults$B
  r = ccresults$cc
  
  K <- length(r)
  U <- Y %*% cbind(A, Null(A))
  V <- X %*% cbind(B, Null(B))
  
  # Initialise counter
  cnt = matrix(0, nrow = 1, ncol = K);
  lW  = matrix(0, nrow = 1,ncol = K);
  
  # For each permutation
  for (p in 1:nP) {
    #print(paste0('Permutation: ', p, '/', nP))
    
    # If user didn't supply a set of permutations,
    # permute randomly both Y and X.
    # Otherwise, use the permutation set to shuffle
    # one side only.
    if (is.logical(Pset)) {
      # First permutation is no permutation
      if (p == 1) {
        idxY <- 1:Ny
        idxX <- 1:Nx
      } else {
        idxY <- sample(Ny)
        idxX <- sample(Nx)
      }
    } else {
      # let's see whether it is a problem that the two are slightly 
      # different types
      idxY <- as.matrix(Pset[,p, drop = FALSE])
      idxX <- 1:Nx
    }
    
    # For each canonical variable up to K canonical components
    for (k in 1:K) {
      #cat(k, ' ')
      # Run canonical correlation analysis (CCA) on permuted data
      # Using only the current and remaining components from U and V
      rperm_all <- cca(Qz %*% U[idxY, k:ncol(U)], Qw %*% V[idxX, k:ncol(V)], R, S)
      
      # Extract canonical correlations from permuted CCA result
      rperm <- rperm_all$cc
      
      # Compute log-Wilks Lambda statistic from permuted canonical correlations
      # This statistic accumulates evidence across components (from k to end)
      lWtmp <- -rev(cumsum(rev(log(1 - rperm^2))))
      
      # Store the log-Wilks Lambda statistic for the current component k
      lW[k] <- lWtmp[1]
    }
    
    # On first permutation, store the observed log-Wilks Lambda values
    if (p == 1) {
      lW1 <- lW
    }
    
    # Count how many times permuted stats exceed or equal observed values
    # (used to compute uncorrected and family-wise error rate p-values)
    cnt <- cnt + (lW >= lW1)
    #cat('\n')
  }
  punc  = cnt/nP
  pfwer = cummax(punc) # pfwer
  cc = r               # canonical correlations
  cwl = A              # canonical weights (left)
  cwr = B              # canonical weights (right)
  cvl = Qz%*%Y%*%A         # canonical variables (left)
  cvr = Qw%*%X%*%B         # canonical variables (right)
  return(list(pfwer = pfwer, 
              cc = r, 
              cwl = cwl,
              cwr = cwr,
              cvl = cvl,
              cvr = cvr))
}
# ------------------------------------------------------------------
# 1.  Create one reproducible fold structure to re-use everywhere
# ------------------------------------------------------------------
create_cv_folds <- function(group_labels,
                            n_folds   = 10,
                            n_repeats = 5,
                            seed      = 123) {
  if (!is.null(seed)) set.seed(seed)
  
  fold_out <- vector("list", n_repeats)
  
  for (rep in seq_len(n_repeats)) {
    ## stratified folds (caret) 
    fold_out[[rep]] <- createFolds(group_labels,
                                   k = n_folds,
                                   list = TRUE,
                                   returnTrain = FALSE)
  }
  return(fold_out)   # list[[rep]][[fold]] gives test indices
}

# ------------------------------------------------------------------
# 2.  Cross-validated aligned CCA that re-uses the fold structure
# ------------------------------------------------------------------
run_cv_aligned <- function(X_mat,
                           Y_mat,
                           folds_list,        # <- fixed folds
                           ref_X_load,
                           ref_Y_load,
                           U_train = NULL,
                           V_train = NULL,
                           U_test = NULL,
                           V_test = NULL,
                           obs_perm = NULL,
                           perm_idx = NULL,
                           k_components) {
  
  n_repeats <- length(folds_list)
  n_folds   <- length(folds_list[[1]])
  #K         <- min(ncol(ref_X_load), ncol(ref_Y_load))
  K         <- k_components
  
  cc_test_all <- matrix(NA, nrow = n_repeats * n_folds, ncol = K)
  cc_train_all <- matrix(NA, nrow = n_repeats * n_folds, ncol = K)
  lW_train_all <- matrix(NA, nrow = n_repeats * n_folds, ncol = K)
  
  if (obs_perm == 'obs'){
    U_test_all <- vector("list", n_repeats * n_folds)
    V_test_all <- vector("list", n_repeats * n_folds)
    U_train_all <- vector("list", n_repeats * n_folds)
    V_train_all <- vector("list", n_repeats * n_folds)
  }

  fold_index  <- 0
  alignment_history <- matrix(0, nrow = n_repeats * n_folds, ncol = K)
  
  for (rep_idx in seq_len(n_repeats)) {
    folds <- folds_list[[rep_idx]]
    
    for (fold_idx in seq_len(n_folds)) {
      fold_index <- fold_index + 1
      
      test_id  <- folds[[fold_idx]]
      train_id <- setdiff(seq_len(nrow(X_mat)), test_id)
      
      if (obs_perm == 'obs'){
        
        # ---------- Imputation ----------
        X_train_raw <- impute_cv(X_mat[train_id, ], X_mat[test_id, ],
                                 'studygroup_num', 'uhr_1stdegree_re', 'train')
        X_test_raw  <- impute_cv(X_mat[train_id, ], X_mat[test_id, ],
                                 'studygroup_num', 'uhr_1stdegree_re', 'test')
        Y_train_raw <- Y_mat[train_id, ]
        Y_test_raw  <- Y_mat[test_id, ]
        
        # ---------- Centre using your centre_cv() ----------
        X_train <- center_cv(X_train_raw)                       # centres & trims const-cols
        X_test  <- center_cv(X_train_raw, X_test_raw, 'apply')  # same means/cols
        Y_train <- center_cv(Y_train_raw)
        Y_test  <- center_cv(Y_train_raw, Y_test_raw, 'apply')
      
        # ---------- CCA on training ----------
        train_result <- permcca_winkler(Y_train, X_train, 1, TRUE)  # No permutation here
        A <- train_result$cwr  # Canonical weights (Y side)
        B <- train_result$cwl  # Canonical weights (X side)
        
        #cca_fold <- cancor(X_train, Y_train)
        cc_train_all[fold_index, ] <- train_result$cc[1:K]
        
        # ----- Align to reference, then reorder columns -----
        #aligned   <- align_cca_pairs(train_result$cwr[, 1:K, drop=FALSE],
        #                             train_result$cwl[, 1:K, drop=FALSE],
        #                             ref_X_load[, 1:K, drop=FALSE],
        #                             ref_Y_load[, 1:K, drop=FALSE])
        #
        #X_load_raw <- aligned$X_aligned            # sign-flipped but still in fold order
        #Y_load_raw <- aligned$Y_aligned
        #reorder    <- aligned$alignment_indices    # e.g. reorder[1] = 5 means fold-col5 → ref-col1
        #alignment_history[fold_index, ] <- reorder
        
        # Reorder to reference order
        #X_load <- X_load_raw[, reorder, drop = FALSE]
        #Y_load <- Y_load_raw[, reorder, drop = FALSE]
        
        # reorder columns to reference order
        
        # ---------- Test-set canonical variates ----------
        #X_test_can <- X_test %*% X_load
        #Y_test_can <- Y_test %*% Y_load
        ## Compute aligned canonical variates on training set
        #X_train_can <- X_train %*% X_load
        #Y_train_can <- Y_train %*% Y_load
        # ---------- Just shortly test without alignment
        X_test_can <- X_test %*% A
        Y_test_can <- Y_test %*% B
        # Compute aligned canonical variates on training set
        X_train_can <- X_train %*% A
        Y_train_can <- Y_train %*% B
        cc_vec <- diag(cor(X_test_can, Y_test_can))[1:K]
        cc_test_all[fold_index, ] <- cc_vec
        # Canonical correlations (aligned)
        cc_train_all[fold_index, ] <- diag(cor(X_train_can, Y_train_can))[1:K]
        lW_train_all[fold_index, ] <- -rev(cumsum(rev(log(1 - diag(cor(X_train_can, Y_train_can))[1:K]^2))))
        
        rownames(X_train_can) <- train_id
        rownames(Y_train_can) <- train_id
        rownames(X_test_can)  <- test_id
        rownames(Y_test_can)  <- test_id
        
        U_train_all[[fold_index]] <- data.frame(X_train_can)
        V_train_all[[fold_index]] <- data.frame(Y_train_can)
        U_test_all[[fold_index]]  <- data.frame(X_test_can)
        V_test_all[[fold_index]]  <- data.frame(Y_test_can)
        
      }
      else if(obs_perm == 'perm'){
        
        U_test_fold <- as.matrix(U_test[[fold_index]])
        V_test_fold <- as.matrix(V_test[[fold_index]])
        U_train_fold <- as.matrix(U_train[[fold_index]])
        V_train_fold <- as.matrix(V_train[[fold_index]])
        
        # ---------- CCA on training ----------
        train_id_perm <- as.character(perm_idx[perm_idx %in% rownames(U_train_fold)])
        test_id_perm  <- as.character(perm_idx[perm_idx %in% rownames(U_test_fold)])
        
        lW  = matrix(0, nrow = 1,ncol = K);
        r_perm  = matrix(0, nrow = 1,ncol = K);
        
        for (k in 1:K) {
          rperm_all <- cca(U_train_fold[as.character(train_id_perm), k:ncol(U_train_fold)],
                           V_train_fold[as.character(train_id), k:ncol(V_train_fold)], 0,0)
          
          # Extract canonical correlations from permuted CCA result
          rperm <- rperm_all$cc[1:K]
            
          # Compute log-Wilks Lambda statistic from permuted canonical correlations
          # This statistic accumulates evidence across components (from k to end)
          lWtmp <- -rev(cumsum(rev(log(1 - rperm^2))))
            
          # Store the log-Wilks Lambda statistic for the current component k
          lW[k] <- lWtmp[1]
          r_perm[k] <- rperm[1]
        }
        # save lW for each fold and permutation
        cc_vec <- diag(cor(U_test_fold[as.character(test_id_perm),], 
                           V_test_fold[as.character(test_id),]))[1:K]
        cc_test_all[fold_index,] <- cc_vec[1:K]
        cc_train_all[fold_index,] <- r_perm[1:K]
        lW_train_all[fold_index,] <- lW[1:K]
      }
    }
  }
  if (obs_perm == 'obs'){
    return(list(cc_test_all = cc_test_all,
                cc_train_all = cc_train_all,
                alignment_history = alignment_history,
                lW_train_all = lW_train_all,
                U_test_all = U_test_all,
                U_train_all = U_train_all,
                V_test_all = V_test_all,
                V_train_all = V_train_all))   # rows = folds × repeats ; cols = components
  }
  else if (obs_perm == 'perm'){
    return(list(cc_test_all = cc_test_all,
                cc_train_all = cc_train_all,
                lW_train_all = lW_train_all,
                alignment_history = alignment_history))   # rows = folds × repeats ; cols = components
  }
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

# ------------------------------------------------------@inserted for alignment: 
# run overall CCA to align the weights:
# # @inserted for alignment: Calculate reference weights from overall CCA (before CV loop)
# Reference weights from full sample
# Step 1: Run full-sample CCA on unshuffled data
# Create fixed folds outside for consistency across permutations
# Build fixed folds once
create_cv_folds_without_group <- function(n_obs,
                                          n_folds = 10,
                                          n_repeats = 5,
                                          seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  
  fold_out <- vector("list", n_repeats)
  
  for (rep in seq_len(n_repeats)) {
    indices <- sample(n_obs)  # shuffle row indices
    folds <- split(indices, cut(seq_along(indices), breaks = n_folds, labels = FALSE))
    
    # Name the folds as Fold01, Fold02, ...
    names(folds) <- sprintf("Fold%02d", seq_len(n_folds))
    
    fold_out[[rep]] <- folds
  }
  
  return(fold_out)  # Same structure as caret::createFolds
}



#folds_list <- create_cv_folds(group_labels,
#                              n_folds   = n_folds,
#                              n_repeats = n_repeats,
#                              seed      = 123)
folds_list <- create_cv_folds_without_group(n_obs = nrow(X_orig),
                              n_folds   = n_folds,
                              n_repeats = n_repeats,
                              seed      = 123)

# Reference CCA (unshuffled)
overall_cca <- permcca_winkler(X_imputed, Y_orig, 1, TRUE)
#overall_cca <- cancor(X_imputed, Y_orig)
#ref_X <- overall_cca$xcoef
#ref_Y <- overall_cca$ycoef
ref_X <- overall_cca$cwl
ref_Y <- overall_cca$cwr

# Observed CV performance:
res_obs_new <- run_cv_aligned(X_orig, Y_orig, folds_list, ref_X, ref_Y,
                              obs_perm = 'obs', k_components = 3)
colMeans(res_obs_new$cc_test_all)
colMeans(res_obs_new$cc_train_all)
colMeans(res_obs_new$lW_train_all)

U_test_obs  <- res_obs_new$U_test_all
U_train_obs <- res_obs_new$U_train_all
V_test_obs  <- res_obs_new$V_test_all
V_train_obs <- res_obs_new$V_train_all

# Permuted performances
n_perm <- 1000 
perm_results <- vector("list", n_perm)
set.seed(123)
perm_idx_list <- replicate(n_perm, sample(seq_len(nrow(X_orig))), simplify = FALSE)

#for (p in 1:n_perm) {
for (p in 1:n_perm) {
  print(paste0("Permutation number: ", p, ' of overall: ', n_perm, 
               ' permutations.', sep = ''))
  perm_results[[p]] <- run_cv_aligned(X_orig, Y_orig, folds_list, ref_X, ref_Y,
                                      U_train = U_train_obs, V_train = V_train_obs,
                                      U_test = U_test_obs, V_test = V_test_obs,
                                      obs_perm = 'perm', perm_idx = perm_idx_list[[p]],
                                      k_components = 3)
}

#X_mat <- X_orig
#Y_mat <- Y_orig
#folds_list <- folds_list
#ref_X_load <- ref_X
#ref_Y_load <- ref_Y
#U_train = U_train_obs 
#V_train = V_train_obs
#U_test = U_test_obs
#V_test = V_test_obs
#obs_perm = 'perm'
#perm_idx = perm_idx_list[[p]]


# Dimensions
K <- ncol(res_obs_new$cc_test_all)

# Initialize matrices
cc_test_null <- matrix(NA, nrow = n_perm, ncol = K)
lW_train_null <- matrix(NA, nrow = n_perm, ncol = K)
cc_train_null <- matrix(NA, nrow = n_perm, ncol = K)

# Fill matrices with means from each permutation
for (p in 1:n_perm) {
  cc_test_null[p, ] <- colMeans(perm_results[[p]]$cc_test_all)
  lW_train_null[p, ] <- colMeans(perm_results[[p]]$lW_train_all)
  cc_train_null[p, ] <- colMeans(perm_results[[p]]$cc_train_all)
}
#2. Then: Compute observed means
cc_test_obs_mean <- colMeans(res_obs_new$cc_test_all)
lW_train_obs_mean <- colMeans(res_obs_new$lW_train_all)
cc_train_obs_mean <- colMeans(res_obs_new$cc_train_all)
#3. Now: Calculate empirical p-values
emp_p_cc <- sapply(1:K, function(k) {
    (sum(cc_test_null[, k] >= cc_test_obs_mean[k]) + 1) / (n_perm + 1)
  })

emp_p_lW <- sapply(1:K, function(k) {
  (sum(lW_train_null[, k] >= lW_train_obs_mean[k]) + 1) / (n_perm + 1)
})

emp_p_train_cc <- sapply(1:K, function(k) {
  (sum(cc_train_null[, k] >= cc_train_obs_mean[k]) + 1) / (n_perm + 1)
})
#Optionally, you can also apply cummax() to enforce monotonicity across components:
emp_p_cc_cummax <- cummax(emp_p_cc)
emp_p_lW_cummax <- cummax(emp_p_lW)
emp_p_cc_train_cummax <- cummax(emp_p_train_cc)
#✅ Optional: Print results
round(data.frame(
  Component = 1:K,
  cc_obs = cc_test_obs_mean,
  p_cc = emp_p_cc,
  p_cc_monotonic = emp_p_cc_cummax,
  lW_obs = lW_train_obs_mean,
  p_lW = emp_p_lW,
  p_lW_monotonic = emp_p_lW_cummax,
  cc_train_obs = cc_train_obs_mean,
  p_cc_train = emp_p_train_cc,
  p_cc_train_monotonic = emp_p_cc_train_cummax
), 4)


# do some diagnostics of your function:
cc_df <- data.frame(
  Component = factor(rep(1:K, each = n_perm)),
  Permuted = as.vector(cc_test_null)
)

cc_obs_df <- data.frame(
  Component = factor(1:K),
  Observed = cc_test_obs_mean
)

ggplot(cc_df, aes(x = Component, y = Permuted)) +
  geom_violin(fill = "lightgray") +
  geom_point(data = cc_obs_df, aes(x = Component, y = Observed), color = "red", size = 2) +
  labs(title = "Observed vs Permuted Canonical Correlations", y = "Canonical Correlation")

cc_train_df <- data.frame(
  Component = factor(rep(1:K, each = n_perm)),
  Permuted = as.vector(cc_train_null)
)

cc_obs_train_df <- data.frame(
  Component = factor(1:K),
  Observed = cc_train_obs_mean
)

ggplot(cc_train_df, aes(x = Component, y = Permuted)) +
  geom_violin(fill = "lightgray") +
  geom_point(data = cc_obs_train_df, aes(x = Component, y = Observed), color = "red", size = 2) +
  labs(title = "Observed vs Permuted Canonical Correlations", y = "Canonical Correlation")

lW_train_df <- data.frame(
  Component = factor(rep(1:K, each = n_perm)),
  Permuted = as.vector(lW_train_null)
)

lW_obs_train_df <- data.frame(
  Component = factor(1:K),
  Observed = lW_train_obs_mean
)

ggplot(lW_train_df, aes(x = Component, y = Permuted)) +
  geom_violin(fill = "lightgray") +
  geom_point(data = lW_obs_train_df, aes(x = Component, y = Observed), color = "red", size = 2) +
  labs(title = "Observed vs Permuted Canonical Correlations", y = "Canonical Correlation")


df_cc_test_obs_mean <- cc_test_obs_mean %>%
  as.data.frame()%>%
  pivot_longer(everything(), names_to = 'components', values_to = 'real_values')%>%
  rownames_to_column()%>%
  mutate(components = paste0('V', rowname, sep = ''))%>%
  dplyr::select(components, real_values)

cc_test_null %>%
  as.data.frame()%>%
  pivot_longer(everything(), names_to = 'components', values_to = 'values')%>%
  left_join(., df_cc_test_obs_mean, by = 'components')%>%
  mutate(components = factor(components, levels = paste0("V", 1:K)))%>%
  ggplot(aes(x = values))+geom_histogram()+
  facet_wrap(.~ components)+
  geom_vline(aes(xintercept = real_values), color = 'red')


df_cc_train_obs_mean <- cc_train_obs_mean %>%
  as.data.frame()%>%
  pivot_longer(everything(), names_to = 'components', values_to = 'real_values')%>%
  rownames_to_column()%>%
  mutate(components = paste0('V', rowname, sep = ''))%>%
  dplyr::select(components, real_values)

cc_train_null %>%
  as.data.frame()%>%
  pivot_longer(everything(), names_to = 'components', values_to = 'values')%>%
  left_join(., df_cc_train_obs_mean, by = 'components')%>%
  mutate(components = factor(components, levels = paste0("V", 1:K)))%>%
  ggplot(aes(x = values))+geom_histogram()+
  facet_wrap(.~ components, scales = 'free')+
  geom_vline(aes(xintercept = real_values), color = 'red')


df_lW_train_obs_mean <- lW_train_obs_mean %>%
  as.data.frame()%>%
  pivot_longer(everything(), names_to = 'components', values_to = 'real_values')%>%
  rownames_to_column()%>%
  mutate(components = paste0('V', rowname, sep = ''))%>%
  dplyr::select(components, real_values)

lW_train_null %>%
  as.data.frame()%>%
  pivot_longer(everything(), names_to = 'components', values_to = 'values')%>%
  left_join(., df_lW_train_obs_mean, by = 'components')%>%
  mutate(components = factor(components, levels = paste0("V", 1:K)))%>%
  ggplot(aes(x = values))+geom_histogram()+
  facet_wrap(.~ components, scales = 'free')+
  geom_vline(aes(xintercept = real_values), color = 'red')



# test only the first two components that were significant in the full sample
# Dimensions
#K <- ncol(res_obs_new$cc_test_all)
K <- 3

# Initialize matrices
cc_test_null <- matrix(NA, nrow = n_perm, ncol = K)
lW_train_null <- matrix(NA, nrow = n_perm, ncol = K)
cc_train_null <- matrix(NA, nrow = n_perm, ncol = K)

# Fill matrices with means from each permutation
for (p in 1:n_perm) {
  cc_test_null[p, ] <- colMeans(perm_results[[p]]$cc_test_all[,1:K])
  lW_train_null[p, ] <- colMeans(perm_results[[p]]$lW_train_all[,1:K])
  cc_train_null[p, ] <- colMeans(perm_results[[p]]$cc_train_all[,1:K])
}
#2. Then: Compute observed means
cc_test_obs_mean <- colMeans(res_obs_new$cc_test_all[,1:K])
lW_train_obs_mean <- colMeans(res_obs_new$lW_train_all[,1:K])
cc_train_obs_mean <- colMeans(res_obs_new$cc_train_all[,1:K])
#3. Now: Calculate empirical p-values
emp_p_cc <- sapply(1:K, function(k) {
  (sum(cc_test_null[, k] >= cc_test_obs_mean[k]) + 1) / (n_perm + 1)
})

emp_p_lW <- sapply(1:K, function(k) {
  (sum(lW_train_null[, k] >= lW_train_obs_mean[k]) + 1) / (n_perm + 1)
})

emp_p_train_cc <- sapply(1:K, function(k) {
  (sum(cc_train_null[, k] >= cc_train_obs_mean[k]) + 1) / (n_perm + 1)
})
#Optionally, you can also apply cummax() to enforce monotonicity across components:
emp_p_cc_cummax <- cummax(emp_p_cc)
emp_p_lW_cummax <- cummax(emp_p_lW)
emp_p_cc_train_cummax <- cummax(emp_p_train_cc)
#✅ Optional: Print results
round(data.frame(
  Component = 1:K,
  cc_obs = cc_test_obs_mean,
  p_cc = emp_p_cc,
  p_cc_monotonic = emp_p_cc_cummax,
  lW_obs = lW_train_obs_mean,
  p_lW = emp_p_lW,
  p_lW_monotonic = emp_p_lW_cummax,
  cc_train_obs = cc_train_obs_mean,
  p_cc_train = emp_p_train_cc,
  p_cc_train_monotonic = emp_p_cc_train_cummax
), 4)

