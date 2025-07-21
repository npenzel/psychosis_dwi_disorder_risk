# ---------------------------------------------------------------------------- #
#                     Run the cross-validated CCA analysis                     #
# ---------------------------------------------------------------------------- #

source('cv_function/center_cv_function.R')
source('cv_function/impute_cv_function.R')
source('cv_function/align_and_singflip_function.R')

# -------------------------------------------------------------
# 2)  CV-CCA with alignment + your centre_cv() inside
# -------------------------------------------------------------
run_cv_cca <- function(X_mat,
                       Y_mat,
                       folds_list,        
                       ref_X_load,
                       ref_Y_load,
                       K = NULL,
                       permute_Y = FALSE) {
  
  n_repeats <- length(folds_list)
  n_folds   <- length(folds_list[[1]])
  if (is.null(K)){
    K         <- min(ncol(ref_X_load), ncol(ref_Y_load))
  }
  
  cc_test_all <- matrix(NA, nrow = n_repeats * n_folds, ncol = K)
  cc_train_all <- matrix(NA, nrow = n_repeats * n_folds, ncol = K)
  fold_index  <- 0
  alignment_history <- matrix(0, nrow = n_repeats * n_folds, ncol = K)
  signflip_matrix <- matrix(FALSE, nrow = n_repeats * n_folds, ncol = K)
  n_signflips <- integer(n_repeats * n_folds)
  
  for (rep_idx in seq_len(n_repeats)) {
    folds <- folds_list[[rep_idx]]
    
    for (fold_idx in seq_len(n_folds)) {
      fold_index <- fold_index + 1
      
      test_id  <- folds[[fold_idx]]
      train_id <- setdiff(seq_len(nrow(X_mat)), test_id)
      
      # ---------- Imputation ----------
      X_train_raw <- impute_cv(X_mat[train_id, ], X_mat[test_id, ],
                               'studygroup_num', 'uhr_1stdegree_re', 'train')
      X_test_raw  <- impute_cv(X_mat[train_id, ], X_mat[test_id, ],
                               'studygroup_num', 'uhr_1stdegree_re', 'test')
      Y_train_raw <- Y_mat[train_id, ]
      Y_test_raw  <- Y_mat[test_id , ]
      
      if (permute_Y) {
        Y_train_raw <- Y_train_raw[sample(nrow(Y_train_raw)), ]
      }
      
      # ---------- Centre using your centre_cv() ----------
      X_train <- center_cv(X_train_raw)                       # centres & trims const-cols
      X_test  <- center_cv(X_train_raw, X_test_raw, 'apply')  # same means/cols
      Y_train <- center_cv(Y_train_raw)
      Y_test  <- center_cv(Y_train_raw, Y_test_raw, 'apply')
      
      # ---------- CCA on training ----------
      train_result <- permcca_winkler(Y_train, X_train, 1, TRUE)  # No permutation here
      A <- train_result$cwr  # Canonical weights (Y side)
      B <- train_result$cwl  # Canonical weights (X side)
      
      cc_train_all[fold_index, ] <- train_result$cc[1:K]
      
      aligned_check <- cca_check_align_and_signflip(X_fold = train_result$cwr[, 1:K, drop=FALSE],
                                                    Y_fold = train_result$cwl[, 1:K, drop=FALSE],
                                                    X_ref = ref_X_load[, 1:K, drop=FALSE],
                                                    Y_ref = ref_Y_load[, 1:K, drop=FALSE])
      
      reorder    <- aligned_check$alignment_indices    # e.g. reorder[1] = 5 means fold-col5 → ref-col1
      alignment_history[fold_index, ] <- reorder
      
      # Optionally: store sign flip diagnostics
      signflip_matrix[fold_index, ] <- aligned_check$flip_flags
      n_signflips[fold_index] <- aligned_check$n_flips
      
      X_test_can <- X_test %*% A
      Y_test_can <- Y_test %*% B
      
      X_train_can <- X_train %*% A
      Y_train_can <- Y_train %*% B
      cc_vec <- diag(cor(X_test_can[,1:K], Y_test_can[,1:K]))
      cc_test_all[fold_index, ] <- cc_vec
      # Canonical correlations (aligned)
      cc_train_all[fold_index, ] <- diag(cor(X_train_can[,1:K], Y_train_can[,1:K]))
    }
  }
  return(list(cc_test_all = cc_test_all,
              cc_train_all = cc_train_all,
              alignment_history = alignment_history,
              signflip_matrix = signflip_matrix,
              n_signflips = n_signflips))
}