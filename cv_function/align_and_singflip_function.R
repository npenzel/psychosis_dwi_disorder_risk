# ---------------------------------------------------------------------------- #
#                     check alignment and sign-flip function                   #
# ---------------------------------------------------------------------------- #

cca_check_align_and_signflip <- function(X_fold, Y_fold, X_ref, Y_ref, verbose = FALSE) {
  K <- ncol(X_ref)
  # Compute similarity matrices (absolute correlation)
  sim_X <- abs(cor(X_fold, X_ref))
  sim_Y <- abs(cor(Y_fold, Y_ref))
  joint_sim <- (sim_X + sim_Y) / 2
  
  # Solve best matching (Hungarian algorithm)
  assignment <- solve_LSAP(joint_sim, maximum = TRUE)
  
  flip_flags <- logical(K)
  flip_details <- data.frame(
    Component = 1:K,
    Assigned = NA_integer_,
    Sign_X = NA_integer_,
    Sign_Y = NA_integer_,
    FlipScore = NA_real_,
    FlipApplied = FALSE
  )
  
  for (k in 1:K) {
    j <- assignment[k]
    
    cor_X <- cor(X_fold[, j], X_ref[, k])
    cor_Y <- cor(Y_fold[, j], Y_ref[, k])
    
    # Determine signs
    sign_X <- sign(cor_X)
    sign_Y <- sign(cor_Y)
    
    # Flip if signs are inconsistent
    flip_applied <- sign_X != sign_Y
    flip_flags[k] <- flip_applied
    
    flip_details[k, ] <- list(
      Component = k,
      Assigned = j,
      Sign_X = sign_X,
      Sign_Y = sign_Y,
      FlipScore = cor_X + cor_Y,  # optional: keep for reference
      FlipApplied = flip_applied
    )
    
    if (verbose) {
      cat(sprintf("Component %d → %d | cor_X = %.3f | cor_Y = %.3f | flip = %s\n",
                  j, k, cor_X, cor_Y, ifelse(flip_applied, "-", "+")))
    }
  }
  
  total_flips <- sum(flip_flags)
  
  return(list(
    alignment_indices = assignment,
    flip_flags = flip_flags,
    n_flips = total_flips,
    flip_details = flip_details
  ))
}