# ---------------------------------------------------------------------------- #
#                     Create CV folds for cross-validation.                    #
# ---------------------------------------------------------------------------- #

# ------------------------------------------------------------------
# 1.  Create one reproducible fold structure to re-use everywhere
# ------------------------------------------------------------------
create_cv_folds <- function(n_obs,
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
