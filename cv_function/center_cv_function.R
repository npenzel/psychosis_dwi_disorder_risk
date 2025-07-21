# ---------------------------------------------------------------------------- #
#                     Center data within cross-validation                      #
# ---------------------------------------------------------------------------- #

# Function to center a data matrix (X) and apply the same centering transformation to a test set (X_test_cv_perm)
# The centering is performed by subtracting the column means of the training data (X) and adjusting for columns that have no variation.
center_cv <- function(X, X_test_cv_perm = NULL, apply_mean = NULL) {
  # Check if apply_mean exists (indicating whether we should apply centering to the test set as well)
  # X = training data to center directly and learn,
  # X_test_cv_perm = test data to apply centering from training data.
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