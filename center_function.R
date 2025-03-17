# center function
center <- function(X) {
  icte <- apply(diff(X, 1, 1)^2, 2, sum) == 0
  X <- sweep(X, 2, colMeans(X), "-")
  X <- X[, !icte]
  return(X)
}
