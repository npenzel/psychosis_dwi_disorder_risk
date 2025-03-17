# CCA function
cca <- function(Y, X, R, S) {
  # Compute CCA
  # Inputs are:
  # Y =
  # X = 
  # R = 
  # S =
  N <- nrow(Y)
  # Compute QR decomposition
  qr_decomp <- qr(Y, LAPACK = TRUE)
  Qy <- qr.Q(qr_decomp)
  Ry <- qr.R(qr_decomp)
  iY <- qr_decomp$pivot
  # same for imaging X variables:
  #[Qx,Rx,iX] = qr(X,0);
  # Compute QR decomposition
  qr_decomp <- qr(X, LAPACK = TRUE)
  Qx <- qr.Q(qr_decomp)
  Rx <- qr.R(qr_decomp)
  iX <- qr_decomp$pivot
  # Calculate rank of matrices Y and X
  #   K  = min(rank(Y),rank(X));
  K <- min(qr(Y)$rank, qr(X)$rank)
  svd_result <- svd(t(Qy) %*% Qx)
  D <- svd_result$d
  L <- svd_result$u
  M <- svd_result$v
  cc <- pmin(pmax(D[1:K], 0), 1)
  A <- solve(Ry) %*% L[, 1:K] * sqrt(N - R)
  B <- solve(Rx) %*% M[, 1:K] * sqrt(N - S)
  A[iY, ] <- A
  B[iX, ] <- B
  return(list(A = A, B = B, cc = cc))
}