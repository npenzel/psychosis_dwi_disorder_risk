# bootstrap-function, echo = FALSE, warning = FALSE, message=FALSE}
bootstrap_winkler <- function(perm_y, perm_x, nf2keep, nIter, critical.value){
  # For the bootstrap function I combined the bootstrap appraoch from:
  # https://github.com/HerveAbdi/data4PCCAR/blob/master/R/Boot4CCA.R
  # but used the permutation test for significant variables from 
  # Winkler, A. M., Renaud, O., Smith, S. M. & Nichols, T. E. 
  # Permutation inference for canonical correlation analysis. NeuroImage 220, 117065 (2020)
  # There MATLAB code was published: 
  # https://github.com/andersonwinkler/PermCCA/blob/master/permcca.m
  .boot.ratio.test <- function(boot.cube, critical.value = 2) {
    boot.cube.mean <- apply(boot.cube, c(1, 2), mean)
    boot.cube.mean_repeat <- array(boot.cube.mean, dim = c(dim(boot.cube)))
    boot.cube.dev <- (boot.cube - boot.cube.mean_repeat)^2
    s.boot <- (apply(boot.cube.dev, c(1, 2), mean))^(1/2)
    boot.ratios <- boot.cube.mean/s.boot
    significant.boot.ratios <- (abs(boot.ratios) > critical.value)
    rownames(boot.ratios) <- rownames(boot.cube)
    rownames(significant.boot.ratios) <- rownames(boot.cube)
    return(list(sig.boot.ratios = significant.boot.ratios, 
                boot.ratios = boot.ratios))
  }
  
  Z = NULL
  Ny = nrow(perm_y)
  Nx = nrow(perm_x)
  if(Ny != Nx){
    print('Y and X do not have the same number of rows.')
  }
  N = Ny 
  rm(Ny, Nx)
  I = diag(N)
  
  # since we do not residualize
  Qz = I
  # when numbers are very small ~ e-15 the sign is different in R
  perm_y  = center(perm_y)
  perm_y <- t(Qz) %*% perm_y # this does not change anything and is related to the residuals i think.
  
  Ny = nrow(perm_y)
  # since I do not write this with the possibility to correct for nuisance
  # variates I specify that R = 0
  R = 0
  W = Z
  Qw = Qz
  perm_x = center(perm_x)
  perm_x <- t(Qz) %*% perm_x
  Nx = nrow(perm_x)
  # since I do not write this with the possibility to correct for nuisance
  # variates I specify that S = 0
  S  = 0
  # select the columns:
  nI = NCOL(perm_y)
  nJ = NCOL(perm_x)
  nN = NROW(perm_y)
  
  # Initial CCA
  perm_y = Qz %*% perm_y
  perm_x = Qw %*% perm_x
  R = R
  S = S
  ccresults = cca(Qz %*% perm_y, Qw %*% perm_x, R, S)
  A = ccresults$A
  B = ccresults$B
  r = ccresults$cc
  
  u = A
  v = B
  
  Lx <- perm_y %*% u
  Ly <- perm_x %*% v 
  Fi <- u * matrix(r, nI, nf2keep, byrow = TRUE)
  Fj <- v * matrix(r, nJ, nf2keep, byrow = TRUE)
  
  fj.boot <- array(NA, dim = c(nJ, nf2keep, nIter))
  dimnames(fj.boot)[1] <- list(colnames(perm_x))
  dimnames(fj.boot)[2] <- list(paste0("Dimension ", 1:nf2keep))
  dimnames(fj.boot)[3] <- list(paste0("Iteration ", 1:nIter))
  fi.boot <- array(NA, dim = c(nI, nf2keep, nIter))
  dimnames(fi.boot)[1] <- list(colnames(perm_y))
  dimnames(fi.boot)[2] <- list(paste0("Dimension ", 1:nf2keep))
  dimnames(fi.boot)[3] <- list(paste0("Iteration ", 1:nIter))
  fi.boot.weight <- array(NA, dim = c(nI, nf2keep, nIter))
  dimnames(fi.boot.weight)[1] <- list(colnames(perm_y))
  dimnames(fi.boot.weight)[2] <- list(paste0("Dimension ", 1:nf2keep))
  dimnames(fi.boot.weight)[3] <- list(paste0("Iteration ", 1:nIter))
  fj.boot.weight <- array(NA, dim = c(nJ, nf2keep, nIter))
  dimnames(fj.boot.weight)[1] <- list(colnames(perm_x))
  dimnames(fj.boot.weight)[2] <- list(paste0("Dimension ", 1:nf2keep))
  dimnames(fj.boot.weight)[3] <- list(paste0("Iteration ", 1:nIter))
  for (ell in 1:nIter) {
    boot.index <- sample(nN, replace = TRUE)
    fi.boot[, , ell] <- t(perm_y[boot.index, ]) %*% Ly[boot.index,]
    fj.boot[, , ell] <- t(perm_x[boot.index, ]) %*% Lx[boot.index,]
    fi.boot.weight[, , ell] <- cor(perm_y[boot.index,],Lx[boot.index,],use="complete.obs")
    fj.boot.weight[, , ell] <- cor(perm_x[boot.index,],Ly[boot.index,],use="complete.obs")
  }
  
  BR.j <- .boot.ratio.test(fj.boot, critical.value)
  BR.i <- .boot.ratio.test(fi.boot, critical.value)
  return.list <- structure(list(bootstrapBrick.i = fi.boot, 
                                bootRatios.i = BR.i$boot.ratios, 
                                bootRatiosSignificant.i = BR.i$sig.boot.ratios, 
                                bootstrapBrick.j = fj.boot, 
                                bootRatios.j = BR.j$boot.ratios,
                                bootWeights.i = fi.boot.weight, 
                                bootWeights.j = fj.boot.weight,
                                bootRatiosSignificant.j = BR.j$sig.boot.ratios, 
                                clin_proj = Lx,
                                dwi_proj = Ly), 
                           class = "bootBrick.ij4plsc")
  return(return.list)
}
