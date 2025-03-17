# Permute CCA based on the method from Winkler et al., XXX
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
  
  # For permutation, I translated the code from Winkler et al., 2020
  # form Matlab to R:
  # There MATLAB code was published: 
  # https://github.com/andersonwinkler/PermCCA/blob/master/permcca.m
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
  ccresults = cca(Qz %*% Y, Qw %*% X, R, S) # before it was cca
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
    # Otherwise, use the permtuation set to shuffle
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
    
    # For each canonical variable
    for (k in 1:K) {
      #cat(k, ' ')
      rperm_all <- cca(Qz %*% U[idxY, k:ncol(U)], Qw %*% V[idxX, k:ncol(V)], R, S)
      rperm <- rperm_all$cc
      lWtmp <- -rev(cumsum(rev(log(1 - rperm^2))))
      lW[k] <- lWtmp[1]
    }
    
    if (p == 1) {
      lW1 <- lW
    }
    
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
