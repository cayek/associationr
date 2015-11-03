#' LFMM with ridge likelihood
#'
#' @export
lfmm.ml.ridge <- function(G, X, K, lambda) {

  mu = apply(G, 2, mean)
  G_ = G - matrix(mu,nrow(G),ncol(G),byrow = TRUE)
  P_X_D_t = diag(nrow(G)) - tcrossprod(X,X) * 1 /(as.numeric(crossprod(X,X)) + lambda)
  eigen.P_X_D_t = eigen(P_X_D_t, symmetric = TRUE)
  sqrt.P_X_D_t = eigen.P_X_D_t$vectors %*% diag(sqrt(eigen.P_X_D_t$values)) %*% t(eigen.P_X_D_t$vectors)
  if(K!=0) {
    # L = sqrt.P_X_D_t^-1 svd_K( sqrt.P_X_D_t G )
    svd.G = svd(sqrt.P_X_D_t %*% G_, nu = K,nv=K)
    PU = svd.G$u %*% diag(svd.G$d[1:K])
    V = svd.G$v
    PL = PU %*% t(V)
    L = solve(sqrt.P_X_D_t,PL)
    U = solve(sqrt.P_X_D_t,PU)
  }else{
    U = NULL
    V = NULL
    L = matrix(0,nrow(G),ncol(G))
  }

  # B = (X^T X + lambda)^-1 X^T ( G - L )
  B = 1 /(as.numeric(crossprod(X,X)) + lambda) * crossprod(X,G_-L)

  ## robust estimation of standard deviation
  mad = median(abs(as.numeric(B)))
  sd = 1.4826 * mad

  zscore = (as.numeric(B) - 0.0) / sd

  # pvalue

  pvalue = 2 * pnorm(abs(zscore),lower.tail = FALSE)

  return(list(pvalue = pvalue, U = U, V = V, mu = mu, B = B, zscore = zscore))


}


#' LFMM with ridge likelihood implementation 2
#'
#' @export
lfmm.ml.ridge2 <- function(G, X, K, lambda) {

  mu = apply(G, 2, mean)
  G_ = G - matrix(mu,nrow(G),ncol(G),byrow = TRUE)
  P_X_D_t = diag(nrow(G)) - tcrossprod(X,X) * 1 /(as.numeric(crossprod(X,X)) + lambda)
  eigen.P_X_D_t = eigen(P_X_D_t, symmetric = TRUE)
  sqrt.P_X_D_t = eigen.P_X_D_t$vectors %*% diag(sqrt(eigen.P_X_D_t$values)) %*% t(eigen.P_X_D_t$vectors)
  if(K!=0) {
    # L = sqrt.P_X_D_t^-1 svd_K( sqrt.P_X_D_t G )
    svd.G = svd(sqrt.P_X_D_t %*% G_, nu = K,nv=K)
    PU = svd.G$u %*% diag(svd.G$d[1:K])
    V = svd.G$v
    PL = PU %*% t(V)
    #L = solve(P_X_D_t,L)
    #U = solve(P_X_D_t,PU)
    U = NULL
  }else{
    U = NULL
    V = NULL
    L = matrix(0,nrow(G),ncol(G))
  }

  # B = (X^T X + lambda)^-1 X^T ( G - L )
  X.norm.square = as.numeric(crossprod(X,X))
  # cat("1/lamba*X^T PL =",sum(1 / lambda * crossprod(X,PL)))
  B = 1 /(X.norm.square + lambda) * crossprod(X,G_) - 1 / sqrt((X.norm.square + lambda)*lambda ) * crossprod(X,PL)

  ## robust estimation of standard deviation
  mad = median(abs(as.numeric(B)))
  sd = 1.4826 * mad

  zscore = (as.numeric(B) - 0.0) / sd

  # pvalue

  pvalue = 2 * pnorm(abs(zscore),lower.tail = FALSE)

  return(list(pvalue = pvalue, U = U, V = V, mu = mu, B = B, zscore = zscore))


}

