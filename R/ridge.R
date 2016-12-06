#' LFMM with ridge likelihood
#'
#' @export
lfmm.ml.ridge <- function(G, X, K, lambda) {

  mu = apply(G, 2, mean)
  G_ = G - matrix(mu,nrow(G),ncol(G),byrow = TRUE)
  P_X_D_t = diag(nrow(G)) - tcrossprod(X,X) * 1 /(as.numeric(crossprod(X,X)) + lambda)
  eigen.P_X_D_t = eigen(P_X_D_t, symmetric = TRUE)
  sqrt.P_X_D_t = eigen.P_X_D_t$vectors %*% diag(sqrt(eigen.P_X_D_t$values),nrow(G),nrow(G)) %*% t(eigen.P_X_D_t$vectors)
  if(K!=0) {
    # L = sqrt.P_X_D_t^-1 svd_K( sqrt.P_X_D_t G )
    svd.G = svd(sqrt.P_X_D_t %*% G_, nu = K,nv=K)
    PU = svd.G$u %*% diag(svd.G$d[1:K],K,K)
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
  if(K!=0) {
    P_X_D_t = diag(nrow(G)) - tcrossprod(X,X) * 1 /(as.numeric(crossprod(X,X)) + lambda)
    eigen.P_X_D_t = eigen(P_X_D_t, symmetric = TRUE)
    sqrt.P_X_D_t = eigen.P_X_D_t$vectors %*% diag(sqrt(eigen.P_X_D_t$values),nrow(G),nrow(G)) %*% t(eigen.P_X_D_t$vectors)
    # L = sqrt.P_X_D_t^-1 svd_K( sqrt.P_X_D_t G )
    svd.G = svd(sqrt.P_X_D_t %*% G_, nu = K,nv=K)
    PU = svd.G$u %*% diag(svd.G$d[1:K],K,K)
    V = svd.G$v
    PL = PU %*% t(V)
    #L = solve(P_X_D_t,L)
    #U = solve(P_X_D_t,PU)
    U = NULL
  }else{
    U = NULL
    V = NULL
    L = matrix(0,nrow(G),ncol(G))
    PL = matrix(0,nrow(G),ncol(G))
  }

  # B = (X^T X + lambda)^-1 X^T ( G - L )
  X.norm.square = as.numeric(crossprod(X,X))
  # cat("1/lamba*X^T PL =",sum(1 / lambda * crossprod(X,PL)))
  if( lambda==0.0 ){
    B = 1 /(X.norm.square) * crossprod(X,G_)
  } else {
    B = 1 /(X.norm.square + lambda) * crossprod(X,G_) - 1 / sqrt((X.norm.square + lambda)*lambda ) * crossprod(X,PL)
  }
  ## robust estimation of standard deviation
  mad = median(abs(as.numeric(B)))
  sd = 1.4826 * mad

  zscore = (as.numeric(B) - 0.0) / sd

  # pvalue

  pvalue = 2 * pnorm(abs(zscore),lower.tail = FALSE)

  return(list(pvalue = pvalue, U = U, V = V, mu = mu, B = B, zscore = zscore))


}


#' LFMM with ridge likelihood implementation 2 + normalization
#'
#' @export
lfmm.ml.ridge2.normalization <- function(G, X, K, lambda) {

  # P = apply(G, 2, mean)
  # aux = sqrt(P*(1-P))
  aux = apply(G, 2, sd)
  aux[aux==0.0] = 1.0
  G_ = G %*% diag(1/aux)
  return(lfmm.ml.ridge2(G_,X,K, lambda))
}



#' LFMM with ridge likelihood implementation 3
#'
#' @export
lfmm.ml.ridge3 <- function(G, X, K, lambda) {

  mu = apply(G, 2, mean)
  G_ = G - matrix(mu,nrow(G),ncol(G),byrow = TRUE)

  # compute sqrt(P_X_D_t)
  eigen.x = eigen(tcrossprod(X,X) * 1 / (as.numeric(crossprod(X,X))), symmetric = TRUE)
  D = eigen.x$values
  D[abs(D)<1e-10] = 0.0
  D[D != 0.0] = sqrt( lambda / (as.numeric(crossprod(X,X)) + lambda) )
  D[D == 0.0] = 1.0
  sqrt.P_X_D_t = eigen.x$vectors %*% diag(D,nrow(G),nrow(G)) %*% t(eigen.x$vectors)
  if(K!=0) {
    # L = sqrt.P_X_D_t^-1 svd_K( sqrt.P_X_D_t G )
    svd.G = svd(sqrt.P_X_D_t %*% G_, nu = K,nv=K)
    PU = svd.G$u %*% diag(svd.G$d[1:K],K,K)
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




#' LFMM with ridge likelihood implementation 2 + LR for score
#'
#' @export
lfmm.ml.ridge2.lr <- function(G, X, K, lambda) {

  model.fit <- function(tozero.index=NULL) {

    notzero.index = !(1:ncol(G) %in% tozero.index)
    mu = apply(G, 2, mean)
    G_ = G - matrix(mu,nrow(G),ncol(G),byrow = TRUE)
    if(K!=0) {
      P_X_D_t = diag(nrow(G)) - tcrossprod(X,X) * 1 /(as.numeric(crossprod(X,X)) + lambda)
      eigen.P_X_D_t = eigen(P_X_D_t, symmetric = TRUE)
      sqrt.P_X_D_t = eigen.P_X_D_t$vectors %*% diag(sqrt(eigen.P_X_D_t$values),nrow(G),nrow(G)) %*% t(eigen.P_X_D_t$vectors)

      # we project the genotipe matrix
      projected.G_ = G_
      projected.G_[,notzero.index] = sqrt.P_X_D_t %*% G_[,notzero.index]

      svd.G = svd(projected.G_, nu = K,nv=K)
      PU = svd.G$u %*% diag(svd.G$d[1:K],K,K)
      V = svd.G$v
      PL = PU %*% t(V)
      U = NULL
    }else{
      U = NULL
      V = NULL
      L = matrix(0,nrow(G),ncol(G))
      PL = matrix(0,nrow(G),ncol(G))
    }

    # B = (X^T X + lambda)^-1 X^T ( G - L )
    X.norm.square = as.numeric(crossprod(X,X))
    # cat("1/lamba*X^T PL =",sum(1 / lambda * crossprod(X,PL)))
    if( lambda==0.0 ){
      B = 1 /(X.norm.square) * crossprod(X,G_)
    } else {
      B = 1 /(X.norm.square + lambda) * crossprod(X,G_) - 1 / sqrt((X.norm.square + lambda)*lambda ) * crossprod(X,PL)
      B[tozero.index] = 0.0
    }

    # compute sigma
    L = solve(sqrt.P_X_D_t, PL)
    e = norm(G_ - L - X %*% B, type = "F")^2 + lambda * norm(B,type = "F")^2
    sigma2 = e / nrow(G) / ncol(G)


    # compute likelihood, RMK : we not add + ncol(G) / 2 * log(2 * pi * sigma2 / lambda) = + ncol(G) / 2 * log(2 * pi * tho2)
    # RMK = 1 / (2 * sigma2) * e = nrow(G) * ncol(G) / 2 so we do not add it to the likelihood
    minus.log.likelihood = nrow(G) * ncol(G) / 2 * log(2 * pi * sigma2)
    return(list(e = e, minus.log.likelihood = minus.log.likelihood, B = B, U = U, V = V, sigma2 = sigma2))
  }

  # Compute H1 model
  res = model.fit()

  ###############################SCORE###################################
  lr.score = 1:ncol(G)
  for(j in 1:ncol(G)) {
    lr.score[j] = model.fit(j)$minus.log.likelihood
    # lr = - 2 * log(L0/L1)
    lr.score[j] = 2 * lr.score[j] - 2 * res$minus.log.likelihood
  }

  pvalue = pchisq(lr.score, df = 1, lower.tail = FALSE)


  return(list(pvalue = pvalue, U = res$U, V = res$V, mu = res$mu, B = res$B))


}





