#' Sample G = mu + UV^T + XB^T + e
#'
#' @export
sample.normal.model <- function(n, L, K, sigma, prop.outlier, c, mean.B = 0.0, sd.mu= 1.0, mean.mu = 0.5) {

  outlier = ((1-prop.outlier)*L+1):L
  nb.outlier = length(outlier)
  c = 0.6 # correlation between U1 and X


  U = MASS::mvrnorm(n, mu = rep(0.0,K), Sigma = 1.0 * diag(K))
  V = MASS::mvrnorm(L, mu = rep(0.0,K), Sigma = 1.0 * diag(K))
  mu = matrix(rep(rnorm(L, mean.mu, sd.mu),n), nrow = n,ncol = L,byrow=TRUE)
  X = matrix(BioCompToolsR::sample_correlated_X(U[,1],c), nrow = n,ncol = 1)
  B = rep(0,L)
  B[outlier] = rnorm(nb.outlier, mean.B,1)
  epsilon = MASS::mvrnorm(n, mu = rep(0.0,L), Sigma = sigma * diag(L))

  G = mu + U %*% t(V) + X %*% t(B) + epsilon


  return(list(G = G, X = X, U = U, V = V, B = B, epsilon = epsilon, outlier = outlier))
}

#' Sample P(G) = logistic(mu + UV^T + XB^T + e)
#'
#' @export
sample.binary.model <- function(n, L, K, sigma, prop.outlier, c, ploidy) {

  logistic = function(x){
    return(1/(1+exp(-x)))
  }

  normal.model = sample.normal.model(n, L, K, sigma, prop.outlier, c, mean.mu = 0.0)
  P_G = logistic(normal.model$G)
  G = apply(P_G, 1:2, function(p) rbinom(1,ploidy,p))

  normal.model$G = G
  return(normal.model)
}


#' Sample Phenotype = Somme(B_j * G_j) + epsilon
#' Same simulation that in "Genomic inflation factors under polygenic inheritance"
#'
#' @export
sample.phenotype <- function(G, squared.h, number.causal) {

  assertthat::assert_that(squared.h != 0.0)

  n = nrow(G)
  L = ncol(G)

  assertthat::assert_that(number.causal < L)

  outlier = sample(1:L,number.causal)
  b = matrix(0,L,1)

  b[outlier] = rnorm(number.causal)

  g = G %*% b
  e.sigma2 = var(g) * (1 - squared.h) / squared.h
  X = g + matrix(rnorm(n, sqrt(e.sigma2)), n, 1)
  return(list(X = X, outlier = outlier, B = b ))
}
