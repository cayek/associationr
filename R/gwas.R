#' X = G * B + epsilon + regularisation ridge
#'
#' @export
gwas.linearmodel.ridge <- function(G,X,parameter) {
  if(is.null(parameter$lambda)) {
    stop("Parameter lambda is mandatory")
  }
  lambda = parameter$lambda
  n = nrow(G)
  L = ncol(G)

  B = solve(crossprod(G,G) + lambda * diag(nrow = L, ncol = L), crossprod(G, X))

  mad = median(abs(as.numeric(B)))
  sd = 1.4826 * mad

  zscore = (as.numeric(B) - 0.0) / sd

  # pvalue

  pvalue = 2 * pnorm(abs(zscore),lower.tail = FALSE)

  return(list(pvalue = pvalue, B = B, zscore = zscore))

}


#' X = mu + G_j * B + epsilon
#'
#' @export
gwas.linearmodel <- function(G,X,parameter) {

n = nrow(G)
L = ncol(G)
B = 1:L
zscore = 1:L
pvalue = 1:L
for(j in 1:L) {
  l=lm(X~G[,j] )
  zscore[j] = summary(l)$coefficients[2,3]
  pvalue[j] = summary(l)$coefficients[2,4]
  B[j] = summary(l)$coefficients[2,1]
}
return(list(pvalue = pvalue, B = B, zscore = zscore))
}
