#' Sampling of linear model parameters using bootstrap method.
#'
#' \code{boot.lm} returns the function which sample linear model parameters.
boot.lm <- function(formula, data, ...) {
  function(nbsample = 1){

    estimator <- function() {
      stats::lm(formula=formula,
         data=data[sample(nrow(data), replace=TRUE),], ...)$coef
    }

    bstrap <- sapply(X=1:nbsample,
                     FUN=function(x) estimator())

    return(bstrap)
  }
}


