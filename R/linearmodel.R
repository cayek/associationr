#' only LM
#'
#' @export
linear_model <- function( G, X ) {

  aux = function(g) {
    l=lm(g~X )
    return(data.frame(tscore = summary(l)$coefficients[2,3],
                      pvalue = summary(l)$coefficients[2,4],
                      B = summary(l)$coefficients[2,1]))
  }
  list = apply(G, 2, aux)
  lm.fit = Reduce( rbind, list)

  return(lm.fit)

}
