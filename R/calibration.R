#' calibrate zscore with mad estimator of variance
#'
#' @export
calibration.mad <- function(zscore) {

  ## robust estimation of standard deviation
  mad = median(abs(zscore))
  sd = 1.4826 * mad

  calibrated.zscore = (zscore - 0.0) / sd

  # pvalue

  calibrated.pvalue = 2 * pnorm(abs(calibrated.zscore),lower.tail = FALSE)

  return(list(calibrated.pvalue = calibrated.pvalue, calibrated.zscore = calibrated.zscore))

}
