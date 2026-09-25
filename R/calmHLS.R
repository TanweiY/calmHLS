#' Full mHLS pipeline
#'
#' A convenience function that runs \code{calc_behavior_scores()} followed by
#' \code{calc_mHLS()}.
#'
#' @inheritParams calc_behavior_scores
#' @param ... Additional arguments passed to \code{calc_mHLS}, e.g. custom domains.
#'
#' @return The output of \code{calc_mHLS}, a data.frame with mHLS and domain scores.
#' @export
#' @examples
#' \dontrun{
#'   final <- calmHLS(beta_lifedf, base_SSt)
#' }

calmHLS <- function(beta, samplesheet) {
  step1 <- calc_behavior_scores(beta, samplesheet)
  step2 <- calc_mHLS(step1)
  return(step2)
}




