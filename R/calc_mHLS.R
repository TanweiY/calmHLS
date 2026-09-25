#' Calculate overall mHLS and domain-specific sub-scores
#'
#' Combines the pre-calculated individual scores and diet CpGs with the final
#' model coefficients to produce the overall mHLS  and raw
#' weighted sums for each lifestyle domain. Missing predictors are imputed using
#' training-set medians.
#'
#' @param scores_and_diet A data.frame as returned by \code{calc_behavior_scores()}.
#'   Must contain \code{ID} and all predictor variables (or a subset thereof).
#' @param domains A named list of character vectors defining the predictors for
#'   each domain. Defaults to internal domain grouping.
#' @param coeffs A data.frame with columns \code{Var} and \code{s1}, the
#'   final model coefficients. Defaults to the internal \code{mHLS_coefficients}.
#' @param medians A data.frame with columns \code{Variable} and \code{Median},
#'   used to impute missing predictors. Defaults to internal \code{train_medians}.
#'
#'
#' @return A data.frame containing:
#'  \item{ID}{Sample identifier.}
#'  \item{mHLS}{Overall DNA methylation-based healthy lifestyle score (linear predictor). Higher values indicate a healthier profile.}
#'  \item{RawScore_noSmoking}{Raw domain score for smoking behavior (higher = healthier, i.e., less smoking exposure).}
#'  \item{RawScore_lowBMI}{Raw domain score for BMI (higher = healthier BMI profile).}
#'  \item{RawScore_lowAlcohol}{Raw domain score for alcohol consumption (higher = healthier, i.e., lower alcohol intake).}
#'  \item{RawScore_goodDiet}{Raw domain score for dietary quality (higher = healthier diet).}
#'  \item{mHLS_category}{Training-set tertile-based category of overall mHLS: **Low**, **Medium**, **High** (higher is healthier).}
#'  \item{noSmoking_category}{Training-set tertile-based category for the smoking domain (higher indicates healthier behavior).}
#'  \item{lowBMI_category}{Training-set tertile-based category for the BMI domain (higher indicates healthier BMI).}
#'  \item{lowAlcohol_category}{Training-set tertile-based category for the alcohol domain (higher indicates lower alcohol intake).}
#'  \item{goodDiet_category}{Training-set tertile-based category for the diet domain (higher indicates healthier diet).}
#'
#'
#' @export
#' @examples
#' \dontrun{
#'   final <- calc_mHLS(scores_and_diet)
#'   head(final)
#' }
#'
calc_mHLS <- function(scores_and_diet) {

  coeffs<-mHLS_coefficients
  medians<-train_medians

  # define default domain groups
  smoking_vars <- c("mSmk_233", "probs_fsmk")
  bmi_vars     <- "mBMI_397"
  alcohol_vars <- "mAlc_450"

  all_predictors <- setdiff(coeffs$Var, "(Intercept)")

  diet_vars <- setdiff(all_predictors, c(smoking_vars, bmi_vars, alcohol_vars))
  domains <- list(noSmoking = smoking_vars, lowBMI = bmi_vars,
                  lowAlcohol = alcohol_vars, goodDiet = diet_vars)

  # extract ID
  ids <- scores_and_diet$ID
  # build calculation matrix from available predictors
  all_target_vars <- setdiff(coeffs$Var, "(Intercept)")
  avail_vars <- intersect(all_target_vars, colnames(scores_and_diet))
  missing_vars <- setdiff(all_target_vars, avail_vars)
  message("Calculating mHLS: ", length(missing_vars),
          " missing features will be imputed using the median values from the training set.")

  calc_mat <- as.matrix(scores_and_diet[, avail_vars, drop = FALSE])

  # impute missing variables with training medians
  median_lookup <- setNames(medians$Median, medians$Variable)
  if (length(missing_vars) > 0) {
    for (mv in missing_vars) {
      calc_mat <- cbind(calc_mat, rep(median_lookup[mv], nrow(calc_mat)))
      colnames(calc_mat)[ncol(calc_mat)] <- mv
    }
  }

  # ensure column order matches weights
  calc_mat <- calc_mat[, all_target_vars, drop = FALSE]

  # weights (without intercept)
  w_all <- coeffs$s1[match(all_target_vars, coeffs$Var)]
  intercept <- coeffs$s1[coeffs$Var == "(Intercept)"]

  # overall mHLS
  mHLS <- as.vector(calc_mat %*% w_all + intercept)

  # domain sub-scores (raw weighted sum, no intercept)
  df_allscores <- data.frame(ID = ids, mHLS = mHLS, stringsAsFactors = FALSE)
  for (dom in names(domains)) {
    vars <- domains[[dom]]
    # keep only those that exist in calc_mat (they will after imputation)
    vars_present <- intersect(vars, colnames(calc_mat))
    if (length(vars_present) == 0) next
    w_dom <- coeffs$s1[match(vars_present, coeffs$Var)]
    raw_score <- as.vector(calc_mat[, vars_present, drop = FALSE] %*% w_dom)
    df_allscores[[paste0("RawScore_", dom)]] <- raw_score
  }

  # mHLS category
  df_allscores$mHLS_category <- cut(df_allscores$mHLS,
    breaks = mHLS_tertile_cutoffs,
    labels = c("Low", "Medium", "High"),
    right = FALSE)

  # Domain categories
  for (dom in names(domain_tertile_cuts)) {
    raw_col <- paste0("RawScore_", dom)
    if (raw_col %in% colnames(df_allscores)) {
      cuts <- domain_tertile_cuts[[dom]]
      df_allscores[[paste0(dom, "_category")]] <- cut(
        df_allscores[[raw_col]],
        breaks = c(-Inf, cuts[1], cuts[2], Inf),
        labels = c("Low", "Medium", "High"),
        right = FALSE
      )
    }
  }


  return(df_allscores)
}





