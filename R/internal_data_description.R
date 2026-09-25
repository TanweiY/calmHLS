#’ Internal package data for calmHLS
#’
#’ @name calmHLS-internal
#’ @keywords internal
#’ @description
#’ Internal datasets stored in R/sysdata.rda and used by the package
#’ functions for calculating methylation-based lifestyle scores, smoking
#’ status probabilities, and mHLS categories.
#’
#’ \describe{
#’   \item{all_required_cpgs}{
#’   Character vector containing all CpGs required across the implemented
#’   methylation-based lifestyle scores.
#’   }
#’
#’   \item{diet_cpgs}{
#’   Character vector of CpGs used in the diet domain score.
#’   }
#’
#’   \item{smk233_coef}{
#’   Data frame containing CpG identifiers and corresponding weights for
#’   calculation of the mSmk-233 smoking score.
#’   }
#’
#’   \item{alcohol_coef}{
#’   Data frame containing CpG identifiers and corresponding weights for
#’   calculation of the mAlc-450 alcohol consumption score.
#’   }
#’
#’   \item{bmi_coef}{
#’   Data frame containing CpG identifiers and corresponding weights for
#’   calculation of the mBMI-397 body mass index score.
#’   }
#’
#’   \item{mHLS_coefficients}{
#’   Data frame containing regression coefficients used to calculate the
#’   overall methylation-based Healthy Lifestyle Score (mHLS).
#’   }
#’
#’   \item{train_medians}{
#’   Data frame of median values derived from the training dataset and used
#’   for missing-value imputation.
#’   }
#’
#’   \item{mHLS_tertile_cutoffs}{
#’   Named numeric vector containing the 33.3rd and 66.7th percentile
#’   cut-points of the overall mHLS in the training dataset.
#’   }
#’
#’   \item{domain_tertile_cuts}{
#’   Named list of numeric vectors containing the tertile cut-points for the
#’   raw domain-specific scores (Smoking, BMI, Alcohol, and Diet) derived
#’   from the training dataset.
#’   }
#’
#’   \item{CS_final_coefs}{
#’   Multinomial logistic regression coefficients for predicting the
#’   probability of being a current smoker.
#’   }
#’
#’   \item{FS_final_coefs}{
#’   Multinomial logistic regression coefficients for predicting the
#’   probability of being a former smoker.
#’   }
#’
#’   \item{NS_final_coefs}{
#’   Multinomial logistic regression coefficients for predicting the
#’   probability of being a never smoker.
#’   }
#’ }
NULL
