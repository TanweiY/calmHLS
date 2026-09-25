#' Calculate individual DNAm-based scores and extract diet CpGs
#'
#' Computes mSmk-233, mAlc-450, mBMI-397, optionally smoking probability
#' (probs_fsmk) via EpiSmokEr, and extracts the diet-related CpG columns.
#'
#' @param beta_df A data.frame with a column \code{ID} and CpG columns (beta values).
#'   Samples in rows, CpGs in columns.
#'
#' @param samplesheet A data.frame with columns \code{ID}, \code{Age}, \code{Gender}
#'   (\code{"Male"/"Female"}), and \code{sex} (1=male, 2=female).
#'
#'
#' @return A data.frame containing:
#'   \item{ID}{Sample identifiers.}
#'   \item{mSmk_233}{DNAm smoking score (233 CpGs).}
#'   \item{probs_fsmk}{Predicted probability of being a former smoker.}
#'   \item{mAlc_450}{DNAm alcohol score (450 CpGs).}
#'   \item{mBMI_397}{DNAm BMI score (397 CpGs).}
#'   \item{...}{All available diet-related CpG columns (from \code{diet_cpgs}).}
#'
#' @export
#' @examples
#' \dontrun{
#'   data(beta_df)        # user's methylation data
#'   data(base_SSt)           # samplesheet with Age, Gender, sex
#'   scores <- calc_behavior_scores(beta_df, base_SSt)
#' }

calc_behavior_scores <- function(beta_df, samplesheet) {

   # validate input
  if (!"ID" %in% colnames(beta_df)) {
    stop("The input data.frame must contain a column named 'ID'.", call. = FALSE)
  }


  if (is.null(samplesheet)) {
    stop("A 'samplesheet' with ID, Age, Gender, and sex is required for calculating probs_fsmk")}

  ids <- beta_df$ID

  # --- determine available CpGs for each score ---
  cpg_avail <- intersect(colnames(beta_df), all_required_cpgs)
  message(length(cpg_avail), " out of ", length(all_required_cpgs),
          " required CpGs available in the beta matrix.")

  beta_sub <- beta_df[, c("ID", cpg_avail), drop = FALSE]
  # convert to matrix for calculations (without ID)
  beta_mat <- as.matrix(beta_sub[, -1, drop = FALSE])
  rownames(beta_mat) <- ids

  # --- 1. probs_fsmk ---
  samplesheet <- samplesheet[match(ids, samplesheet$ID), ]
  rownames(samplesheet)<-samplesheet$ID

  beta_t <- t(beta_mat)  # rows CpGs, columns samples
  #  expects samplesheet with ID, Age, Gender, sex
  epi_res <- SSt(dataset = beta_t,
                 samplesheet = samplesheet,
                 ref.CS = CS_final_coefs,
                 ref.FS = FS_final_coefs,
                 ref.NS = NS_final_coefs)

  # --- 2. mSmk-233 ---
  smk_cpgs <- intersect(colnames(beta_sub), smk233_coef$CpG)
  message("Calculating mSmk-233: ", length(smk_cpgs), " out of ", nrow(smk233_coef),
          " CpGs available.")
  smk_scores <- smk233_coef[match(smk_cpgs, smk233_coef$CpG), ]
  mSmk_233 <- as.vector(beta_mat[, smk_scores$CpG, drop = FALSE] %*% smk_scores$Weight)

  # --- 3. mAlc-450 ---
  alc_cpgs <- intersect(colnames(beta_sub), alcohol_coef$CpG)
  message("Calculating mAlc-450: ", length(alc_cpgs), " out of ", nrow(alcohol_coef),
          " CpGs available.")
  alc_scores <- alcohol_coef[match(alc_cpgs, alcohol_coef$CpG), ]
  mAlc_450 <- as.vector(beta_mat[, alc_scores$CpG, drop = FALSE] %*% alc_scores$Weight)

  # --- 4. mBMI-397 ---
  bmi_cpgs <- intersect(colnames(beta_sub), bmi_coef$CpG)
  message("Calculating mBMI-397: ", length(bmi_cpgs), " out of ", nrow(bmi_coef),
          " CpGs available.")
  bmi_scores <- bmi_coef[match(bmi_cpgs, bmi_coef$CpG), ]
  mBMI_397 <- as.vector(beta_mat[, bmi_scores$CpG, drop = FALSE] %*% bmi_scores$Weight) + 1.45985838513289

  # --- assemble result ---
  out <- data.frame(ID = ids,
                    mSmk_233 = mSmk_233,
                    probs_fsmk = epi_res[, "probs_FS"],
                    mAlc_450 = mAlc_450,
                    mBMI_397 = mBMI_397,
                    row.names = NULL,
                    stringsAsFactors = FALSE)

  # --- 5. Extract diet CpGs ---
  diet_avail <- intersect(colnames(beta_sub), diet_cpgs)
  message( length(diet_avail), " out of ", length(diet_cpgs),
           " diet-related CpGs available.")

  if (length(diet_avail) > 0) {
    diet_df <- beta_sub[, c("ID", diet_avail), drop = FALSE]
    # remove ID to avoid duplicate column
    out <- merge(out, diet_df, by = "ID", all.x = TRUE)
  }

  return(out)

}
