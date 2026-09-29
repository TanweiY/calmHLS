# calmHLS

**Calculate DNA Methylation‑Based Healthy Lifestyle Scores (mHLS)**

`calmHLS` provides a simple R interface to compute individual DNA methylation (DNAm)‑based scores for smoking, alcohol consumption, and BMI, extract diet‑related CpGs, and combine them into an overall **Healthy Lifestyle Score (mHLS)** along with domain‑specific sub‑scores. The package uses pre‑trained weights and coefficients

---

## Installation

The package can be installed directly from GitHub.  

```r

# Install calmHLS 
remotes::install_github("TanweiY/calmHLS")
```
---

## Required input data

- **Methylation data (`beta`)**: a `data.frame` with a column **`ID`** (unique sample identifier) and columns for CpG sites (beta values between 0 and 1). Samples in rows, CpGs in columns.
- **Sample sheet (`samplesheet`)**: a `data.frame` with at least the columns **`ID`**, **`Age`** (numeric), **`Gender`** (`"Male"` or `"Female"`), and **`sex`** (1 = male, 2 = female). 

The package includes internal tables of the 1231 CpGs required for the scores; you will see a message indicating how many of them are present in your data.

---

## Example Usage with internal example datasets

The typical workflow consists of two steps:

1. **Calculate individual behaviour scores and extract diet CpGs**
2. **Compute the overall mHLS and domain sub‑scores** – these are automatically categorised into tertile‑based **Low / Medium / High** groups using cut‑points derived from the original training set.

A convenience wrapper combines both steps into one call.

```r
library(calmHLS)

# ---- Load example data (shipped with the package) ----
data("example_beta")          # Methylation matrix
data("example_samplesheet")   # Sample sheet

# ---- Option A: Step by step ----
behavior_scores <- calc_behavior_scores(example_beta, example_samplesheet)
final_scores    <- calc_mHLS(behavior_scores)
head(final_scores)

# ---- Option B: All-in-one wrapper ----
final_scores <- calmHLS(example_beta, example_samplesheet)
head(final_scores)
```

---

## Output

## Output

The final `data.frame` contains the following columns:

| Column | Description |
|--------|-------------|
| `ID` | Sample identifier |
| `mHLS` | Overall DNAm-based Healthy Lifestyle Score (linear predictor). Higher values indicate a healthier lifestyle profile. |
| `RawScore_noSmoking` | Raw domain sub-score for smoking behavior; higher values indicate healthier (less smoking exposure). |
| `RawScore_lowBMI` | Raw domain sub-score for BMI; higher values indicate a healthier BMI profile. |
| `RawScore_lowAlcohol` | Raw domain sub-score for alcohol consumption; higher values indicate lower alcohol intake (healthier profile). |
| `RawScore_goodDiet` | Raw domain sub-score for diet quality; higher values indicate a healthier dietary pattern. |
| `mHLS_category` | Training-set tertile-based categories of overall mHLS: **Low**, **Medium**, **High** (higher is healthier). |
| `noSmoking_category` | Training-set tertile-based category for the smoking domain (higher indicates healthier behavior). |
| `lowBMI_category` | Training-set tertile-based category for the BMI domain (higher indicates healthier BMI). |
| `lowAlcohol_category` | Training-set tertile-based category for the alcohol domain (higher indicates lower alcohol intake). |
| `goodDiet_category` | Training-set tertile-based category for the diet domain (higher indicates healthier diet). |

### Category definition

The category thresholds are defined using the 33.3% and 66.7% percentiles of each score in the original training dataset. These cutoffs are stored internally and consistently applied to all new datasets.

### Additional output

The intermediate object returned by `calc_behavior_scores()` also includes additional components, such as:
- Individual model components (e.g., `mSmk_233`, `probs_fsmk`, `mAlc_450`, `mBMI_397`)
- CpG-level variables used in the diet-related score calculation

---

## Package features

- **Automatic imputation**: missing predictors (e.g., diet CpGs absent from your methylation data) are filled with training‑set medians when computing the mHLS.
- **Tertile‑based categorisation**: both the overall mHLS and each domain sub‑score are automatically classified into **Low / Medium / High** groups using fixed training‑set cut‑points.
- **Transparent messaging**: the number of available CpGs for each score and the number of imputed features are reported.
---

## Dependencies

- **R ≥ 3.5**
- **stats** (base R)

---

## License
This software is proprietary and **not open source**. All rights reserved.

The `calmHLS` package is made available for **non-commercial academic and research use only**. Commercial use, modification, redistribution, and derivative works are not permitted without prior written permission from the copyright holder.

No patent license, express or implied, is granted. The copyright holder reserves all patent rights.

See the [`LICENSE`](LICENSE) file for the full terms.

---

## Citation

If you use `calmHLS` in your work, please cite the relevant DNAm score publications (to be added).

## R shiny app
You can also calculate the scores on the website without R: https://tanweiyuan.shinyapps.io/mHLS_calculator/
