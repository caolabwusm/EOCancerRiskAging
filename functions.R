# Publication-ready helper functions for
# Generational Shifts in Early-Onset Cancer Risk

suppressPackageStartupMessages({
  library(dplyr)
  library(rlang)
  library(stringr)
  library(purrr)
  library(survival)
  library(broom)
  library(tidyr)
  library(BioAge)
})

load_core_packages <- function() {
  pkgs <- c(
    "tidyverse", "data.table", "lubridate", "survival", "broom",
    "mgcv", "BioAge", "openxlsx", "patchwork", "table1"
  )
  invisible(lapply(pkgs, require, character.only = TRUE))
}

assert_has_columns <- function(data, cols, data_name = deparse(substitute(data))) {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "%s is missing required columns: %s",
      data_name,
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

standardize_vec <- function(x) {
  as.numeric((x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE))
}

add_quantile_group <- function(x, probs = c(0, 1/3, 2/3, 1), labels = NULL) {
  if (is.null(labels)) labels <- seq_len(length(probs) - 1)
  cut(x, breaks = quantile(x, probs = probs, na.rm = TRUE), include.lowest = TRUE, labels = labels)
}

remove_sd_outliers <- function(data, cols, threshold = 5) {
  assert_has_columns(data, cols)
  keep <- rowSums(abs(scale(data[, cols])) < threshold) == length(cols)
  data[keep, , drop = FALSE]
}

compute_age_acceleration <- function(data, age_col, ba_col, prefix) {
  age_col <- enquo(age_col)
  ba_col  <- enquo(ba_col)

  model <- stats::lm(
    stats::reformulate(attr(terms(stats::as.formula(paste0("~", quo_name(age_col)))), "term.labels"),
                       response = quo_name(ba_col)),
    data = data
  )

  aa <- residuals(model)
  std <- standardize_vec(aa)

  data %>%
    mutate(
      !!paste0("AA_", prefix) := aa,
      !!paste0("Standardized_AA_", prefix) := std,
      !!paste0("Tertiles_Standardized_AA_", prefix) := add_quantile_group(aa),
      !!paste0("Quartiles_Standardized_AA_", prefix) := add_quantile_group(aa, probs = c(0, 0.25, 0.5, 0.75, 1))
    )
}

prepare_phenoage_data <- function(data, biomarker_cols, age_col = "age") {
  assert_has_columns(data, c(age_col, biomarker_cols))
  data %>%
    filter(if_all(all_of(c(age_col, biomarker_cols)), ~ !is.na(.x))) %>%
    remove_sd_outliers(cols = biomarker_cols, threshold = 5)
}

compute_phenoage_from_biomarkers <- function(data, age_col = "age", biomarker_cols) {
  assert_has_columns(data, c(age_col, biomarker_cols))
  fit <- BioAge::phenoage_nhanes(biomarker_cols)
  out <- BioAge::phenoage_calc(
    data = data,
    biomarkers = biomarker_cols,
    fit = fit$fit,
    orig = TRUE
  )$data
  out
}

compute_kdm_from_biomarkers <- function(data, sex_col = "sex", female_label = "Female", male_label = "Male", biomarker_cols) {
  assert_has_columns(data, c(sex_col, biomarker_cols, "age"))
  fit <- BioAge::kdm_nhanes(biomarker_cols)

  female_df <- data %>% filter(.data[[sex_col]] == female_label)
  male_df   <- data %>% filter(.data[[sex_col]] == male_label)

  female_out <- BioAge::kdm_calc(
    data = female_df,
    biomarkers = biomarker_cols,
    fit = fit$fit$female,
    s_ba2 = fit$fit$female$s_b2
  )$data

  male_out <- BioAge::kdm_calc(
    data = male_df,
    biomarkers = biomarker_cols,
    fit = fit$fit$male,
    s_ba2 = fit$fit$male$s_b2
  )$data

  bind_rows(female_out, male_out)
}

make_event_age <- function(data, baseline_age_col, followup_time_col, out_col) {
  data %>% mutate(!!out_col := .data[[baseline_age_col]] + .data[[followup_time_col]])
}

build_base_covariates <- function(dataset = c("ukbb", "aou")) {
  dataset <- match.arg(dataset)
  if (dataset == "ukbb") {
    c(
      "sex", "race_new", "socStatus_new", "education_new", "bmiCont",
      "smokingCat_new", "alcoPerDayCat_new", "metTotalCat",
      "healthy_diet_new", "beforeBaselineCOPD", "db_new", "beforeBaselineCVD",
      paste0("genetic_principal_components_f22009_0_", 1:10)
    )
  } else {
    c(
      "sex", "race_ethnicity", "social_deprivation_index", "education_new",
      "bmiCont", "smokingCat_new", "alcohol_status", "physical_activity_cat",
      "cancerFamilyHistory", "beforeBaselineCOPD", "db_new", "beforeBaselineCVD"
    )
  }
}

build_formula <- function(time_start, time_stop, event, exposure, covariates = NULL) {
  rhs <- c(exposure, covariates)
  stats::as.formula(
    paste0(
      "survival::Surv(", time_start, ", ", time_stop, ", ", event, ") ~ ",
      paste(rhs, collapse = " + ")
    )
  )
}

run_single_cox <- function(data, formula, ties = "efron") {
  survival::coxph(formula, data = data, ties = ties) %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE)
}

run_cox_series <- function(data, event_vars, time_stop_var, exposure, covariates, time_start_var = "ageBaseline") {
  purrr::map_dfr(event_vars, function(event_var) {
    form <- build_formula(
      time_start = time_start_var,
      time_stop  = time_stop_var,
      event      = event_var,
      exposure   = exposure,
      covariates = covariates
    )

    run_single_cox(data = data, formula = form) %>%
      mutate(outcome = event_var, exposure = exposure)
  })
}

add_bh_fdr <- function(results, p_col = "p.value", group_cols = NULL, out_col = "p_fdr") {
  if (is.null(group_cols)) {
    results[[out_col]] <- p.adjust(results[[p_col]], method = "BH")
    return(results)
  }

  results %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(!!out_col := p.adjust(.data[[p_col]], method = "BH")) %>%
    ungroup()
}

extract_exposure_row <- function(results, exposure_pattern) {
  results %>% filter(str_detect(term, fixed(exposure_pattern)))
}

make_birth_cohort <- function(x, breaks, labels) {
  cut(x, breaks = breaks, labels = labels, right = TRUE, include.lowest = TRUE)
}

fit_birth_cohort_linear <- function(data, outcome, birth_var = "birthCohort", adjust = c("ageBaseline", "ageSquare")) {
  stats::glm(
    stats::as.formula(paste(outcome, "~", paste(c(birth_var, adjust), collapse = " + "))),
    family = gaussian(),
    data = data
  )
}

fit_birth_year_gam <- function(data, outcome, birth_year = "birthYear", adjust = c("ageBaseline", "ageSquare")) {
  mgcv::gam(
    stats::as.formula(paste(outcome, "~ s(", birth_year, ") +", paste(adjust, collapse = " + "))),
    data = data
  )
}

save_xlsx_list <- function(x, file) {
  openxlsx::write.xlsx(x, file = file, overwrite = TRUE)
}

# -------------------------------
# Publication-specific helpers
# -------------------------------

get_manuscript_outcomes <- function() {
  list(
    eo_primary = c(
      "dt_AllSolidCancerEO55", "dt_CNSEO55", "dt_HeadneckEO55",
      "dt_ThyroidEO55", "dt_LungEO55", "dt_BreastEO55",
      "dt_MelanomaEO55", "dt_GastrointestinalEO55", "dt_ColorectalEO55",
      "dt_GIexcludeCRCEO55", "dt_UterusEO55", "dt_OvarianEO55", "dt_ProstateEO55"
    ),
    eo_main_figure = c(
      "dt_AllSolidCancerEO55", "dt_LungEO55", "dt_GastrointestinalEO55",
      "dt_ColorectalEO55", "dt_GIexcludeCRCEO55", "dt_UterusEO55"
    ),
    lo_selected = c(
      "dt_LungLO55", "dt_GastrointestinalLO55", "dt_ColorectalLO55",
      "dt_GIexcludeCRCLO55", "dt_UterusLO55"
    )
  )
}

get_organ_clock_vars <- function() {
  c(
    "standardized_Brain", "standardized_Pituitary", "standardized_Salivary",
    "standardized_Thyroid", "standardized_Esophagus", "standardized_Lung",
    "standardized_Heart", "standardized_Artery", "standardized_Liver",
    "standardized_Stomach", "standardized_Pancreas", "standardized_Kidney",
    "standardized_Intestine", "standardized_Adrenal", "standardized_Immune",
    "standardized_Skin", "standardized_Muscle", "standardized_Adipose",
    "standardized_Female", "standardized_Male", "standardized_Organismal",
    "standardized_Multi.organ", "standardized_Conventional"
  )
}
