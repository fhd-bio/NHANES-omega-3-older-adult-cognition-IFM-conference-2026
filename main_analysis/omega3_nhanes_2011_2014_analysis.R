#!/usr/bin/env Rscript

# NHANES 2011-2014 omega-3 and cognition analysis
#
# Constructs all analytic variables from the raw G and H cycle XPT files.
# Primary definitions:
#   1. Omega-3 supplement use: dedicated past-30-day supplement interview.
#   2. Cognition: complete four-component battery.
#   3. Dietary EPA+DHA: reliable day-one recall, analyzed with WTDRD1.
#
# Sensitivity analyses include an alternative recall-day supplement definition,
# exclude self-reported stroke, add lifestyle covariates, and mutually adjust
# dietary intake for supplement use and total energy.
#
# Run from the project directory:
#   Rscript output/omega3_nhanes_2011_2014_analysis.R upload r_output

required_packages <- c("haven", "survey", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install these R packages before running: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(haven)
  library(survey)
  library(dplyr)
})

options(survey.lonely.psu = "adjust")

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "upload"
output_dir <- if (length(args) >= 2) args[[2]] else "r_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cycles <- c("G", "H")
continuous_covariates <- c("RIDAGEYR", "INDFMPIR", "BMXBMI", "phq9_total")
categorical_covariates <- c("RIAGENDR", "RIDRETH1", "DMDEDUC2", "DIQ010")
cognitive_components <- c("cerad_total", "CFDCSR", "CFDAST", "CFDDS")

omega_pattern <- paste(
  "fish oil", "salmon oil", "marine oil", "cod liver", "krill",
  "flax", "flaxseed", "linseed", "algae", "algal", "vegetarian dha",
  "vegan omega", "lovaza", "vascepa", "epanova", "omega ?-?3",
  sep = "|"
)

resolve_xpt <- function(stem, cycle) {
  pattern <- paste0(
    "^", stem, "_", cycle,
    "( ?\\([0-9]+\\))?\\.xpt$"
  )
  candidates <- list.files(
    input_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE
  )
  if (length(candidates) == 0) {
    stop("Missing input file for ", stem, "_", cycle, call. = FALSE)
  }
  if (length(candidates) > 1) {
    hashes <- unname(tools::md5sum(candidates))
    if (length(unique(hashes)) == 1) {
      message(
        "Multiple copies found for ", stem, "_", cycle,
        "; using ", basename(candidates[[1]])
      )
    } else {
      stop(
        "Multiple non-identical candidates found for ", stem, "_", cycle,
        ": ", paste(basename(candidates), collapse = ", "),
        call. = FALSE
      )
    }
  }
  candidates[[1]]
}

read_cycle_file <- function(stem, cycle) {
  haven::read_xpt(resolve_xpt(stem, cycle))
}

replace_special <- function(x, codes = c(7, 9)) {
  x[x %in% codes] <- NA
  x
}

omega_ids <- function(data) {
  names <- ifelse(is.na(data$DSDSUPP), "", as.character(data$DSDSUPP))
  unique(data$SEQN[grepl(omega_pattern, names, ignore.case = TRUE)])
}

classify_formulation <- function(name) {
  if (is.na(name)) return(NA_character_)
  rules <- list(
    prescription = "lovaza|vascepa|epanova|omega ?-?3 ?-?acid ?ethyl",
    krill = "krill",
    algal = "algae|algal|vegetarian dha|vegan omega",
    cod_liver = "cod liver",
    plant = "flax|flaxseed|linseed",
    fish_oil = "fish oil|salmon oil|marine oil|omega ?-?3"
  )
  for (label in names(rules)) {
    if (grepl(rules[[label]], name, ignore.case = TRUE)) return(label)
  }
  NA_character_
}

collapse_formulation <- function(user, labels) {
  if (user == 0) return("Non-user")
  parsed <- if (is.na(labels)) character(0) else strsplit(labels, ",")[[1]]
  if (
    length(parsed) > 0 &&
    all(parsed %in% c("fish_oil", "cod_liver"))
  ) return("Fish oil/cod liver oil")
  if (identical(parsed, "krill")) return("Krill")
  if (identical(parsed, "plant")) return("Plant-based")
  if (length(parsed) > 0) return("Other/mixed")
  "Unclassified"
}

build_cycle <- function(cycle) {
  demo <- read_cycle_file("DEMO", cycle)
  cfq <- read_cycle_file("CFQ", cycle)
  diet <- read_cycle_file("DR1TOT", cycle)
  bmx <- read_cycle_file("BMX", cycle)
  diq <- read_cycle_file("DIQ", cycle)
  dpq <- read_cycle_file("DPQ", cycle)
  mcq <- read_cycle_file("MCQ", cycle)
  smq <- read_cycle_file("SMQ", cycle)
  alq <- read_cycle_file("ALQ", cycle)
  paq <- read_cycle_file("PAQ", cycle)
  ds1 <- read_cycle_file("DS1IDS", cycle)
  ds2 <- read_cycle_file("DS2IDS", cycle)
  ds30 <- read_cycle_file("DSQIDS", cycle)
  ds30_total <- read_cycle_file("DSQTOT", cycle)

  recall_products <- bind_rows(
    select(ds1, SEQN, DSDSUPP),
    select(ds2, SEQN, DSDSUPP)
  )
  ids_24h <- omega_ids(recall_products)
  ids_30d <- omega_ids(ds30)

  omega_products <- ds30 %>%
    filter(SEQN %in% ids_30d) %>%
    transmute(
      SEQN,
      DSDSUPP = as.character(DSDSUPP),
      formulation = vapply(DSDSUPP, classify_formulation, character(1))
    )

  formulation_by_person <- omega_products %>%
    filter(!is.na(formulation)) %>%
    group_by(SEQN) %>%
    summarise(
      omega_formulations_30d = paste(sort(unique(formulation)), collapse = ","),
      .groups = "drop"
    )

  cognition <- cfq %>%
    select(SEQN, CFDCST1, CFDCST2, CFDCST3, CFDCSR, CFDAST, CFDDS) %>%
    mutate(
      cerad_total = ifelse(
        complete.cases(CFDCST1, CFDCST2, CFDCST3),
        CFDCST1 + CFDCST2 + CFDCST3,
        NA_real_
      ),
      cerad_total_partial = rowSums(
        cbind(CFDCST1, CFDCST2, CFDCST3), na.rm = TRUE
      ),
      cerad_partial_n = rowSums(
        !is.na(cbind(CFDCST1, CFDCST2, CFDCST3))
      ),
      cerad_total_partial = ifelse(
        cerad_partial_n >= 1, cerad_total_partial, NA_real_
      )
    )

  accepted_components <- c(
    "cerad_total_partial", "CFDCSR", "CFDAST", "CFDDS"
  )
  for (variable in accepted_components) {
    z_name <- paste0(variable, "_z_accepted")
    cognition[[z_name]] <- as.numeric(scale(cognition[[variable]]))
  }
  cognition$global_cognition_accepted <- rowMeans(
    cognition[paste0(accepted_components, "_z_accepted")],
    na.rm = TRUE
  )
  cognition$global_cognition_accepted[
    rowSums(!is.na(cognition[paste0(accepted_components, "_z_accepted")])) == 0
  ] <- NA_real_

  dpq_items <- sprintf("DPQ0%d0", 1:9)
  depression <- dpq %>% select(SEQN, all_of(dpq_items))
  answered <- as.data.frame(depression[dpq_items])
  answered[] <- lapply(answered, replace_special)
  answered_count <- rowSums(!is.na(answered))
  depression$phq9_total <- ifelse(
    answered_count >= 7,
    rowSums(answered, na.rm = TRUE) * 9 / answered_count,
    NA_real_
  )

  smoking <- smq %>%
    transmute(
      SEQN,
      smoking_status = case_when(
        SMQ020 == 2 ~ "never",
        SMQ020 == 1 & SMQ040 == 3 ~ "former",
        SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "current",
        TRUE ~ NA_character_
      )
    )

  alcohol <- alq %>%
    transmute(
      SEQN,
      current_alcohol = case_when(
        ALQ101 == 2 | ALQ120Q == 0 ~ 0,
        dplyr::between(ALQ120Q, 1, 366) ~ 1,
        TRUE ~ NA_real_
      )
    )

  activity <- paq %>%
    transmute(
      SEQN,
      recreational_activity = case_when(
        PAQ650 == 1 | PAQ665 == 1 ~ 1,
        PAQ650 == 2 & PAQ665 == 2 ~ 0,
        TRUE ~ NA_real_
      )
    )

  data <- demo %>%
    select(
      SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, DMDEDUC2, INDFMPIR,
      SDMVPSU, SDMVSTRA, WTMEC2YR
    ) %>%
    inner_join(cognition, by = "SEQN") %>%
    left_join(
      diet %>% select(
        SEQN, WTDRD1, DR1DRSTZ, DR1TKCAL, DR1TP205, DR1TP226
      ),
      by = "SEQN"
    ) %>%
    left_join(bmx %>% select(SEQN, BMXBMI), by = "SEQN") %>%
    left_join(diq %>% select(SEQN, DIQ010), by = "SEQN") %>%
    left_join(depression %>% select(SEQN, phq9_total), by = "SEQN") %>%
    left_join(mcq %>% select(SEQN, MCQ160F), by = "SEQN") %>%
    left_join(smoking, by = "SEQN") %>%
    left_join(alcohol, by = "SEQN") %>%
    left_join(activity, by = "SEQN") %>%
    left_join(ds30_total %>% select(SEQN, DSD010), by = "SEQN") %>%
    left_join(formulation_by_person, by = "SEQN") %>%
    filter(RIDAGEYR >= 60) %>%
    mutate(
      cycle = cycle,
      omega_user_24h = as.integer(SEQN %in% ids_24h),
      omega_user_30d = as.integer(SEQN %in% ids_30d),
      dietary_epa_dha_g = DR1TP205 + DR1TP226,
      log_dietary_epa_dha = log1p(dietary_epa_dha_g),
      DMDEDUC2 = replace_special(DMDEDUC2),
      DIQ010 = replace_special(DIQ010),
      MCQ160F = replace_special(MCQ160F),
      WTMEC4YR = WTMEC2YR / 2,
      WTDRD1_4YR = WTDRD1 / 2
    )

  list(data = data, products = omega_products)
}

message("Reading and combining NHANES cycles...")
built <- lapply(cycles, build_cycle)
data <- bind_rows(lapply(built, `[[`, "data"))
products <- bind_rows(lapply(built, `[[`, "products"))

# Standardize the corrected components across the combined 2011-2014 sample.
complete_battery <- complete.cases(data[cognitive_components])
for (variable in cognitive_components) {
  z_name <- paste0(variable, "_z")
  reference <- data[[variable]][complete_battery]
  data[[z_name]] <- (data[[variable]] - mean(reference)) / sd(reference)
}
data$global_cognition_complete <- rowMeans(
  data[paste0(cognitive_components, "_z")],
  na.rm = FALSE
)

# Factor coding is set once so model matrices are stable across analyses.
data <- data %>%
  mutate(
    some_college_or_higher = as.integer(DMDEDUC2 %in% c(4, 5)),
    current_smoker = as.integer(smoking_status == "current"),
    RIAGENDR = factor(RIAGENDR),
    RIDRETH1 = factor(RIDRETH1),
    DMDEDUC2 = factor(DMDEDUC2),
    DIQ010 = factor(DIQ010),
    smoking_status = factor(
      smoking_status, levels = c("never", "former", "current")
    ),
    omega_group_30d = factor(
      omega_user_30d, levels = c(0, 1), labels = c("Non-user", "User")
    ),
    formulation_group = factor(
      mapply(
        collapse_formulation, omega_user_30d, omega_formulations_30d,
        USE.NAMES = FALSE
      ),
      levels = c(
        "Fish oil/cod liver oil", "Krill", "Plant-based", "Other/mixed",
        "Non-user", "Unclassified"
      )
    )
  )

complete_rows <- function(frame, variables) {
  frame[complete.cases(frame[variables]), , drop = FALSE]
}

common_variables <- c(continuous_covariates, categorical_covariates)

accepted <- complete_rows(
  data,
  c(
    "global_cognition_accepted", "omega_user_24h", "WTMEC4YR",
    common_variables
  )
)

supplement_primary <- complete_rows(
  data[complete_battery, ],
  c(
    "global_cognition_complete", "omega_user_30d", "WTMEC4YR",
    common_variables
  )
)

dietary_primary <- complete_rows(
  data[
    complete_battery & !is.na(data$DR1DRSTZ) & data$DR1DRSTZ == 1,
  ],
  c(
    "global_cognition_complete", "log_dietary_epa_dha",
    "dietary_epa_dha_g", "DR1TKCAL", "WTDRD1_4YR",
    "omega_user_30d", common_variables
  )
)

make_design <- function(frame, weight_variable) {
  svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = as.formula(paste0("~", weight_variable)),
    nest = TRUE,
    data = frame
  )
}

accepted_design <- make_design(accepted, "WTMEC4YR")
supplement_design <- make_design(supplement_primary, "WTMEC4YR")
dietary_design <- make_design(dietary_primary, "WTDRD1_4YR")

adjustment_terms <- paste(
  c(continuous_covariates, categorical_covariates), collapse = " + "
)

model_accepted <- svyglm(
  as.formula(
    paste(
      "global_cognition_accepted ~ omega_user_24h +", adjustment_terms
    )
  ),
  design = accepted_design
)

model_supplement <- svyglm(
  as.formula(
    paste(
      "global_cognition_complete ~ omega_user_30d +", adjustment_terms
    )
  ),
  design = supplement_design
)

# Exploratory comparison among omega-3 users only. Rare single formulations
# are pooled as "Other/mixed" to avoid attempting to estimate singleton groups.
formulation_users <- supplement_primary %>%
  filter(omega_user_30d == 1) %>%
  mutate(formulation_group = droplevels(formulation_group))
if (any(formulation_users$formulation_group == "Unclassified")) {
  stop("At least one omega-3 user could not be classified", call. = FALSE)
}
formulation_design <- make_design(formulation_users, "WTMEC4YR")
model_formulation <- svyglm(
  as.formula(
    paste(
      "global_cognition_complete ~ formulation_group +", adjustment_terms
    )
  ),
  design = formulation_design
)

model_dietary <- svyglm(
  as.formula(
    paste(
      "global_cognition_complete ~ log_dietary_epa_dha +",
      adjustment_terms
    )
  ),
  design = dietary_design
)

model_joint_energy <- svyglm(
  as.formula(
    paste(
      "global_cognition_complete ~ log_dietary_epa_dha +",
      "omega_user_30d + DR1TKCAL +", adjustment_terms
    )
  ),
  design = dietary_design
)

stroke_free <- supplement_primary[
  is.na(supplement_primary$MCQ160F) | supplement_primary$MCQ160F != 1,
]
stroke_design <- make_design(stroke_free, "WTMEC4YR")
model_stroke_free <- svyglm(
  as.formula(
    paste(
      "global_cognition_complete ~ omega_user_30d +", adjustment_terms
    )
  ),
  design = stroke_design
)

lifestyle_variables <- c(
  "smoking_status", "current_alcohol", "recreational_activity"
)
lifestyle <- complete_rows(
  supplement_primary, lifestyle_variables
)
lifestyle_design <- make_design(lifestyle, "WTMEC4YR")
model_lifestyle <- svyglm(
  as.formula(
    paste(
      "global_cognition_complete ~ omega_user_30d +",
      "current_alcohol + recreational_activity + smoking_status +",
      adjustment_terms
    )
  ),
  design = lifestyle_design
)

tidy_model <- function(fitted_model, model_name, focus_terms = NULL) {
  coefficient_table <- as.data.frame(summary(fitted_model)$coefficients)
  coefficient_table$term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL
  names(coefficient_table)[1:4] <- c(
    "estimate", "standard_error", "statistic", "p_value"
  )
  # Calculate model-level values before mutate(). If the output column is also
  # named "model", dplyr otherwise masks the fitted model object with that
  # newly created character column.
  model_n <- length(residuals(fitted_model))
  model_residual_df <- fitted_model$df.residual
  critical <- qt(0.975, df = model_residual_df)
  coefficient_table <- coefficient_table %>%
    mutate(
      model = .env$model_name,
      n = .env$model_n,
      residual_df = .env$model_residual_df,
      ci_low = estimate - critical * standard_error,
      ci_high = estimate + critical * standard_error
    ) %>%
    select(
      model, term, n, estimate, standard_error, ci_low, ci_high,
      statistic, p_value, residual_df
    )
  if (!is.null(focus_terms)) {
    coefficient_table <- coefficient_table %>%
      filter(term %in% focus_terms)
  }
  coefficient_table
}

key_models <- bind_rows(
  tidy_model(model_accepted, "Original 24-hour definition", "omega_user_24h"),
  tidy_model(
    model_supplement, "Primary 30-day supplement model", "omega_user_30d"
  ),
  tidy_model(
    model_dietary, "Primary dietary model", "log_dietary_epa_dha"
  ),
  tidy_model(
    model_joint_energy, "Dietary + energy + supplement",
    c("log_dietary_epa_dha", "omega_user_30d")
  ),
  tidy_model(
    model_stroke_free, "Stroke-excluded sensitivity", "omega_user_30d"
  ),
  tidy_model(
    model_lifestyle, "Lifestyle-adjusted sensitivity", "omega_user_30d"
  )
)

formulation_terms <- grep(
  "^formulation_group", names(coef(model_formulation)), value = TRUE
)
formulation_results <- tidy_model(
  model_formulation, "Exploratory formulation model", formulation_terms
)
formulation_beta <- coef(model_formulation)[formulation_terms]
formulation_covariance <- vcov(model_formulation)[
  formulation_terms, formulation_terms, drop = FALSE
]
formulation_df1 <- length(formulation_terms)
formulation_f <- as.numeric(
  t(formulation_beta) %*% solve(formulation_covariance, formulation_beta) /
    formulation_df1
)
formulation_test_output <- data.frame(
  test = "Overall formulation difference among omega-3 users",
  reference = "Fish oil/cod liver oil",
  f_statistic = formulation_f,
  numerator_df = formulation_df1,
  denominator_df = model_formulation$df.residual,
  p_value = pf(
    formulation_f, formulation_df1, model_formulation$df.residual,
    lower.tail = FALSE
  )
)
formulation_counts <- formulation_users %>%
  count(formulation_group, name = "unweighted_n") %>%
  arrange(desc(unweighted_n))

# Domain-specific exploratory models, with Benjamini-Hochberg FDR correction.
domain_labels <- c(
  cerad_total_z = "CERAD immediate learning",
  CFDCSR_z = "CERAD delayed recall",
  CFDAST_z = "Animal Fluency",
  CFDDS_z = "DSST"
)

domain_results <- list()
for (outcome in names(domain_labels)) {
  supplement_domain <- svyglm(
    as.formula(
      paste(outcome, "~ omega_user_30d +", adjustment_terms)
    ),
    design = supplement_design
  )
  dietary_domain <- svyglm(
    as.formula(
      paste(outcome, "~ log_dietary_epa_dha +", adjustment_terms)
    ),
    design = dietary_design
  )
  domain_results[[paste0(outcome, "_supp")]] <- tidy_model(
    supplement_domain, "Supplement use", "omega_user_30d"
  ) %>% mutate(domain = domain_labels[[outcome]])
  domain_results[[paste0(outcome, "_diet")]] <- tidy_model(
    dietary_domain, "Dietary EPA+DHA", "log_dietary_epa_dha"
  ) %>% mutate(domain = domain_labels[[outcome]])
}
domain_results <- bind_rows(domain_results) %>%
  group_by(model) %>%
  mutate(fdr_q_value = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  select(
    model, domain, term, n, estimate, standard_error, ci_low, ci_high,
    p_value, fdr_q_value, residual_df
  )

# Restricted cubic spline with four weighted knots, matching the Python audit.
weighted_quantile <- function(x, probabilities, weights) {
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  x <- x[keep]
  weights <- weights[keep]
  order_index <- order(x)
  x <- x[order_index]
  weights <- weights[order_index]
  cumulative <- (cumsum(weights) - 0.5 * weights) / sum(weights)
  approx(
    x = cumulative, y = x, xout = probabilities,
    method = "linear", rule = 2
  )$y
}

rcs_basis <- function(x, knots) {
  if (length(knots) < 3 || any(diff(knots) <= 0)) {
    stop("Restricted cubic spline knots must be unique.", call. = FALSE)
  }
  last <- knots[length(knots)]
  penultimate <- knots[length(knots) - 1]
  scale_value <- (last - knots[1])^2
  positive_cube <- function(value) pmax(value, 0)^3
  columns <- list(x)
  for (knot in knots[seq_len(length(knots) - 2)]) {
    term <- positive_cube(x - knot)
    term <- term -
      ((last - knot) / (last - penultimate)) *
      positive_cube(x - penultimate)
    term <- term +
      ((penultimate - knot) / (last - penultimate)) *
      positive_cube(x - last)
    columns[[length(columns) + 1]] <- term / scale_value
  }
  do.call(cbind, columns)
}

spline_knots <- weighted_quantile(
  dietary_primary$log_dietary_epa_dha,
  c(0.05, 0.35, 0.65, 0.95),
  dietary_primary$WTDRD1_4YR
)
spline_matrix <- rcs_basis(
  dietary_primary$log_dietary_epa_dha, spline_knots
)
colnames(spline_matrix) <- c("diet_rcs_linear", "diet_rcs_nl1", "diet_rcs_nl2")
dietary_spline_data <- bind_cols(
  dietary_primary, as.data.frame(spline_matrix)
)
dietary_spline_design <- make_design(dietary_spline_data, "WTDRD1_4YR")
model_spline <- svyglm(
  as.formula(
    paste(
      "global_cognition_complete ~ diet_rcs_linear + diet_rcs_nl1 +",
      "diet_rcs_nl2 + omega_user_30d + DR1TKCAL +", adjustment_terms
    )
  ),
  design = dietary_spline_design
)
spline_nonlinear_terms <- c("diet_rcs_nl1", "diet_rcs_nl2")
spline_beta <- coef(model_spline)[spline_nonlinear_terms]
spline_covariance <- vcov(model_spline)[
  spline_nonlinear_terms, spline_nonlinear_terms, drop = FALSE
]
spline_f <- as.numeric(
  t(spline_beta) %*% solve(spline_covariance) %*% spline_beta /
    length(spline_nonlinear_terms)
)
spline_p <- pf(
  spline_f,
  df1 = length(spline_nonlinear_terms),
  df2 = model_spline$df.residual,
  lower.tail = FALSE
)
spline_test_output <- data.frame(
  test = "Restricted cubic spline nonlinearity",
  statistic = spline_f,
  numerator_df = length(spline_nonlinear_terms),
  denominator_df = model_spline$df.residual,
  p_value = spline_p,
  knot_5 = spline_knots[1],
  knot_35 = spline_knots[2],
  knot_65 = spline_knots[3],
  knot_95 = spline_knots[4]
)

# Survey-weighted two-group baseline table with standardized mean differences.
weighted_continuous_row <- function(variable, label) {
  group_designs <- list(
    nonuser = subset(supplement_design, omega_user_30d == 0),
    user = subset(supplement_design, omega_user_30d == 1)
  )
  means <- vapply(
    group_designs,
    function(design) as.numeric(coef(svymean(
      as.formula(paste0("~", variable)), design, na.rm = TRUE
    ))),
    numeric(1)
  )
  variances <- vapply(
    group_designs,
    function(design) as.numeric(svyvar(
      as.formula(paste0("~", variable)), design, na.rm = TRUE
    )),
    numeric(1)
  )
  smd <- (means[["user"]] - means[["nonuser"]]) /
    sqrt((variances[["user"]] + variances[["nonuser"]]) / 2)
  data.frame(
    characteristic = label,
    level = "",
    nonuser = sprintf("%.1f (%.1f)", means[["nonuser"]], sqrt(variances[["nonuser"]])),
    user = sprintf("%.1f (%.1f)", means[["user"]], sqrt(variances[["user"]])),
    standardized_difference = abs(unname(smd))
  )
}

weighted_binary_row <- function(variable, level_value, label) {
  indicator_name <- paste0("baseline_indicator_", variable, "_", level_value)
  supplement_design$variables[[indicator_name]] <-
    as.numeric(
      as.character(supplement_design$variables[[variable]]) ==
        as.character(level_value)
    )
  group_designs <- list(
    nonuser = subset(supplement_design, omega_user_30d == 0),
    user = subset(supplement_design, omega_user_30d == 1)
  )
  proportions <- vapply(
    group_designs,
    function(design) as.numeric(coef(svymean(
      as.formula(paste0("~", indicator_name)), design, na.rm = TRUE
    ))),
    numeric(1)
  )
  pooled <- (proportions[["user"]] + proportions[["nonuser"]]) / 2
  smd <- if (pooled > 0 && pooled < 1) {
    (proportions[["user"]] - proportions[["nonuser"]]) /
      sqrt(pooled * (1 - pooled))
  } else {
    NA_real_
  }
  data.frame(
    characteristic = label,
    level = as.character(level_value),
    nonuser = sprintf("%.1f%%", 100 * proportions[["nonuser"]]),
    user = sprintf("%.1f%%", 100 * proportions[["user"]]),
    standardized_difference = abs(unname(smd))
  )
}

baseline_table <- bind_rows(
  weighted_continuous_row("RIDAGEYR", "Age, years"),
  weighted_binary_row("RIAGENDR", 2, "Female sex"),
  weighted_binary_row(
    "some_college_or_higher", 1,
    "Some college/associate degree or higher"
  ),
  weighted_continuous_row("INDFMPIR", "Family income-to-poverty ratio"),
  weighted_continuous_row("BMXBMI", "BMI, kg/m²"),
  weighted_binary_row("DIQ010", 1, "Diagnosed diabetes"),
  weighted_continuous_row("phq9_total", "PHQ-9 score"),
  weighted_binary_row("current_smoker", 1, "Current smoker"),
  weighted_binary_row(
    "recreational_activity", 1,
    "Any moderate/vigorous recreational activity"
  )
)

unweighted_counts <- supplement_primary %>%
  count(omega_group_30d, name = "unweighted_n")

sample_flow <- data.frame(
  stage = c(
    "Adults aged 60 years or older with cognition file linkage",
    "Original complete-case analysis",
    "Complete four-component battery and primary covariates",
    "Reliable day-one dietary recall and dietary model covariates",
    "Complete lifestyle sensitivity-analysis covariates"
  ),
  n = c(
    nrow(data),
    nrow(accepted),
    nrow(supplement_primary),
    nrow(dietary_primary),
    nrow(lifestyle)
  )
)

# Adjusted spline contrast relative to 100 mg/day.
grid_mg <- seq(0, 1000, by = 5)
grid_log_g <- log1p(grid_mg / 1000)
grid_basis <- rcs_basis(grid_log_g, spline_knots)
reference_basis <- as.numeric(rcs_basis(log1p(0.1), spline_knots))
spline_terms <- c("diet_rcs_linear", "diet_rcs_nl1", "diet_rcs_nl2")
spline_contrast <- sweep(grid_basis, 2, reference_basis, FUN = "-")
spline_coefficients <- coef(model_spline)[spline_terms]
spline_term_covariance <- vcov(model_spline)[
  spline_terms, spline_terms, drop = FALSE
]
spline_estimate <- as.numeric(spline_contrast %*% spline_coefficients)
spline_variance <- rowSums(
  (spline_contrast %*% spline_term_covariance) * spline_contrast
)
spline_standard_error <- sqrt(pmax(spline_variance, 0))
spline_critical <- qt(0.975, df = model_spline$df.residual)
spline_lower <- spline_estimate - spline_critical * spline_standard_error
spline_upper <- spline_estimate + spline_critical * spline_standard_error

png(
  file.path(output_dir, "r_dietary_spline.png"),
  width = 1800, height = 1150, res = 220
)
plot(
  grid_mg, spline_estimate, type = "n",
  ylim = range(c(spline_lower, spline_upper), finite = TRUE),
  xlab = "Dietary EPA + DHA (mg/day)",
  ylab = "Adjusted difference in cognitive z-score\n(reference: 100 mg/day)",
  main = sprintf(
    "Exploratory dose-response (P for nonlinearity = %.3f)", spline_p
  ),
  bty = "l"
)
polygon(
  c(grid_mg, rev(grid_mg)),
  c(spline_lower, rev(spline_upper)),
  col = "#BFD7EACC", border = NA
)
lines(grid_mg, spline_estimate, col = "#1261A0", lwd = 3)
abline(h = 0, col = "#444444", lty = 2)
abline(v = 100, col = "#888888", lty = 3)
dev.off()

write.csv(
  key_models, file.path(output_dir, "r_key_models.csv"), row.names = FALSE
)
write.csv(
  formulation_results,
  file.path(output_dir, "r_formulation_exploratory.csv"), row.names = FALSE
)
write.csv(
  formulation_test_output,
  file.path(output_dir, "r_formulation_overall_test.csv"), row.names = FALSE
)
write.csv(
  formulation_counts,
  file.path(output_dir, "r_formulation_counts.csv"), row.names = FALSE
)
write.csv(
  domain_results, file.path(output_dir, "r_domain_models.csv"),
  row.names = FALSE
)
write.csv(
  spline_test_output, file.path(output_dir, "r_spline_test.csv"),
  row.names = FALSE
)
write.csv(
  baseline_table, file.path(output_dir, "r_baseline_weighted.csv"),
  row.names = FALSE
)
write.csv(
  unweighted_counts, file.path(output_dir, "r_group_counts.csv"),
  row.names = FALSE
)
write.csv(
  sample_flow, file.path(output_dir, "r_sample_flow.csv"), row.names = FALSE
)
write.csv(
  products, file.path(output_dir, "r_omega3_products_30d.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    accepted = model_accepted,
    supplement_primary = model_supplement,
    formulation_exploratory = model_formulation,
    dietary_primary = model_dietary,
    joint_energy = model_joint_energy,
    stroke_free = model_stroke_free,
    lifestyle = model_lifestyle,
    spline = model_spline
  ),
  file.path(output_dir, "r_fitted_models.rds")
)

capture.output(
  {
    cat("CORRECTED NHANES 2011-2014 OMEGA-3 ANALYSIS\n\n")
    cat("Sample flow:\n")
    print(sample_flow, row.names = FALSE)
    cat("\nGroup counts:\n")
    print(unweighted_counts, row.names = FALSE)
    cat("\nKey models:\n")
    print(key_models, row.names = FALSE)
    cat("\nExploratory formulation subgroup counts:\n")
    print(formulation_counts, row.names = FALSE)
    cat("\nExploratory formulation coefficients:\n")
    print(formulation_results, row.names = FALSE)
    cat("\nOverall exploratory formulation test:\n")
    print(formulation_test_output, row.names = FALSE)
    cat("\nSpline nonlinearity test:\n")
    print(spline_test_output, row.names = FALSE)
    cat("\nR session:\n")
    print(sessionInfo())
  },
  file = file.path(output_dir, "r_analysis_report.txt")
)

message("Analysis complete. Results written to: ", normalizePath(output_dir))
