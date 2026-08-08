#!/usr/bin/env Rscript

# Prespecified secondary analyses dated 2026-07-30.
#
# This script uses the analytic dataset produced by:
#   main_analysis/omega3_analysis.py
#
# Its Taylor-linearized WLS implementation is algebraically equivalent to the
# Gaussian svyglm calculation for this stratified one-stage NHANES design. It
# first verifies the implementation against the primary-model benchmark.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) {
  args[[1]]
} else if (nzchar(Sys.getenv("OMEGA3_SECONDARY_INPUT"))) {
  Sys.getenv("OMEGA3_SECONDARY_INPUT")
} else {
  "tmp/reanalysis/rebuilt_analytic_data.csv"
}
output_dir <- if (length(args) >= 2) {
  args[[2]]
} else if (nzchar(Sys.getenv("OMEGA3_SECONDARY_OUTPUT"))) {
  Sys.getenv("OMEGA3_SECONDARY_OUTPUT")
} else {
  "secondary_analysis/output"
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(input_file, check.names = FALSE)

factor_variables <- c("RIAGENDR", "RIDRETH1", "DMDEDUC2", "DIQ010")
for (variable in factor_variables) {
  data[[variable]] <- factor(data[[variable]])
}
data$cycle <- factor(
  data$cycle,
  levels = c("G", "H"),
  labels = c("2011-2012", "2013-2014")
)
data$WTMEC4YR <- data$WTMEC2YR / 2

complete_battery <- complete.cases(
  data[c("cerad_total", "CFDCSR", "CFDAST", "CFDDS")]
)
core_variables <- c(
  "global_cognition_complete", "omega_user_30d",
  "RIDAGEYR", "INDFMPIR", "BMXBMI", "phq9_total",
  factor_variables, "SDMVPSU", "SDMVSTRA", "WTMEC2YR", "cycle"
)
primary <- data[
  complete_battery & complete.cases(data[core_variables]),
  ,
  drop = FALSE
]

adjustment_terms <- paste(
  c(
    "RIDAGEYR", "INDFMPIR", "BMXBMI", "phq9_total",
    factor_variables
  ),
  collapse = " + "
)
primary_formula <- as.formula(
  paste(
    "global_cognition_complete ~ omega_user_30d +",
    adjustment_terms
  )
)

survey_lm_taylor <- function(
  formula, frame, weight, strata = "SDMVSTRA", psu = "SDMVPSU",
  stratum_prefix = NULL
) {
  model_frame <- model.frame(formula, data = frame, na.action = na.fail)
  y <- model.response(model_frame)
  x <- model.matrix(attr(model_frame, "terms"), model_frame)
  w <- frame[[weight]]

  xtwx_inverse <- solve(crossprod(x, w * x))
  beta <- as.vector(xtwx_inverse %*% crossprod(x, w * y))
  names(beta) <- colnames(x)
  residual <- as.vector(y - x %*% beta)
  score <- x * as.vector(w * residual)

  stratum_id <- as.character(frame[[strata]])
  if (!is.null(stratum_prefix)) {
    stratum_id <- interaction(
      frame[[stratum_prefix]], stratum_id, drop = TRUE, lex.order = TRUE
    )
  }
  cluster_id <- interaction(
    stratum_id, frame[[psu]], drop = TRUE, lex.order = TRUE
  )

  cluster_scores <- rowsum(score, cluster_id, reorder = TRUE)
  cluster_strata <- tapply(
    as.character(stratum_id), cluster_id, function(value) value[[1]]
  )
  cluster_strata <- cluster_strata[rownames(cluster_scores)]

  meat <- matrix(0, nrow = ncol(x), ncol = ncol(x))
  for (stratum in unique(cluster_strata)) {
    stratum_scores <- cluster_scores[cluster_strata == stratum, , drop = FALSE]
    number_psus <- nrow(stratum_scores)
    if (number_psus < 2) {
      stop("A stratum has fewer than two PSUs.")
    }
    centered <- sweep(
      stratum_scores, 2, colMeans(stratum_scores), FUN = "-"
    )
    meat <- meat +
      (number_psus / (number_psus - 1)) * crossprod(centered)
  }

  variance <- xtwx_inverse %*% meat %*% xtwx_inverse
  standard_error <- sqrt(diag(variance))
  names(standard_error) <- colnames(x)

  number_psus <- nrow(cluster_scores)
  number_strata <- length(unique(cluster_strata))
  design_df <- number_psus - number_strata
  residual_df <- design_df - ncol(x) + 1

  list(
    coefficients = beta,
    standard_errors = standard_error,
    variance = variance,
    design_df = design_df,
    residual_df = residual_df,
    n = nrow(frame),
    rank = ncol(x)
  )
}

tidy_term <- function(model, term, label, reference = "t") {
  estimate <- unname(model$coefficients[[term]])
  standard_error <- unname(model$standard_errors[[term]])
  statistic <- estimate / standard_error

  if (reference == "normal") {
    critical <- qnorm(0.975)
    p_value <- 2 * pnorm(abs(statistic), lower.tail = FALSE)
    inference_df <- Inf
  } else {
    if (model$residual_df <= 0) {
      stop("Residual survey degrees of freedom are not positive.")
    }
    critical <- qt(0.975, df = model$residual_df)
    p_value <- 2 * pt(
      abs(statistic), df = model$residual_df, lower.tail = FALSE
    )
    inference_df <- model$residual_df
  }

  data.frame(
    analysis = label,
    term = term,
    n = model$n,
    estimate = estimate,
    standard_error = standard_error,
    ci_low = estimate - critical * standard_error,
    ci_high = estimate + critical * standard_error,
    statistic = statistic,
    p_value = p_value,
    reference_distribution = reference,
    inference_df = inference_df,
    design_df = model$design_df,
    residual_design_df = model$residual_df
  )
}

# Verify the calculation against the frozen primary model.
frozen_primary <- survey_lm_taylor(
  primary_formula, primary, "WTMEC4YR", stratum_prefix = "cycle"
)
expected_beta <- 0.150966049029382
expected_se <- 0.0351884788505344
if (
  abs(frozen_primary$coefficients[["omega_user_30d"]] - expected_beta) > 1e-10 ||
  abs(frozen_primary$standard_errors[["omega_user_30d"]] - expected_se) > 1e-10
) {
  stop(
    sprintf(
      paste(
        "Primary-model verification failed; secondary analyses were not run.",
        "Observed beta %.15f and SE %.15f."
      ),
      frozen_primary$coefficients[["omega_user_30d"]],
      frozen_primary$standard_errors[["omega_user_30d"]]
    )
  )
}

# Analysis 1: active comparator.
primary$supplement_group <- ifelse(
  primary$omega_user_30d == 1,
  "Omega-3 supplement",
  ifelse(
    primary$omega_user_30d == 0 & primary$DSD010 == 1,
    "Other supplements only",
    ifelse(
      primary$omega_user_30d == 0 & primary$DSD010 == 2,
      "No supplements",
      NA_character_
    )
  )
)

active <- primary[
  primary$supplement_group %in%
    c("Other supplements only", "Omega-3 supplement"),
  ,
  drop = FALSE
]
active$active_omega3 <- as.integer(
  active$supplement_group == "Omega-3 supplement"
)
active_formula <- as.formula(
  paste(
    "global_cognition_complete ~ active_omega3 +",
    adjustment_terms
  )
)
active_model <- survey_lm_taylor(
  active_formula, active, "WTMEC4YR", stratum_prefix = "cycle"
)
active_result <- tidy_term(
  active_model,
  "active_omega3",
  "Omega-3 supplement vs other supplements only",
  reference = "t"
)

active_counts <- aggregate(
  SEQN ~ supplement_group, data = active, FUN = length
)
names(active_counts)[2] <- "unweighted_n"
active_weight_totals <- aggregate(
  WTMEC4YR ~ supplement_group, data = active, FUN = sum
)
active_weight_totals$weighted_percent <- with(
  active_weight_totals, 100 * WTMEC4YR / sum(WTMEC4YR)
)
active_counts <- merge(
  active_counts,
  active_weight_totals[c("supplement_group", "weighted_percent")],
  by = "supplement_group",
  sort = FALSE
)

# Analysis 2: cycle-specific replication.
cycle_results <- list()
for (cycle_level in levels(primary$cycle)) {
  cycle_data <- primary[primary$cycle == cycle_level, , drop = FALSE]
  cycle_model <- survey_lm_taylor(
    primary_formula, cycle_data, "WTMEC2YR"
  )
  cycle_results[[cycle_level]] <- tidy_term(
    cycle_model,
    "omega_user_30d",
    paste0("Any omega-3 supplement vs non-use: ", cycle_level),
    reference = "normal"
  )
}
cycle_results <- do.call(rbind, cycle_results)
rownames(cycle_results) <- NULL

interaction_formula <- as.formula(
  paste(
    "global_cognition_complete ~ omega_user_30d * cycle +",
    adjustment_terms
  )
)
interaction_model <- survey_lm_taylor(
  interaction_formula, primary, "WTMEC4YR", stratum_prefix = "cycle"
)
interaction_term <- grep(
  "^omega_user_30d:cycle", names(interaction_model$coefficients), value = TRUE
)
if (length(interaction_term) != 1) {
  stop("Could not uniquely identify the omega-3-by-cycle interaction term.")
}
interaction_result <- tidy_term(
  interaction_model,
  interaction_term,
  "Omega-3 use by survey-cycle interaction",
  reference = "t"
)

all_results <- rbind(active_result, cycle_results, interaction_result)

write.csv(
  all_results,
  file.path(output_dir, "secondary_analysis_results.csv"),
  row.names = FALSE
)
write.csv(
  active_counts,
  file.path(output_dir, "active_comparator_counts.csv"),
  row.names = FALSE
)

capture.output(
  {
    cat("PRESPECIFIED SECONDARY ANALYSES — 30 JULY 2026\n\n")
    cat("Primary-model verification:\n")
    cat(
      sprintf(
        "beta = %.15f; SE = %.15f; exact verification passed\n\n",
        frozen_primary$coefficients[["omega_user_30d"]],
        frozen_primary$standard_errors[["omega_user_30d"]]
      )
    )
    cat("Active-comparator group sizes:\n")
    print(active_counts, row.names = FALSE)
    cat("\nSpecified secondary results:\n")
    print(all_results, row.names = FALSE, digits = 6)
    cat("\nR session:\n")
    print(sessionInfo())
  },
  file = file.path(output_dir, "secondary_analysis_report.txt")
)

message(
  "Secondary analyses complete: ",
  normalizePath(output_dir, mustWork = TRUE)
)
