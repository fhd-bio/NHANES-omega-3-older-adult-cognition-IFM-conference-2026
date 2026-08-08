# IFM 2026 Poster: Omega-3 Exposure and Cognition in Older Adults

This repository contains the analysis code, outputs, and supporting documentation for an IFM 2026 poster examining the association of omega-3 supplement use and dietary EPA+DHA intake with cognitive performance among US adults aged 60 years or older.

The analysis uses data from the 2011–2012 and 2013–2014 cycles of the National Health and Nutrition Examination Survey (NHANES).

## Study overview

The study evaluates two omega-3 exposures:

1. Past-30-day use of identified omega-3 supplements.
2. Day-one dietary intake of eicosapentaenoic acid and docosahexaenoic acid, expressed as combined EPA+DHA intake.

The outcome is a global cognition z-score derived from four cognitive measures:

- CERAD immediate learning;
- CERAD delayed recall;
- Animal Fluency;
- Digit Symbol Substitution Test.

The primary cognitive outcome requires complete data for all four components.

All regression models account for the complex NHANES survey design and are adjusted for:

- age;
- sex;
- race and ethnicity;
- educational attainment;
- family income-to-poverty ratio;
- body mass index;
- diagnosed diabetes;
- PHQ-9 score.

## Repository structure

### `01_main_analysis`

- `omega3_nhanes_2011_2014_corrected.R`  
  Rebuilds the analytic variables from the raw NHANES XPT files and runs the primary supplement, dietary, and sensitivity analyses.

- `omega3_reanalysis.py`  
  Provides an independent Python reconstruction and audit of the main analyses. It uses the past-30-day supplement exposure, requires a complete cognitive battery, and applies the Day-one dietary recall weight to dietary EPA+DHA models.

### `02_secondary_analysis`

- `secondary_analysis_plan_2026-07-30.md`  
  Dated analysis plan for the active-comparator and cycle-consistency analyses.

- `run_secondary_analyses.R`  
  Runs the active-comparator, cycle-specific, omega-3-by-cycle interaction, and exploratory formulation analyses.

### `03_outputs`

- `r_key_models.csv`  
  Estimates from the primary supplement, dietary, and sensitivity models.

- `r_sample_flow.csv`  
  Participant counts used in the poster sample-flow figure.

- `baseline_characteristics_poster.tsv`  
  Survey-weighted participant characteristics used on the poster.

- `active_comparator_counts.csv`  
  Group counts for the active-comparator analysis.

- `secondary_analysis_results.csv`  
  Full results from the secondary analyses.

- `secondary_analysis_summary.md`  
  Concise interpretation and poster-reporting summary.

- `python_key_results.json`  
  Machine-readable results from the Python reanalysis for comparison with the R outputs.

- `formulation_counts.csv`  
  Unweighted subgroup sizes for the exploratory omega-3 formulation analysis.

- `formulation_exploratory_results.csv`  
  Adjusted formulation contrasts and the overall survey-design Wald test.

## Analytic sample

A total of 3,472 adults aged 60 years or older were assessed for eligibility.

For the supplement analysis:

- 881 participants were excluded because of incomplete cognition data, missing supplement exposure, or missing core covariates;
- 2,591 participants remained in the primary supplement model.

For the dietary analysis:

- a further 138 participants lacked a reliable Day-one dietary recall or other required dietary-model data;
- 2,453 participants remained in the dietary EPA+DHA model.

## Exposure definitions

### Omega-3 supplement use

Supplement exposure is based on identified omega-3-containing products reported during the dedicated past-30-day dietary supplement interview.

Participants are classified as:

- omega-3 supplement users;
- omega-3 non-users.

### Dietary EPA+DHA

Dietary exposure is based on combined EPA and DHA intake from the Day-one 24-hour dietary recall.

Because this measure reflects intake from a single recall day, it is treated as a one-day intake estimate rather than a measure of habitual long-term intake.

The supplement and dietary coefficients represent different exposure definitions and should not be directly compared as though they were measured on the same scale.

## Survey design and weighting

All analyses incorporate the NHANES complex survey design, including:

- sampling weights;
- strata;
- primary sampling units.

The pooled R analyses divide the two-year weights by two to create combined four-year weights.

- Supplement models use the MEC examination weight: `WTMEC2YR / 2`.
- Dietary EPA+DHA models use the Day-one dietary recall weight: `WTDRD1 / 2`.

The Python audit uses `WTDRD1` directly for the dietary model. Dividing every survey weight by the same constant changes only the common scale of the weights and does not change the fitted regression estimates.

## Main results

### Omega-3 supplement use

Past-30-day omega-3 supplement users had higher adjusted global cognition scores than non-users:

- β = 0.151;
- 95% CI: 0.076 to 0.226;
- p < 0.001.

### Dietary EPA+DHA

Higher Day-one dietary EPA+DHA intake was positively associated with global cognition:

- β = 0.217;
- 95% CI: 0.031 to 0.403;
- p = 0.025.

## Secondary analyses

### Active-comparator analysis

To reduce, although not eliminate, the influence of general supplement-taking behaviour, 499 omega-3 supplement users were compared with 1,211 users of non-omega-3 supplements.

The association remained positive:

- β = 0.149;
- 95% CI: 0.064 to 0.234;
- p = 0.002.

### Cycle consistency

The supplement association was examined separately in each NHANES cycle:

- 2011–2012: β = 0.157;
- 2013–2014: β = 0.139.

The pooled omega-3-by-cycle interaction test showed no evidence that the association differed between cycles:

- p for interaction = 0.662.

### Exploratory formulation analysis

Omega-3 formulation was explored among the 499 supplement users.

Fish oil or cod-liver oil was used as the reference category, while very rare products were pooled into an other/mixed group.

The overall formulation comparison was not statistically significant:

- Wald F(3,14) = 0.72;
- p = 0.557.

The non-fish-oil formulation groups were small and produced wide confidence intervals. No formulation-specific conclusions are therefore drawn.

## Baseline characteristics

The poster baseline table compares:

- 499 omega-3 supplement users;
- 2,092 omega-3 non-users.

It reports survey-weighted means or percentages together with absolute standardized mean differences.

Differences in education, income, smoking, diabetes, and physical activity indicate that residual confounding and healthy-user bias remain plausible despite multivariable adjustment.

## Interpretation

The findings show a positive cross-sectional association between omega-3 exposure and global cognition in older US adults.

They do not establish that omega-3 supplements or dietary EPA+DHA caused better cognitive performance.

The observational design, single-day dietary assessment, potential exposure misclassification, residual confounding, and healthy-user bias limit causal interpretation.

The results should not be used to recommend omega-3 supplementation solely for cognitive protection.

## Context with recent evidence

The NHANES findings are interpreted alongside recent evidence from:

- the KLOSCAD longitudinal cohort;
- the ADNI cohort;
- the PreventE4 randomized trial.

These studies do not provide a consistent causal picture. The poster therefore presents the NHANES findings as observational evidence that requires cautious interpretation.

## Reproduction

Raw NHANES XPT files are not included in this repository because they are publicly available and would substantially increase repository size.

The R and Python scripts rebuild the analysis using the required 2011–2012 (`_G`) and 2013–2014 (`_H`) NHANES component files.

Local input and output directory paths may need to be edited before running the scripts.

## Important note

This repository supports a conference poster and should not be interpreted as a clinical guideline or treatment recommendation.
All analysis code in this repository was generated with the assistance of generative AI tools. The study design, analytic decisions, variable definitions, interpretation of results, and final reporting were determined and reviewed by the authors. The generated code and statistical outputs were checked and validated by the authors, who remain fully responsible for the accuracy and integrity of the work.
