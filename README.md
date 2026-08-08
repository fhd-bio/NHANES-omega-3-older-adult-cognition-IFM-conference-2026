# IFM 2026 Poster: Omega-3 Exposure and Cognition in Older Adults

This repository contains the analysis code, outputs, and supporting documentation for an IFM 2026 poster examining associations of omega-3 supplement use and dietary EPA+DHA intake with cognitive performance among US adults aged 60 years or older.

The analysis uses the 2011–2012 and 2013–2014 cycles of the National Health and Nutrition Examination Survey (NHANES) and has a cross-sectional design.

## Study overview

The study evaluates two omega-3 exposures:

- Past-30-day use of identified omega-3 supplements.
- Day-one dietary intake of eicosapentaenoic acid and docosahexaenoic acid, expressed as combined EPA+DHA intake.

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

### `main_analysis/`

- [`omega3_nhanes_2011_2014_analysis.R`](main_analysis/omega3_nhanes_2011_2014_analysis.R) — reconstructs the analytic variables from the raw NHANES XPT files and runs the main R analyses.
- [`omega3_analysis.py`](main_analysis/omega3_analysis.py) — provides an independent Python reconstruction and audit of the main analyses. It uses the past-30-day supplement exposure, requires a complete cognitive battery, applies the Day-one dietary recall weight to dietary EPA+DHA models, and writes the analytic CSV required by the secondary R script.

### `secondary_analysis/`

- [`README_secondary_analysis_plan.md`](secondary_analysis/README_secondary_analysis_plan.md) — dated analysis plan for the active-comparator and cycle-consistency analyses.
- [`secondary_analysis.R`](secondary_analysis/secondary_analysis.R) — runs the active-comparator, cycle-specific, and omega-3-by-cycle interaction analyses. This script must be run after `main_analysis/omega3_analysis.py` because it reads the analytic CSV created by the Python reconstruction.

### `outputs/`

- [`r_key_models.csv`](outputs/r_key_models.csv) — estimates from the main R analyses.
- [`r_sample_flow.csv`](outputs/r_sample_flow.csv) — participant counts used in the poster sample-flow figure.
- [`baseline_characteristics_poster.tsv`](outputs/baseline_characteristics_poster.tsv) — survey-weighted participant characteristics used on the poster.
- [`active_comparator_counts.csv`](outputs/active_comparator_counts.csv) — group counts for the active-comparator analysis.
- [`secondary_analysis_results.csv`](outputs/secondary_analysis_results.csv) — results from the active-comparator and cycle-consistency analyses.
- [`python_key_results.json`](outputs/python_key_results.json) — machine-readable results from the independent Python reconstruction.
- [`formulation_counts.csv`](outputs/formulation_counts.csv) — unweighted subgroup sizes for the exploratory formulation analysis.
- [`formulation_specific_exploratory_results.csv`](outputs/formulation_specific_exploratory_results.csv) — adjusted exploratory formulation contrasts and the overall survey-design Wald test.

### Other files

- [`requirements.txt`](requirements.txt) — pinned Python dependencies.
- [`Data README.md`](Data%20README.md) — required NHANES files and official download links.
- [`R_DEPENDENCIES.md`](R_DEPENDENCIES.md) — R dependency and installation instructions.
- `VALIDATION.md` — validation status, environment, checks, and final clean-run record; to be added after the final clean R run.

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

All analyses incorporate the NHANES complex survey design, including sampling weights, strata, and primary sampling units.

The pooled R analyses divide the two-year weights by two to create combined four-year weights:

- supplement models use `WTMEC2YR / 2`;
- dietary EPA+DHA models use `WTDRD1 / 2`.

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

Omega-3 formulation was explored among the 499 supplement users. Fish oil or cod-liver oil was used as the reference category, while rare products were pooled into an other/mixed group.

The overall formulation comparison was not statistically significant:

- Wald F(3,14) = 0.72;
- p = 0.557.

The non-fish-oil formulation groups were small and produced wide confidence intervals. No formulation-specific conclusions are drawn.

## Baseline characteristics

The poster baseline table compares 499 omega-3 supplement users with 2,092 non-users and reports survey-weighted means or percentages with absolute standardized mean differences.

Differences in education, income, smoking, diabetes, and physical activity indicate that residual confounding and healthy-user bias remain plausible despite multivariable adjustment.

## Interpretation

The findings show positive cross-sectional associations between omega-3 exposure and global cognition in older US adults.

They do not establish that omega-3 supplements or dietary EPA+DHA caused better cognitive performance. The cross-sectional design, single-day dietary assessment, possible exposure misclassification, residual confounding, and healthy-user bias limit causal interpretation.

The results should not be used to recommend omega-3 supplementation solely for cognitive protection.

## Context with recent evidence

The NHANES findings are interpreted alongside recent evidence from the KLOSCAD longitudinal cohort, ADNI cohort, and PreventE4 randomized trial. These studies do not provide a consistent causal picture, so the poster presents the NHANES findings as observational evidence requiring cautious interpretation.

## Dependencies

### Python

Python 3.12 is recommended. From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

### R

The main R pipeline requires `haven`, `survey`, and `dplyr`. The secondary R script uses base R only. See [`R_DEPENDENCIES.md`](R_DEPENDENCIES.md).

## Data preparation

Raw NHANES XPT files are not committed because they are publicly available and would substantially increase repository size.

Download the 26 required files into `tmp/raw_nhanes/` as described in [`Data README.md`](Data%20README.md).

## Run the analysis

Run all commands from the repository root.

### 1. Create the Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

### 2. Reconstruct the data and run the independent Python analysis

```bash
python3 main_analysis/omega3_analysis.py
```

This writes results to `tmp/reanalysis/`, including:

- `tmp/reanalysis/key_results.json`;
- `tmp/reanalysis/rebuilt_analytic_data.csv`.

### 3. Run the secondary R analyses

This step must follow Step 2 because it consumes the Python-created analytic CSV.

```bash
Rscript secondary_analysis/secondary_analysis.R \
  tmp/reanalysis/rebuilt_analytic_data.csv \
  tmp/secondary_output
```

### 4. Run the independent main R pipeline

```bash
Rscript main_analysis/omega3_nhanes_2011_2014_analysis.R \
  tmp/raw_nhanes \
  tmp/r_output
```

The main R and Python pipelines each reconstruct the analytic variables directly from the original public NHANES component files. In this repository, “reconstruct” means merging participant-level component files, deriving the exposures, outcome, covariates, survey weights, and analytic samples, and then refitting the reported models.

### 5. Review validation outputs

Compare the fresh outputs with the committed files in `outputs/` and complete the final clean-run section in `VALIDATION.md` after that file is added.

## Important note

This repository supports a conference poster and should not be interpreted as a clinical guideline or treatment recommendation.

All analysis code in this repository was generated with the assistance of generative AI tools. The study design, analytic decisions, variable definitions, interpretation of results, and final reporting were determined and reviewed by the authors. The authors remain fully responsible for the accuracy and integrity of the work.
