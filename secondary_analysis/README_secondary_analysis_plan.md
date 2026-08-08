**Project:** Omega-3 supplementation, dietary EPA+DHA, and cognition in NHANES 2011–2014  
**Plan date:** 30 July 2026 
**Status at signing:** Written after all primary and previously reported sensitivity models were frozen and before inspecting the results of the analyses below.

## Scope
1. An active-comparator analysis of omega-3 supplement users versus users of other supplements.
2. Cycle-specific replication in NHANES 2011–2012 and 2013–2014, plus a formal omega-3-use-by-cycle interaction test.


## Frozen primary definitions

- **Population:** NHANES participants aged 60 years or older.
- **Outcome:** The frozen complete four-component global cognition z-score, formed from CERAD immediate learning, CERAD delayed recall, Animal Fluency, and DSST. Component standardization remains based on the combined 2011–2014 complete-battery sample.
- **Omega-3 exposure:** Any identified omega-3-containing supplement reported in the dedicated past-30-day dietary-supplement interview.
- **Core adjustment set:** age, sex, race/ethnicity, education, family income-to-poverty ratio, BMI, diagnosed diabetes, and PHQ-9 score.
- **Design variables:** NHANES strata, PSUs, and the MEC examination weight. Combined-cycle analyses use WTMEC4YR (WTMEC2YR/2); separate-cycle analyses use WTMEC2YR. Scaling the weights by a cycle-constant does not change coefficients or design-based standard errors.
- **Missing data:** complete-case analysis for the outcome, exposure, core adjustment set, weight, strata, and PSU variables.

## Analysis 1: active comparator

The primary active-comparator model will be restricted to:

- **Omega-3 users:** at least one identified omega-3-containing product in the past-30-day supplement file.
- **Other-supplement users:** reported any dietary-supplement use in the past 30 days (DSD010 = 1) but had no identified omega-3-containing product.

Participants reporting no dietary-supplement use, an ambiguous/refused response, or missing supplement-use status will not enter this comparison.

A survey-weighted linear regression will estimate the adjusted mean difference in global cognition for omega-3 users relative to other-supplement users, using the frozen adjustment set and combined-cycle survey design. The coefficient, design-based standard error, 95% confidence interval, p-value, unweighted group sizes, and survey-weighted group percentages will be retained.

The active-comparator result will replace the stroke-exclusion row in the main poster forest plot if the model is estimable and both comparison groups contain at least 100 unweighted participants. Statistical significance is not a criterion for inclusion.

## Analysis 2: cycle-specific replication and interaction

The frozen primary supplement model (omega-3 users versus non-users of omega-3, regardless of other supplement use) will be fitted separately in:

- NHANES 2011–2012.
- NHANES 2013–2014.

Each model will use its cycle-specific MEC weight and the same outcome and adjustment set. Cycle-specific coefficients and Taylor-linearized design-based standard errors will be reported. Because the fully adjusted single-cycle models exhaust the default residual survey degrees of freedom, their 95% confidence intervals and descriptive p-values will use the standard Normal reference distribution. These cycle-specific p-values will not be used to declare replication or heterogeneity.

Formal heterogeneity will be tested in the pooled 2011–2014 survey design by adding survey cycle, omega-3 use, and their interaction to the frozen primary model. The omega-3-by-cycle interaction coefficient and its design-based Wald p-value, using the pooled model's residual survey degrees of freedom, will be the inferential test of between-cycle difference.

Interpretation will focus on the direction, magnitude, confidence-interval overlap, and interaction test—not on whether one cycle crosses p = 0.05 and the other does not.

## Reporting commitment

All specified estimates will be kept and reported whether statistically significant or not. Results will be labelled as prespecified secondary analyses dated 30 July 2026 and will not be described as confirmatory evidence of causality.

