# Validation record

## Status

**Passed on 8 August 2026.**

The public repository was cloned into a fresh Google Colab session and the Python reconstruction, independent main R pipeline, and prespecified secondary R analyses were run from their documented inputs. All three workflows completed successfully and reproduced the committed results to numerical precision.

## Repository version tested

- Repository: [`fhd-bio/NHANES-omega-3-older-adult-cognition-IFM-conference-2026`](https://github.com/fhd-bio/NHANES-omega-3-older-adult-cognition-IFM-conference-2026)
- Branch: `main`
- Commit: [`a43aafafbf105076fe4c58420a479f713d4d1769`](https://github.com/fhd-bio/NHANES-omega-3-older-adult-cognition-IFM-conference-2026/commit/a43aafafbf105076fe4c58420a479f713d4d1769)
- Executed notebook: [`Validation_check.ipynb`](Validation_check.ipynb)

## Validation procedure

1. Cloned the public repository into a fresh Google Colab session.
2. Installed the Python and R dependencies documented by the repository.
3. Downloaded all 26 required public NHANES 2011–2012 and 2013–2014 XPT component files.
4. Confirmed that 26 XPT files were present and that none was smaller than 10 KB.
5. Ran the independent Python reconstruction:

   ```bash
   python3 main_analysis/omega3_analysis.py
   ```

6. Ran the independent main R pipeline from the raw XPT files:

   ```bash
   Rscript main_analysis/omega3_nhanes_2011_2014_analysis.R \
     tmp/raw_nhanes \
     tmp/r_output
   ```

7. Ran the prespecified secondary R analyses using the analytic CSV produced by Python:

   ```bash
   Rscript secondary_analysis/secondary_analysis.R \
     tmp/reanalysis/rebuilt_analytic_data.csv \
     tmp/secondary_output
   ```

8. Compared the clean-run samples and estimates with the committed outputs and the values reported on the poster.

## Reproduced primary results

| Analysis | n | Estimate | 95% CI | p-value |
|---|---:|---:|---:|---:|
| Past-30-day omega-3 supplement use | 2,591 | 0.150966 | 0.076370 to 0.225562 | 0.000562 |
| Day-one dietary EPA+DHA | 2,453 | 0.216779 | 0.030940 to 0.402618 | 0.025002 |

The independent R and Python estimates agreed to floating-point precision and reproduce the rounded poster values of β = 0.151 and β = 0.217.

## Reproduced secondary results

| Analysis | n | Estimate | 95% CI | p-value |
|---|---:|---:|---:|---:|
| Omega-3 users versus users of other supplements | 1,710 | 0.148892 | 0.063613 to 0.234171 | 0.001937 |
| Supplement association, 2011–2012 | 1,191 | 0.157450 | 0.038218 to 0.276682 | 0.009648 |
| Supplement association, 2013–2014 | 1,400 | 0.139452 | 0.071031 to 0.207872 | 0.0000648 |
| Omega-3-by-cycle interaction | 2,591 | -0.031831 | -0.184459 to 0.120797 | 0.661500 |

The active-comparator sample contained 499 omega-3 supplement users and 1,211 users of other supplements. The secondary script's internal benchmark check also confirmed the primary supplement estimate and standard error exactly before fitting the secondary models.

## Output comparison

The clean-run outputs were compared with the corresponding files in `outputs/`:

- sample-flow counts and active-comparator counts were byte-for-byte identical;
- all corresponding analyses, groups, sample sizes, estimates, standard errors, confidence intervals, test statistics, and p-values agreed;
- where files were not byte-for-byte identical, the differences were limited to CSV quoting or floating-point representation at approximately the 13th decimal place or later;
- no obsolete lifestyle/alcohol sensitivity model appeared in the clean-run outputs.

These negligible machine-precision differences do not affect any rounded value, inference, interpretation, table, or figure.

## Validation environment

- Google Colab, Python 3 kernel
- Python 3.12
- R 4.6.0 on Ubuntu 22.04.5 LTS
- `haven` 2.5.5
- `survey` 4.5
- `dplyr` 1.2.1

Complete R session information is recorded in the generated `r_analysis_report.txt` and `secondary_analysis_report.txt` files. The notebook contains the executed commands and retained outputs. The package installer reported Colab dependency-conflict warnings after successfully installing the pinned Python packages; these warnings did not prevent any analysis from running or alter the reproduced results.

## Data and repository notes

The 26 raw NHANES XPT files were downloaded directly from the public CDC/NCHS source for validation and are not redistributed in this repository. Participant-level derived data, the fitted-model RDS file, and the complete validation ZIP are also not committed.

This validation confirms reproducible execution of the tested repository revision and agreement of the generated results with the committed outputs. It does not convert the cross-sectional findings into causal evidence or remove the limitations described in the poster and README.
