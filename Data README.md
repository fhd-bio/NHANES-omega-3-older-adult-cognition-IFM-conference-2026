# Required NHANES data files

The raw NHANES SAS transport files are public and are not committed to this repository.

Place all files in:

```text
tmp/raw_nhanes/
```

The analyses require 13 components from each of two cycles, for 26 XPT files in total.

| Component | 2011–2012 | 2013–2014 | Purpose |
|---|---|---|---|
| Demographics and weights | `DEMO_G.XPT` | `DEMO_H.XPT` | Age, demographic covariates, survey design, MEC weights |
| Cognitive functioning | `CFQ_G.XPT` | `CFQ_H.XPT` | CERAD, Animal Fluency, and DSST outcomes |
| Day-one dietary totals | `DR1TOT_G.XPT` | `DR1TOT_H.XPT` | EPA, DHA, energy, dietary recall status and weight |
| Body measures | `BMX_G.XPT` | `BMX_H.XPT` | BMI |
| Diabetes | `DIQ_G.XPT` | `DIQ_H.XPT` | Diagnosed diabetes |
| Depression screener | `DPQ_G.XPT` | `DPQ_H.XPT` | PHQ-9 items |
| Medical conditions | `MCQ_G.XPT` | `MCQ_H.XPT` | Stroke variable used in the retained sensitivity analysis |
| Smoking | `SMQ_G.XPT` | `SMQ_H.XPT` | Smoking characteristic for the baseline table |
| Physical activity | `PAQ_G.XPT` | `PAQ_H.XPT` | Recreational activity characteristic for the baseline table |
| Day-one supplement products | `DS1IDS_G.XPT` | `DS1IDS_H.XPT` | Earlier recall-based exposure replication |
| Day-two supplement products | `DS2IDS_G.XPT` | `DS2IDS_H.XPT` | Earlier recall-based exposure replication |
| Past-30-day supplement products | `DSQIDS_G.XPT` | `DSQIDS_H.XPT` | Primary omega-3 supplement exposure and formulation |
| Past-30-day supplement totals | `DSQTOT_G.XPT` | `DSQTOT_H.XPT` | Any-supplement indicator for the active comparator |

## Official sources

- [NHANES 2011–2012 data, documentation, and codebooks](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2011)
- [NHANES 2013–2014 data, documentation, and codebooks](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2013)
- [CDC NHANES dataset and documentation tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/datasets.aspx)

## Copy-and-paste download commands

Run these commands from the repository root. They download the exact public XPT files from the CDC/NCHS servers.

```bash
mkdir -p tmp/raw_nhanes

for stem in DEMO CFQ DR1TOT BMX DIQ DPQ MCQ SMQ PAQ DS1IDS DS2IDS DSQIDS DSQTOT
do
  curl -fL "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/${stem}_G.xpt" \
    -o "tmp/raw_nhanes/${stem}_G.XPT"
  curl -fL "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/${stem}_H.xpt" \
    -o "tmp/raw_nhanes/${stem}_H.XPT"
done
```

Confirm that 26 files were downloaded:

```bash
find tmp/raw_nhanes -maxdepth 1 -type f -iname '*.xpt' | wc -l
```

The expected result is:

```text
26
```

## File-integrity checks

A failed download can sometimes leave an HTML error page with an `.XPT` filename. Before running the analyses:

```bash
find tmp/raw_nhanes -maxdepth 1 -type f -iname '*.xpt' -size -10k -print
```

This command should return no files. The R and Python scripts will also fail if a required component cannot be parsed.

## Naming conventions

The `_G` suffix denotes NHANES 2011–2012, and `_H` denotes NHANES 2013–2014. The R script accepts upper- or lower-case `.xpt` extensions. The Python script expects the filenames shown in the table.

