"""Python implementation of the NHANES 2011-2014 omega-3 cognition analysis."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import f as f_dist
from scipy.stats import t as t_dist


RAW = Path("tmp/raw_nhanes")
OUT = Path("tmp/reanalysis")
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = ("G", "H")
CONTINUOUS_COVARIATES = ["RIDAGEYR", "INDFMPIR", "BMXBMI", "phq9_total"]
CATEGORICAL_COVARIATES = ["RIAGENDR", "RIDRETH1", "DMDEDUC2", "DIQ010"]
COGNITIVE_COMPONENTS = ["cerad_total", "CFDCSR", "CFDAST", "CFDDS"]

OMEGA_PATTERN = re.compile(
    r"fish oil|salmon oil|marine oil|cod liver|krill|flax|flaxseed|linseed|"
    r"algae|algal|vegetarian dha|vegan omega|lovaza|vascepa|epanova|"
    r"omega ?-?3",
    re.I,
)


def read_xpt(stem: str, cycle: str) -> pd.DataFrame:
    return pd.read_sas(
        RAW / f"{stem}_{cycle}.XPT", format="xport", encoding="latin1"
    )


def clean_special_codes(series: pd.Series, codes: tuple[int, ...]) -> pd.Series:
    return series.replace({code: np.nan for code in codes})


def omega_ids(frame: pd.DataFrame) -> set[float]:
    names = frame["DSDSUPP"].fillna("").astype(str)
    return set(frame.loc[names.str.contains(OMEGA_PATTERN), "SEQN"])


def formulation(name: object) -> str | None:
    if pd.isna(name):
        return None
    text = str(name)
    rules = [
        ("prescription", r"lovaza|vascepa|epanova|omega ?-?3 ?-?acid ?ethyl"),
        ("krill", r"krill"),
        ("algal", r"algae|algal|vegetarian dha|vegan omega"),
        ("cod_liver", r"cod liver"),
        ("plant", r"flax|flaxseed|linseed"),
        ("fish_oil", r"fish oil|salmon oil|marine oil|omega ?-?3"),
    ]
    for label, pattern in rules:
        if re.search(pattern, text, re.I):
            return label
    return None


def collapsed_formulation(user: int, labels: object) -> str:
    """Collapse product labels into estimable, poster-relevant groups."""
    if int(user) == 0:
        return "Non-user"
    parsed = set(str(labels).split(",")) if pd.notna(labels) else set()
    if parsed and parsed.issubset({"fish_oil", "cod_liver"}):
        return "Fish oil/cod liver oil"
    if parsed == {"krill"}:
        return "Krill"
    if parsed == {"plant"}:
        return "Plant-based"
    if parsed:
        return "Other/mixed"
    return "Unclassified"


def build_cycle(cycle: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    demo = read_xpt("DEMO", cycle)
    cfq = read_xpt("CFQ", cycle)
    diet = read_xpt("DR1TOT", cycle)
    bmx = read_xpt("BMX", cycle)
    diq = read_xpt("DIQ", cycle)
    dpq = read_xpt("DPQ", cycle)
    mcq = read_xpt("MCQ", cycle)
    smq = read_xpt("SMQ", cycle)
    alq = read_xpt("ALQ", cycle)
    paq = read_xpt("PAQ", cycle)
    ds1 = read_xpt("DS1IDS", cycle)
    ds2 = read_xpt("DS2IDS", cycle)
    ds30 = read_xpt("DSQIDS", cycle)
    ds30_total = read_xpt("DSQTOT", cycle)

    day1_day2 = pd.concat(
        [ds1[["SEQN", "DSDSUPP"]], ds2[["SEQN", "DSDSUPP"]]],
        ignore_index=True,
    )
    ids_24h = omega_ids(day1_day2)
    ids_30d = omega_ids(ds30)

    products = ds30.loc[
        ds30["SEQN"].isin(ids_30d), ["SEQN", "DSDSUPP"]
    ].copy()
    products["formulation"] = products["DSDSUPP"].map(formulation)
    products = products.dropna(subset=["formulation"])
    form_by_person = (
        products.groupby("SEQN")["formulation"]
        .agg(lambda x: ",".join(sorted(set(x))))
        .rename("omega_formulations_30d")
        .reset_index()
    )

    cognition = cfq[
        ["SEQN", "CFDCST1", "CFDCST2", "CFDCST3", "CFDCSR", "CFDAST", "CFDDS"]
    ].copy()
    cognition["cerad_total"] = cognition[
        ["CFDCST1", "CFDCST2", "CFDCST3"]
    ].sum(axis=1, min_count=3)

    # Reproduce the accepted outcome exactly for auditing only.
    cognition["cerad_total_partial"] = cognition[
        ["CFDCST1", "CFDCST2", "CFDCST3"]
    ].sum(axis=1, min_count=1)
    accepted_parts = ["cerad_total_partial", "CFDCSR", "CFDAST", "CFDDS"]
    for variable in accepted_parts:
        cognition[f"{variable}_z_accepted"] = (
            cognition[variable] - cognition[variable].mean()
        ) / cognition[variable].std()
    cognition["global_cognition_accepted"] = cognition[
        [f"{variable}_z_accepted" for variable in accepted_parts]
    ].mean(axis=1)

    dpq_items = [f"DPQ0{i}0" for i in range(1, 10)]
    depression = dpq[["SEQN"] + dpq_items].copy()
    answered = depression[dpq_items].replace({7: np.nan, 9: np.nan})
    count = answered.notna().sum(axis=1)
    depression["phq9_total"] = np.where(
        count >= 7,
        answered.sum(axis=1, min_count=1) * 9 / count,
        np.nan,
    )

    smoking = smq[["SEQN", "SMQ020", "SMQ040"]].copy()
    smoking["smoking_status"] = np.select(
        [
            smoking["SMQ020"].eq(2),
            smoking["SMQ020"].eq(1) & smoking["SMQ040"].eq(3),
            smoking["SMQ020"].eq(1) & smoking["SMQ040"].isin([1, 2]),
        ],
        ["never", "former", "current"],
        default=None,
    )

    alcohol = alq[["SEQN", "ALQ101", "ALQ120Q"]].copy()
    alcohol["current_alcohol"] = np.select(
        [
            alcohol["ALQ101"].eq(2) | alcohol["ALQ120Q"].eq(0),
            alcohol["ALQ120Q"].between(1, 366),
        ],
        [0.0, 1.0],
        default=np.nan,
    )

    activity = paq[["SEQN", "PAQ650", "PAQ665"]].copy()
    activity["recreational_activity"] = np.select(
        [
            activity["PAQ650"].eq(1) | activity["PAQ665"].eq(1),
            activity["PAQ650"].eq(2) & activity["PAQ665"].eq(2),
        ],
        [1.0, 0.0],
        default=np.nan,
    )

    selected_demo = [
        "SEQN",
        "RIDAGEYR",
        "RIAGENDR",
        "RIDRETH1",
        "DMDEDUC2",
        "INDFMPIR",
        "SDMVPSU",
        "SDMVSTRA",
        "WTMEC2YR",
    ]
    selected_diet = [
        "SEQN",
        "WTDRD1",
        "DR1DRSTZ",
        "DR1TKCAL",
        "DR1TP205",
        "DR1TP226",
    ]
    data = demo[selected_demo].merge(cognition, on="SEQN", how="inner")
    for extra in [
        diet[selected_diet],
        bmx[["SEQN", "BMXBMI"]],
        diq[["SEQN", "DIQ010"]],
        depression[["SEQN", "phq9_total"]],
        mcq[["SEQN", "MCQ160F"]],
        smoking[["SEQN", "smoking_status"]],
        alcohol[["SEQN", "current_alcohol"]],
        activity[["SEQN", "recreational_activity"]],
        ds30_total[["SEQN", "DSD010"]],
        form_by_person,
    ]:
        data = data.merge(extra, on="SEQN", how="left")

    data = data.loc[data["RIDAGEYR"] >= 60].copy()
    data["cycle"] = cycle
    data["omega_user_24h"] = data["SEQN"].isin(ids_24h).astype(int)
    data["omega_user_30d"] = data["SEQN"].isin(ids_30d).astype(int)
    data["dietary_epa_dha_g"] = data["DR1TP205"] + data["DR1TP226"]
    data["log_dietary_epa_dha"] = np.log1p(data["dietary_epa_dha_g"])

    data["DMDEDUC2"] = clean_special_codes(data["DMDEDUC2"], (7, 9))
    data["DIQ010"] = clean_special_codes(data["DIQ010"], (7, 9))
    data["MCQ160F"] = clean_special_codes(data["MCQ160F"], (7, 9))
    return data, products


def weighted_quantile(
    values: np.ndarray, quantiles: np.ndarray, weights: np.ndarray
) -> np.ndarray:
    order = np.argsort(values)
    values = np.asarray(values)[order]
    weights = np.asarray(weights)[order]
    cumulative = np.cumsum(weights) - 0.5 * weights
    cumulative = cumulative / weights.sum()
    return np.interp(quantiles, cumulative, values)


def rcs_basis(x: np.ndarray, knots: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    knots = np.asarray(knots, dtype=float)
    if len(knots) < 3 or np.any(np.diff(knots) <= 0):
        raise ValueError(f"Restricted cubic spline knots must be unique: {knots}")
    last = knots[-1]
    penultimate = knots[-2]
    scale = (last - knots[0]) ** 2
    columns = [x]
    for knot in knots[:-2]:
        term = np.maximum(x - knot, 0) ** 3
        term -= (
            (last - knot) / (last - penultimate)
        ) * np.maximum(x - penultimate, 0) ** 3
        term += (
            (penultimate - knot) / (last - penultimate)
        ) * np.maximum(x - last, 0) ** 3
        columns.append(term / scale)
    return np.column_stack(columns)


def make_design(
    data: pd.DataFrame,
    continuous: list[str],
    categorical: list[str],
    matrices: dict[str, np.ndarray] | None = None,
) -> tuple[np.ndarray, list[str]]:
    columns = [np.ones(len(data))]
    names = ["Intercept"]
    for variable in continuous:
        columns.append(data[variable].to_numpy(dtype=float))
        names.append(variable)
    for variable in categorical:
        dummies = pd.get_dummies(
            data[variable].astype(str),
            prefix=variable,
            drop_first=True,
            dtype=float,
        )
        for name in dummies:
            columns.append(dummies[name].to_numpy())
            names.append(name)
    for prefix, matrix in (matrices or {}).items():
        for index in range(matrix.shape[1]):
            columns.append(matrix[:, index])
            names.append(f"{prefix}_{index + 1}")
    return np.column_stack(columns), names


def survey_wls(
    data: pd.DataFrame,
    outcome: str,
    continuous: list[str],
    categorical: list[str],
    weight: str,
    matrices: dict[str, np.ndarray] | None = None,
) -> dict:
    x, names = make_design(data, continuous, categorical, matrices)
    y = data[outcome].to_numpy(dtype=float)
    w = data[weight].to_numpy(dtype=float)
    bread = np.linalg.inv(x.T @ (w[:, None] * x))
    beta = bread @ (x.T @ (w * y))
    residual = y - x @ beta
    scores = x * (w * residual)[:, None]
    meat = np.zeros((x.shape[1], x.shape[1]))
    psu_total = 0
    strata_total = 0
    reset = data.reset_index(drop=True)
    for _, indices in reset.groupby("SDMVSTRA").groups.items():
        indices = np.asarray(list(indices))
        psu_values = reset.loc[indices, "SDMVPSU"].to_numpy()
        psu_scores = np.vstack(
            [
                scores[indices[psu_values == psu]].sum(axis=0)
                for psu in np.unique(psu_values)
            ]
        )
        count = len(psu_scores)
        if count > 1:
            psu_scores -= psu_scores.mean(axis=0)
            meat += count / (count - 1) * (psu_scores.T @ psu_scores)
            psu_total += count
            strata_total += 1
    covariance = bread @ meat @ bread
    standard_error = np.sqrt(np.diag(covariance))
    design_df = psu_total - strata_total
    residual_df = design_df - len(beta) + 1
    critical = t_dist.ppf(0.975, residual_df)
    p_values = 2 * t_dist.sf(np.abs(beta / standard_error), residual_df)
    table = pd.DataFrame(
        {
            "term": names,
            "estimate": beta,
            "standard_error": standard_error,
            "ci_low": beta - critical * standard_error,
            "ci_high": beta + critical * standard_error,
            "p_value": p_values,
        }
    )
    return {
        "table": table,
        "beta": beta,
        "covariance": covariance,
        "names": names,
        "design_df": design_df,
        "residual_df": residual_df,
        "n": len(data),
    }


def complete_cases(data: pd.DataFrame, variables: list[str]) -> pd.DataFrame:
    return data.dropna(subset=variables).copy()


def summarize_term(model: dict, term: str) -> dict:
    row = model["table"].loc[model["table"]["term"] == term].iloc[0]
    return {
        "n": model["n"],
        "estimate": row["estimate"],
        "standard_error": row["standard_error"],
        "ci_low": row["ci_low"],
        "ci_high": row["ci_high"],
        "p_value": row["p_value"],
        "residual_df": model["residual_df"],
    }


def wald_test(model: dict, terms: list[str]) -> dict:
    """Joint survey-design Wald F test for a set of model coefficients."""
    indices = [model["names"].index(term) for term in terms]
    beta = model["beta"][indices]
    covariance = model["covariance"][np.ix_(indices, indices)]
    df1 = len(indices)
    statistic = float(beta.T @ np.linalg.pinv(covariance) @ beta / df1)
    return {
        "f_statistic": statistic,
        "df1": df1,
        "df2": model["residual_df"],
        "p_value": float(f_dist.sf(statistic, df1, model["residual_df"])),
    }


def run() -> None:
    built = [build_cycle(cycle) for cycle in CYCLES]
    data = pd.concat([item[0] for item in built], ignore_index=True)
    products = pd.concat([item[1] for item in built], ignore_index=True)

    complete_battery = data[COGNITIVE_COMPONENTS].notna().all(axis=1)
    for variable in COGNITIVE_COMPONENTS:
        mean = data.loc[complete_battery, variable].mean()
        sd = data.loc[complete_battery, variable].std()
        data[f"{variable}_z"] = (data[variable] - mean) / sd
    data["global_cognition_complete"] = data[
        [f"{variable}_z" for variable in COGNITIVE_COMPONENTS]
    ].mean(axis=1)

    common = CONTINUOUS_COVARIATES + CATEGORICAL_COVARIATES
    accepted = complete_cases(
        data,
        ["global_cognition_accepted", "omega_user_24h", "WTMEC2YR"] + common,
    )
    complete = complete_cases(
        data.loc[complete_battery],
        ["global_cognition_complete", "omega_user_30d", "WTMEC2YR"] + common,
    )
    complete["formulation_group"] = complete.apply(
        lambda row: collapsed_formulation(
            row["omega_user_30d"], row["omega_formulations_30d"]
        ),
        axis=1,
    )
    formulation_users = complete.loc[complete["omega_user_30d"].eq(1)].copy()
    if formulation_users["formulation_group"].eq("Unclassified").any():
        raise ValueError("At least one omega-3 user could not be classified")
    formulation_indicators = {
        "formulation_krill": "Krill",
        "formulation_plant_based": "Plant-based",
        "formulation_other_mixed": "Other/mixed",
    }
    for variable, label in formulation_indicators.items():
        formulation_users[variable] = formulation_users[
            "formulation_group"
        ].eq(label).astype(int)
    dietary = complete_cases(
        data.loc[complete_battery & data["DR1DRSTZ"].eq(1)],
        [
            "global_cognition_complete",
            "log_dietary_epa_dha",
            "dietary_epa_dha_g",
            "DR1TKCAL",
            "WTDRD1",
            "omega_user_30d",
        ]
        + common,
    )

    models: dict[str, dict] = {}
    models["accepted_replication"] = survey_wls(
        accepted,
        "global_cognition_accepted",
        ["omega_user_24h"] + CONTINUOUS_COVARIATES,
        CATEGORICAL_COVARIATES,
        "WTMEC2YR",
    )
    models["supplement_30d_complete"] = survey_wls(
        complete,
        "global_cognition_complete",
        ["omega_user_30d"] + CONTINUOUS_COVARIATES,
        CATEGORICAL_COVARIATES,
        "WTMEC2YR",
    )
    models["formulation_exploratory"] = survey_wls(
        formulation_users,
        "global_cognition_complete",
        list(formulation_indicators) + CONTINUOUS_COVARIATES,
        CATEGORICAL_COVARIATES,
        "WTMEC2YR",
    )
    models["diet_correct_weight"] = survey_wls(
        dietary,
        "global_cognition_complete",
        ["log_dietary_epa_dha"] + CONTINUOUS_COVARIATES,
        CATEGORICAL_COVARIATES,
        "WTDRD1",
    )
    models["joint_energy"] = survey_wls(
        dietary,
        "global_cognition_complete",
        [
            "log_dietary_epa_dha",
            "omega_user_30d",
            "DR1TKCAL",
        ]
        + CONTINUOUS_COVARIATES,
        CATEGORICAL_COVARIATES,
        "WTDRD1",
    )

    lifestyle_vars = [
        "smoking_status",
        "current_alcohol",
        "recreational_activity",
    ]
    lifestyle = complete_cases(complete, lifestyle_vars)
    models["supplement_lifestyle"] = survey_wls(
        lifestyle,
        "global_cognition_complete",
        [
            "omega_user_30d",
            "current_alcohol",
            "recreational_activity",
        ]
        + CONTINUOUS_COVARIATES,
        CATEGORICAL_COVARIATES + ["smoking_status"],
        "WTMEC2YR",
    )

    quantiles = np.array([0.05, 0.35, 0.65, 0.95])
    knots = weighted_quantile(
        dietary["log_dietary_epa_dha"].to_numpy(),
        quantiles,
        dietary["WTDRD1"].to_numpy(),
    )
    spline = rcs_basis(dietary["log_dietary_epa_dha"].to_numpy(), knots)
    models["diet_spline"] = survey_wls(
        dietary,
        "global_cognition_complete",
        ["omega_user_30d", "DR1TKCAL"] + CONTINUOUS_COVARIATES,
        CATEGORICAL_COVARIATES,
        "WTDRD1",
        matrices={"diet_rcs": spline},
    )
    spline_model = models["diet_spline"]
    nonlinear_names = [name for name in spline_model["names"] if name.startswith("diet_rcs_")][1:]
    nonlinear_indices = [spline_model["names"].index(name) for name in nonlinear_names]
    nonlinear_beta = spline_model["beta"][nonlinear_indices]
    nonlinear_cov = spline_model["covariance"][np.ix_(nonlinear_indices, nonlinear_indices)]
    nonlinear_f = float(
        nonlinear_beta.T @ np.linalg.inv(nonlinear_cov) @ nonlinear_beta
        / len(nonlinear_indices)
    )
    nonlinear_p = float(
        f_dist.sf(
            nonlinear_f,
            len(nonlinear_indices),
            spline_model["residual_df"],
        )
    )
    formulation_test = wald_test(
        models["formulation_exploratory"], list(formulation_indicators)
    )
    formulation_counts = (
        formulation_users["formulation_group"]
        .value_counts()
        .rename_axis("formulation_group")
        .reset_index(name="unweighted_n")
    )

    key_results = {
        "accepted_24h_partial_battery": summarize_term(
            models["accepted_replication"], "omega_user_24h"
        ),
        "corrected_30d_complete_battery": summarize_term(
            models["supplement_30d_complete"], "omega_user_30d"
        ),
        "diet_correct_weight": summarize_term(
            models["diet_correct_weight"], "log_dietary_epa_dha"
        ),
        "joint_energy_diet": summarize_term(
            models["joint_energy"], "log_dietary_epa_dha"
        ),
        "joint_energy_supplement": summarize_term(
            models["joint_energy"], "omega_user_30d"
        ),
        "lifestyle_supplement": summarize_term(
            models["supplement_lifestyle"], "omega_user_30d"
        ),
        "exploratory_formulation": {
            "population": "Past-30-day omega-3 supplement users only",
            "reference": "Fish oil/cod liver oil",
            "subgroup_counts": dict(
                zip(
                    formulation_counts["formulation_group"],
                    formulation_counts["unweighted_n"].astype(int),
                )
            ),
            "overall_wald_test": formulation_test,
            "contrasts_vs_reference": {
                label: summarize_term(
                    models["formulation_exploratory"], variable
                )
                for variable, label in formulation_indicators.items()
            },
        },
        "spline_non_linearity": {
            "knots_log_grams": knots.tolist(),
            "f_statistic": nonlinear_f,
            "df1": len(nonlinear_indices),
            "df2": spline_model["residual_df"],
            "p_value": nonlinear_p,
        },
        "sample_flow": {
            "eligible_age_60_plus": len(data),
            "accepted_complete_case": len(accepted),
            "complete_battery_complete_case": len(complete),
            "diet_reliable_complete_case": len(dietary),
            "lifestyle_complete_case": len(lifestyle),
            "accepted_partial_battery_count": int(
                (~accepted[COGNITIVE_COMPONENTS].notna().all(axis=1)).sum()
            ),
            "accepted_24h_users": int(accepted["omega_user_24h"].sum()),
            "accepted_30d_users": int(accepted["omega_user_30d"].sum()),
        },
    }

    with (OUT / "key_results.json").open("w") as stream:
        json.dump(key_results, stream, indent=2)
    data.to_csv(OUT / "rebuilt_analytic_data.csv", index=False)
    products.to_csv(OUT / "omega3_products_30d.csv", index=False)
    formulation_counts.to_csv(OUT / "formulation_counts.csv", index=False)
    for name, model in models.items():
        model["table"].to_csv(OUT / f"{name}.csv", index=False)

    # Adjusted spline curve: effect relative to 100 mg/day at fixed covariates.
    grid_mg = np.linspace(0, 1000, 201)
    grid_log_g = np.log1p(grid_mg / 1000)
    grid_basis = rcs_basis(grid_log_g, knots)
    reference_basis = rcs_basis(np.array([np.log1p(0.1)]), knots)[0]
    spline_indices = [
        spline_model["names"].index(name)
        for name in spline_model["names"]
        if name.startswith("diet_rcs_")
    ]
    contrast = grid_basis - reference_basis
    spline_beta = spline_model["beta"][spline_indices]
    spline_cov = spline_model["covariance"][np.ix_(spline_indices, spline_indices)]
    estimate = contrast @ spline_beta
    variance = np.einsum("ij,jk,ik->i", contrast, spline_cov, contrast)
    critical = t_dist.ppf(0.975, spline_model["residual_df"])
    lower = estimate - critical * np.sqrt(variance)
    upper = estimate + critical * np.sqrt(variance)

    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=200)
    ax.fill_between(grid_mg, lower, upper, color="#bfd7ea", alpha=0.8, linewidth=0)
    ax.plot(grid_mg, estimate, color="#1261a0", linewidth=2.5)
    ax.axhline(0, color="#444444", linewidth=1, linestyle="--")
    ax.axvline(100, color="#888888", linewidth=1, linestyle=":")
    ax.set(
        xlabel="Dietary EPA + DHA (mg/day)",
        ylabel="Adjusted difference in cognitive z-score\n(reference: 100 mg/day)",
        title=f"Exploratory dose-response (P non-linearity = {nonlinear_p:.3f})",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "dietary_spline.png", transparent=False)
    plt.close(fig)

    print(json.dumps(key_results, indent=2))


if __name__ == "__main__":
    run()
