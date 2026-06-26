#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg

from run_targeted_phewas import bh_fdr, normalize_bool, prepare_design
from run_survival_phewas import age_years, build_event_dates, normalize_sex, prepare_covariates


EXPOSURE_Z = "somatic_burden_z"
DAYS_PER_YEAR = 365.25


def standardize_exposure(exposure, exposure_col):
    if exposure_col not in exposure.columns:
        raise ValueError(f"Exposure column not found: {exposure_col}")
    out = exposure[["IID", exposure_col]].copy()
    out[exposure_col] = pd.to_numeric(out[exposure_col], errors="coerce")
    out = out.dropna(subset=[exposure_col])
    sd = out[exposure_col].std()
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Exposure has zero or non-finite standard deviation.")
    out[EXPOSURE_Z] = (out[exposure_col] - out[exposure_col].mean()) / sd
    return out[["IID", exposure_col, EXPOSURE_Z]]


def run_logistic(exposure, phenotypes, summary, covariates, age_cols, pc_cols, sex_col):
    base_cols = age_cols + pc_cols
    extra_cols = ["ever_smoked"] if "ever_smoked" in covariates.columns else []
    df = phenotypes.merge(exposure[["IID", EXPOSURE_Z]], on="IID", how="inner")
    df = df.merge(covariates[["IID", sex_col] + base_cols + extra_cols], on="IID", how="inner")

    rows = []
    for _, meta in summary.iterrows():
        phenotype = meta["phenotype_id"]
        if not normalize_bool(meta.get("include", False)) or phenotype not in df.columns:
            rows.append({
                "analysis": "logistic_any_time",
                "phenotype_id": phenotype,
                "label": meta.get("label", ""),
                "phecodes": meta.get("phecodes", ""),
                "sex_restriction": meta.get("sex_restriction", ""),
                "category": meta.get("category", ""),
                "status": meta.get("status", "excluded"),
                "n": 0,
                "cases": meta.get("cases", np.nan),
                "controls": meta.get("controls", np.nan),
                "beta": np.nan,
                "se": np.nan,
                "or": np.nan,
                "ci_lower": np.nan,
                "ci_upper": np.nan,
                "p": np.nan,
            })
            continue

        covar_cols = base_cols + extra_cols
        if str(meta.get("sex_restriction", "both")).lower() == "both":
            covar_cols.append(sex_col)

        model_df = df[[phenotype, EXPOSURE_Z] + covar_cols].dropna()
        y = pd.to_numeric(model_df[phenotype], errors="coerce")
        x = prepare_design(model_df, EXPOSURE_Z, covar_cols)
        ok = y.notna() & x.notna().all(axis=1)
        y = y.loc[ok].astype(float)
        x = x.loc[ok].astype(float)
        cases = int(y.sum())
        controls = int(len(y) - cases)

        row = {
            "analysis": "logistic_any_time",
            "phenotype_id": phenotype,
            "label": meta.get("label", ""),
            "phecodes": meta.get("phecodes", ""),
            "sex_restriction": meta.get("sex_restriction", ""),
            "category": meta.get("category", ""),
            "status": "ok",
            "n": len(y),
            "cases": cases,
            "controls": controls,
            "beta": np.nan,
            "se": np.nan,
            "or": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p": np.nan,
        }
        try:
            fit = sm.Logit(y, x).fit(disp=False, maxiter=100)
            beta = fit.params[EXPOSURE_Z]
            se = fit.bse[EXPOSURE_Z]
            row.update({
                "beta": beta,
                "se": se,
                "or": math.exp(beta),
                "ci_lower": math.exp(beta - 1.96 * se),
                "ci_upper": math.exp(beta + 1.96 * se),
                "p": fit.pvalues[EXPOSURE_Z],
            })
        except Exception as exc:
            row["status"] = f"failed: {exc}"
        rows.append(row)

    out = pd.DataFrame(rows)
    tested = out["p"].notna().sum()
    out["bonferroni_p"] = np.where(out["p"].notna(), np.minimum(out["p"] * tested, 1), np.nan)
    out["fdr_bh"] = bh_fdr(out["p"].to_numpy()) if len(out) else []
    return out.sort_values(["p", "phenotype_id"], na_position="last")


def fit_phreg(model_df, exposure_col, covariate_cols):
    x = prepare_covariates(model_df, [exposure_col] + covariate_cols)
    ok = (
        model_df["duration_age"].notna()
        & model_df["entry_age"].notna()
        & model_df["status"].notna()
        & x.notna().all(axis=1)
    )
    ok &= model_df["duration_age"].gt(model_df["entry_age"])
    y = model_df.loc[ok, "duration_age"].astype(float)
    entry = model_df.loc[ok, "entry_age"].astype(float)
    status = model_df.loc[ok, "status"].astype(int)
    x = x.loc[ok].astype(float)
    fit = PHReg(y, x, status=status, entry=entry, ties="breslow").fit()
    beta = fit.params[list(x.columns).index(exposure_col)]
    se = fit.bse[list(x.columns).index(exposure_col)]
    p = fit.pvalues[list(x.columns).index(exposure_col)]
    result = {
        "n": len(y),
        "events": int(status.sum()),
        "censored": int(len(status) - status.sum()),
        "beta": beta,
        "se": se,
        "hr": math.exp(beta),
        "ci_lower": math.exp(beta - 1.96 * se),
        "ci_upper": math.exp(beta + 1.96 * se),
        "p": p,
    }
    values = [result["beta"], result["se"], result["hr"], result["ci_lower"], result["ci_upper"], result["p"]]
    if not all(np.isfinite(value) for value in values):
        raise ValueError("PHReg returned non-finite exposure estimate.")
    return result


def fit_phreg_with_fallbacks(model_df, exposure_col, full_covariates, pc_cols, extra_cols, include_sex):
    fallback_sets = [("full_pc16", full_covariates)]
    pc4 = [col for col in pc_cols[:4] if col in full_covariates]
    reduced = pc4 + [col for col in extra_cols if col in full_covariates]
    if include_sex and "sex_normalized" in full_covariates:
        reduced.append("sex_normalized")
    fallback_sets.append(("reduced_pc4", reduced))
    minimal = [col for col in extra_cols if col in full_covariates]
    if include_sex and "sex_normalized" in full_covariates:
        minimal.append("sex_normalized")
    fallback_sets.append(("minimal", minimal))

    errors = []
    for model_name, covariates in fallback_sets:
        try:
            result = fit_phreg(model_df, exposure_col, covariates)
            result["covariate_model"] = model_name
            return result
        except Exception as exc:
            errors.append(f"{model_name}: {exc}")
    raise RuntimeError("; ".join(errors))


def run_post_draw_survival(
    exposure,
    target_map,
    evidence,
    summary,
    survival,
    age_at_sequencing,
    covariates,
    min_code_count,
    pc_cols,
    sex_col,
):
    covariates = covariates.copy()
    covariates["sex_normalized"] = covariates[sex_col].map(normalize_sex)
    extra_cols = ["ever_smoked"] if "ever_smoked" in covariates.columns else []

    age = age_at_sequencing.rename(columns={"ID": "IID"}).copy()
    age["age_at_sequencing"] = pd.to_numeric(age["age_at_sequencing"], errors="coerce")
    base = (
        exposure[["IID", EXPOSURE_Z]]
        .merge(survival, on="IID", how="inner")
        .merge(age[["IID", "age_at_sequencing"]], on="IID", how="inner")
        .merge(covariates[["IID", "sex_normalized"] + pc_cols + extra_cols], on="IID", how="inner")
    )
    for col in ["birth_date", "observation_start_date", "observation_end_date", "death_date"]:
        base[col] = pd.to_datetime(base[col], errors="coerce")
    base = base.dropna(subset=["birth_date", "observation_start_date", "observation_end_date", "age_at_sequencing"])
    base["observation_start_age"] = age_years(base["observation_start_date"], base["birth_date"])
    base["obs_end_age"] = age_years(base["observation_end_date"], base["birth_date"])
    base["death_age"] = age_years(base["death_date"], base["birth_date"])
    base["entry_age"] = base[["observation_start_age", "age_at_sequencing"]].max(axis=1)
    base = base[base["obs_end_age"].gt(base["entry_age"])].copy()

    target_map = target_map[target_map["map_status"].eq("ok")]
    event_dates = build_event_dates(evidence, target_map, min_code_count, set(base["IID"]))

    rows = []
    for _, meta in summary.iterrows():
        phenotype = meta["phenotype_id"]
        covar_cols = pc_cols + extra_cols
        sex_restriction = str(meta.get("sex_restriction", "both")).lower()
        model_base = base.copy()
        if sex_restriction in {"female", "male"}:
            model_base = model_base[model_base["sex_normalized"].eq(sex_restriction)].copy()
        else:
            covar_cols.append("sex_normalized")

        phenotype_events = event_dates[event_dates["phenotype_id"].eq(phenotype)][["IID", "event_date"]]
        model_df = model_base.merge(phenotype_events, on="IID", how="left")
        model_df["event_age"] = age_years(model_df["event_date"], model_df["birth_date"])
        model_df["prevalent_before_draw"] = model_df["event_age"].notna() & model_df["event_age"].le(model_df["entry_age"])
        prevalent_excluded = int(model_df["prevalent_before_draw"].sum())
        model_df = model_df[~model_df["prevalent_before_draw"]].copy()
        censor_age = model_df[["obs_end_age", "death_age"]].min(axis=1, skipna=True)
        model_df["status"] = model_df["event_age"].notna() & model_df["event_age"].le(censor_age)
        model_df["duration_age"] = np.where(model_df["status"], model_df["event_age"], censor_age)

        include_sex = "sex_normalized" in covar_cols
        row = {
            "analysis": "post_blood_draw_survival",
            "phenotype_id": phenotype,
            "label": meta.get("label", ""),
            "outcome_type": "phecode",
            "sex_restriction": sex_restriction,
            "covariate_model": "not_fit",
            "prevalent_before_draw_excluded": prevalent_excluded,
            "n": np.nan,
            "events": np.nan,
            "censored": np.nan,
            "beta": np.nan,
            "se": np.nan,
            "hr": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p": np.nan,
            "status": "ok",
        }
        try:
            row.update(fit_phreg_with_fallbacks(model_df, EXPOSURE_Z, covar_cols, pc_cols, extra_cols, include_sex))
        except Exception as exc:
            row["status"] = f"failed: {exc}"
        rows.append(row)

    mortality = base.copy()
    mortality["status"] = mortality["death_age"].notna() & mortality["death_age"].le(mortality["obs_end_age"])
    mortality["duration_age"] = np.where(mortality["status"], mortality["death_age"], mortality["obs_end_age"])
    row = {
        "analysis": "post_blood_draw_survival",
        "phenotype_id": "mortality",
        "label": "All-cause mortality",
        "outcome_type": "mortality",
        "sex_restriction": "both",
        "covariate_model": "not_fit",
        "prevalent_before_draw_excluded": 0,
        "n": np.nan,
        "events": np.nan,
        "censored": np.nan,
        "beta": np.nan,
        "se": np.nan,
        "hr": np.nan,
        "ci_lower": np.nan,
        "ci_upper": np.nan,
        "p": np.nan,
        "status": "ok",
    }
    try:
        row.update(
            fit_phreg_with_fallbacks(
                mortality,
                EXPOSURE_Z,
                pc_cols + extra_cols + ["sex_normalized"],
                pc_cols,
                extra_cols,
                True,
            )
        )
    except Exception as exc:
        row["status"] = f"failed: {exc}"
    rows.append(row)

    out = pd.DataFrame(rows)
    tested = out["p"].notna().sum()
    out["bonferroni_p"] = np.where(out["p"].notna(), np.minimum(out["p"] * tested, 1), np.nan)
    out["fdr_bh"] = bh_fdr(out["p"].to_numpy()) if len(out) else []
    return out.sort_values(["p", "phenotype_id"], na_position="last")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exposure", required=True)
    parser.add_argument("--exposure-col", default="V3")
    parser.add_argument("--phenotypes", required=True)
    parser.add_argument("--phenotype-summary", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--target-map", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--survival-base", required=True)
    parser.add_argument("--age-at-sequencing", required=True)
    parser.add_argument("--min-code-count", type=int, default=1)
    parser.add_argument("--age-cols", required=True)
    parser.add_argument("--pc-cols", required=True)
    parser.add_argument("--sex-col", default="sex")
    parser.add_argument("--logistic-out", required=True)
    parser.add_argument("--survival-out", required=True)
    args = parser.parse_args()

    exposure = pd.read_csv(args.exposure, sep="\t", dtype={"IID": str})
    exposure = standardize_exposure(exposure, args.exposure_col)
    phenotypes = pd.read_csv(args.phenotypes, sep="\t", dtype={"IID": str})
    summary = pd.read_csv(args.phenotype_summary, sep="\t", dtype={"phenotype_id": str})
    covariates = pd.read_csv(args.covariates, sep="\t", dtype={"IID": str})
    target_map = pd.read_csv(args.target_map, sep="\t", dtype=str)
    evidence = pd.read_csv(args.evidence, sep="\t", dtype=str)
    survival = pd.read_csv(args.survival_base, sep="\t", dtype={"IID": str})
    age_at_sequencing = pd.read_csv(args.age_at_sequencing, sep="\t", dtype={"ID": str})

    age_cols = [x for x in args.age_cols.split(",") if x]
    pc_cols = [x for x in args.pc_cols.split(",") if x]

    logistic = run_logistic(exposure, phenotypes, summary, covariates, age_cols, pc_cols, args.sex_col)
    Path(args.logistic_out).parent.mkdir(parents=True, exist_ok=True)
    logistic.to_csv(args.logistic_out, sep="\t", index=False)

    survival_out = run_post_draw_survival(
        exposure,
        target_map,
        evidence,
        summary,
        survival,
        age_at_sequencing,
        covariates,
        args.min_code_count,
        pc_cols,
        args.sex_col,
    )
    Path(args.survival_out).parent.mkdir(parents=True, exist_ok=True)
    survival_out.to_csv(args.survival_out, sep="\t", index=False)

    print(f"logistic phenotypes tested: {logistic['p'].notna().sum()}")
    print(f"post-blood-draw survival outcomes tested: {survival_out['p'].notna().sum()}")
    print(f"exposure N: {len(exposure)}")


if __name__ == "__main__":
    main()
