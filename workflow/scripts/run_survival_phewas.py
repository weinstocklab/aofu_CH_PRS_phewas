#!/usr/bin/env python3
import argparse
import math

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg

from run_targeted_phewas import bh_fdr, normalize_bool


DAYS_PER_YEAR = 365.25


def normalize_sex(value):
    value = str(value).strip().lower()
    if value.startswith("f"):
        return "female"
    if value.startswith("m"):
        return "male"
    return "unknown"


def age_years(date, birth_date):
    return (date - birth_date).dt.days / DAYS_PER_YEAR


def prepare_covariates(df, covariate_cols):
    x = df[covariate_cols].copy()
    for col in list(x.columns):
        if x[col].dtype == object or str(x[col].dtype) == "bool":
            x = pd.get_dummies(x, columns=[col], drop_first=True, dtype=float)
    return x.apply(pd.to_numeric, errors="coerce")


def build_event_dates(evidence, target_map, min_code_count, cohort_iids):
    mapped = evidence.merge(
        target_map[["phenotype_id", "vocabulary_id", "concept_code"]],
        on=["vocabulary_id", "concept_code"],
        how="inner",
    )
    mapped = mapped[mapped["IID"].isin(cohort_iids)].copy()
    mapped["code_date"] = pd.to_datetime(mapped["code_date"], errors="coerce")
    mapped = mapped.dropna(subset=["code_date"])
    mapped = mapped.drop_duplicates(["IID", "phenotype_id", "code_date"])
    mapped = mapped.sort_values(["IID", "phenotype_id", "code_date"])
    mapped["code_order"] = mapped.groupby(["IID", "phenotype_id"]).cumcount() + 1
    event_dates = mapped[mapped["code_order"].eq(min_code_count)][["IID", "phenotype_id", "code_date"]]
    return event_dates.rename(columns={"code_date": "event_date"})


def fit_phreg(model_df, prs_col, covariate_cols):
    x = prepare_covariates(model_df, [prs_col] + covariate_cols)
    ok = model_df["duration_age"].notna() & model_df["entry_age"].notna() & model_df["status"].notna() & x.notna().all(axis=1)
    ok &= model_df["duration_age"].gt(model_df["entry_age"])
    y = model_df.loc[ok, "duration_age"].astype(float)
    entry = model_df.loc[ok, "entry_age"].astype(float)
    status = model_df.loc[ok, "status"].astype(int)
    x = x.loc[ok].astype(float)
    fit = PHReg(y, x, status=status, entry=entry, ties="breslow").fit()
    beta = fit.params[list(x.columns).index(prs_col)]
    se = fit.bse[list(x.columns).index(prs_col)]
    p = fit.pvalues[list(x.columns).index(prs_col)]
    return {
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True)
    parser.add_argument("--target-map", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--phenotype-summary", required=True)
    parser.add_argument("--survival-base", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--prs-col", default="PRSFNN_out_final_Z")
    parser.add_argument("--pc-cols", required=True)
    parser.add_argument("--sex-col", default="sex")
    parser.add_argument("--min-code-count", type=int, default=1)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    scores = pd.read_csv(args.scores, sep="\t", dtype={"IID": str})
    target_map = pd.read_csv(args.target_map, sep="\t", dtype=str)
    target_map = target_map[target_map["map_status"].eq("ok")]
    evidence = pd.read_csv(args.evidence, sep="\t", dtype=str)
    summary = pd.read_csv(args.phenotype_summary, sep="\t", dtype={"phenotype_id": str})
    survival = pd.read_csv(args.survival_base, sep="\t", dtype={"IID": str})
    covariates = pd.read_csv(args.covariates, sep="\t", dtype={"IID": str})

    pc_cols = [x for x in args.pc_cols.split(",") if x]
    covariates["sex_normalized"] = covariates[args.sex_col].map(normalize_sex)
    base = (
        scores[["IID", args.prs_col]]
        .merge(survival, on="IID", how="inner")
        .merge(covariates[["IID", "sex_normalized"] + pc_cols], on="IID", how="inner")
    )
    for col in ["birth_date", "observation_start_date", "observation_end_date", "death_date"]:
        base[col] = pd.to_datetime(base[col], errors="coerce")
    base = base.dropna(subset=["birth_date", "observation_start_date", "observation_end_date"])
    base["entry_age"] = age_years(base["observation_start_date"], base["birth_date"])
    base["obs_end_age"] = age_years(base["observation_end_date"], base["birth_date"])
    base["death_age"] = age_years(base["death_date"], base["birth_date"])
    base = base[base["obs_end_age"].gt(base["entry_age"])].copy()

    event_dates = build_event_dates(evidence, target_map, args.min_code_count, set(base["IID"]))

    rows = []
    for _, meta in summary.iterrows():
        phenotype = meta["phenotype_id"]
        covar_cols = pc_cols.copy()
        sex_restriction = str(meta.get("sex_restriction", "both")).lower()
        model_base = base.copy()
        if sex_restriction in {"female", "male"}:
            model_base = model_base[model_base["sex_normalized"].eq(sex_restriction)].copy()
        else:
            covar_cols.append("sex_normalized")

        phenotype_events = event_dates[event_dates["phenotype_id"].eq(phenotype)][["IID", "event_date"]]
        model_df = model_base.merge(phenotype_events, on="IID", how="left")
        model_df["event_age"] = age_years(model_df["event_date"], model_df["birth_date"])
        model_df["prevalent_at_entry"] = model_df["event_age"].notna() & model_df["event_age"].le(model_df["entry_age"])
        model_df = model_df[~model_df["prevalent_at_entry"]].copy()
        censor_age = model_df[["obs_end_age", "death_age"]].min(axis=1, skipna=True)
        model_df["status"] = model_df["event_age"].notna() & model_df["event_age"].le(censor_age)
        model_df["duration_age"] = np.where(model_df["status"], model_df["event_age"], censor_age)

        row = {
            "phenotype_id": phenotype,
            "label": meta.get("label", ""),
            "outcome_type": "phecode",
            "sex_restriction": sex_restriction,
            "prevalent_excluded": int((model_base.merge(phenotype_events, on="IID", how="left").assign(event_age=lambda d: age_years(d["event_date"], d["birth_date"]))["event_age"] <= model_base.merge(phenotype_events, on="IID", how="left").assign(event_age=lambda d: age_years(d["event_date"], d["birth_date"]))["entry_age"]).sum()),
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
            row.update(fit_phreg(model_df, args.prs_col, covar_cols))
        except Exception as exc:
            row["status"] = f"failed: {exc}"
        rows.append(row)

    mortality = base.copy()
    mortality["status"] = mortality["death_age"].notna() & mortality["death_age"].le(mortality["obs_end_age"])
    mortality["duration_age"] = np.where(mortality["status"], mortality["death_age"], mortality["obs_end_age"])
    row = {
        "phenotype_id": "mortality",
        "label": "All-cause mortality",
        "outcome_type": "mortality",
        "sex_restriction": "both",
        "prevalent_excluded": 0,
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
        row.update(fit_phreg(mortality, args.prs_col, pc_cols + ["sex_normalized"]))
    except Exception as exc:
        row["status"] = f"failed: {exc}"
    rows.append(row)

    out = pd.DataFrame(rows)
    tested = out["p"].notna().sum()
    out["bonferroni_p"] = np.where(out["p"].notna(), np.minimum(out["p"] * tested, 1), np.nan)
    out["fdr_bh"] = bh_fdr(out["p"].to_numpy()) if len(out) else []
    out = out.sort_values(["p", "phenotype_id"], na_position="last")
    out.to_csv(args.out, sep="\t", index=False)
    print(f"tested {tested} survival outcomes")


if __name__ == "__main__":
    main()
