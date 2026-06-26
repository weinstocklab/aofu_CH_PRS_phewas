#!/usr/bin/env python3
import argparse
import math

import numpy as np
import pandas as pd
import statsmodels.api as sm


def is_binary(series):
    values = set(series.dropna().unique())
    return values.issubset({0, 1, 0.0, 1.0, False, True})


def prepare_design(df, score_col, covariate_cols):
    x = df[[score_col] + covariate_cols].copy()
    for col in x.columns:
        if x[col].dtype == object or str(x[col].dtype) == "bool":
            if set(x[col].dropna().astype(str).str.lower().unique()).issubset({"true", "false"}):
                x[col] = x[col].astype(str).str.lower().map({"true": 1.0, "false": 0.0})
            else:
                x = pd.get_dummies(x, columns=[col], drop_first=True, dtype=float)
    x = x.apply(pd.to_numeric, errors="coerce")
    return sm.add_constant(x, has_constant="add")


def fit_one(df, phenotype, score_col, covariate_cols, min_cases):
    y_raw = df[phenotype]
    model_df = df[[phenotype, score_col] + covariate_cols].dropna()
    n = len(model_df)
    if n == 0:
        return None

    y = model_df[phenotype]
    binary = is_binary(y)
    cases = int(y.sum()) if binary else math.nan
    controls = int(n - cases) if binary else math.nan
    if binary and (cases < min_cases or controls < min_cases):
        return None

    x = prepare_design(model_df, score_col, covariate_cols)
    y = pd.to_numeric(y, errors="coerce")
    ok = y.notna() & x.notna().all(axis=1)
    y = y.loc[ok]
    x = x.loc[ok]
    if len(y) == 0:
        return None

    try:
        if binary:
            result = sm.Logit(y.astype(float), x.astype(float)).fit(disp=False, maxiter=100)
            beta = result.params[score_col]
            se = result.bse[score_col]
            p = result.pvalues[score_col]
            effect = math.exp(beta)
            effect_type = "OR"
        else:
            result = sm.OLS(y.astype(float), x.astype(float)).fit()
            beta = result.params[score_col]
            se = result.bse[score_col]
            p = result.pvalues[score_col]
            effect = beta
            effect_type = "BETA"
    except Exception as exc:
        return {
            "phenotype": phenotype,
            "model": "logistic" if binary else "linear",
            "n": len(y),
            "cases": cases,
            "controls": controls,
            "beta": np.nan,
            "se": np.nan,
            "p": np.nan,
            "effect": np.nan,
            "effect_type": "OR" if binary else "BETA",
            "status": f"failed: {exc}",
        }

    return {
        "phenotype": phenotype,
        "model": "logistic" if binary else "linear",
        "n": len(y),
        "cases": cases,
        "controls": controls,
        "beta": beta,
        "se": se,
        "p": p,
        "effect": effect,
        "effect_type": effect_type,
        "status": "ok",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True)
    parser.add_argument("--phenotypes", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--covariate-cols", required=True)
    parser.add_argument("--min-cases", type=int, default=50)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    scores = pd.read_csv(args.scores, sep="\t", dtype={"IID": str})
    score_cols = [col for col in scores.columns if col.endswith("_Z")]
    if len(score_cols) != 1:
        raise ValueError(f"Expected exactly one *_Z column in score file, found {score_cols}")
    score_col = score_cols[0]

    covariates = pd.read_csv(args.covariates, sep="\t", dtype={"IID": str})
    phenotypes = pd.read_csv(args.phenotypes, sep="\t", dtype={"IID": str})
    covariate_cols = [col for col in args.covariate_cols.split(",") if col]

    missing = [col for col in covariate_cols if col not in covariates.columns]
    if missing:
        raise ValueError(f"Missing covariate columns: {missing}")

    df = phenotypes.merge(scores[["IID", score_col]], on="IID", how="inner")
    df = df.merge(covariates[["IID"] + covariate_cols], on="IID", how="inner")

    phenotype_cols = [col for col in phenotypes.columns if col not in {"FID", "IID"}]
    results = []
    for phenotype in phenotype_cols:
        result = fit_one(df, phenotype, score_col, covariate_cols, args.min_cases)
        if result is not None:
            results.append(result)

    out = pd.DataFrame(results)
    if not out.empty and "p" in out:
        out["bonferroni_p"] = out["p"] * out["p"].notna().sum()
        out["bonferroni_p"] = out["bonferroni_p"].clip(upper=1)
        out = out.sort_values("p", na_position="last")
    out.to_csv(args.out, sep="\t", index=False)
    print(f"tested {len(out)} phenotypes")


if __name__ == "__main__":
    main()
