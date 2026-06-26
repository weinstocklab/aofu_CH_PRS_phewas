#!/usr/bin/env python3
import argparse
import math

import numpy as np
import pandas as pd
import patsy
import scipy.stats as st
import statsmodels.api as sm

from run_targeted_phewas import bh_fdr, normalize_bool, prepare_design


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True)
    parser.add_argument("--phenotypes", required=True)
    parser.add_argument("--phenotype-summary", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--prs-col", default="PRSFNN_out_final_Z")
    parser.add_argument("--age-cols", required=True)
    parser.add_argument("--pc-cols", required=True)
    parser.add_argument("--sex-col", default="sex")
    parser.add_argument("--spline-df", type=int, default=4)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    scores = pd.read_csv(args.scores, sep="\t", dtype={"IID": str})
    phenotypes = pd.read_csv(args.phenotypes, sep="\t", dtype={"IID": str})
    summary = pd.read_csv(args.phenotype_summary, sep="\t", dtype={"phenotype_id": str})
    covariates = pd.read_csv(args.covariates, sep="\t", dtype={"IID": str})

    age_cols = [x for x in args.age_cols.split(",") if x]
    pc_cols = [x for x in args.pc_cols.split(",") if x]
    base_cols = age_cols + pc_cols

    df = phenotypes.merge(scores[["IID", args.prs_col]], on="IID", how="inner")
    df = df.merge(covariates[["IID", args.sex_col] + base_cols], on="IID", how="inner")

    rows = []
    for _, meta in summary.iterrows():
        phenotype = meta["phenotype_id"]
        covar_cols = base_cols.copy()
        if str(meta.get("sex_restriction", "both")).lower() == "both":
            covar_cols.append(args.sex_col)

        row = {
            "phenotype_id": phenotype,
            "label": meta.get("label", ""),
            "phecodes": meta.get("phecodes", ""),
            "sex_restriction": meta.get("sex_restriction", ""),
            "category": meta.get("category", ""),
            "n": 0,
            "cases": meta.get("cases", np.nan),
            "controls": meta.get("controls", np.nan),
            "linear_beta": np.nan,
            "linear_or": np.nan,
            "linear_p": np.nan,
            "spline_df": args.spline_df,
            "lrt_df": np.nan,
            "lrt_stat": np.nan,
            "spline_lrt_p": np.nan,
            "status": "not_run",
        }

        if not normalize_bool(meta.get("include", False)) or phenotype not in df.columns:
            row["status"] = meta.get("status", "excluded")
            rows.append(row)
            continue

        model_df = df[[phenotype, args.prs_col] + covar_cols].dropna()
        y = pd.to_numeric(model_df[phenotype], errors="coerce")
        x_linear = prepare_design(model_df, args.prs_col, covar_cols)
        ok = y.notna() & x_linear.notna().all(axis=1)
        y = y.loc[ok].astype(float)
        model_df = model_df.loc[ok].copy()
        x_linear = x_linear.loc[ok].astype(float)
        row["n"] = len(y)
        row["cases"] = int(y.sum())
        row["controls"] = int(len(y) - y.sum())

        try:
            linear_fit = sm.Logit(y, x_linear).fit(disp=False, maxiter=100)
            linear_beta = linear_fit.params[args.prs_col]
            row["linear_beta"] = linear_beta
            row["linear_or"] = math.exp(linear_beta)
            row["linear_p"] = linear_fit.pvalues[args.prs_col]

            spline_basis = patsy.dmatrix(
                f"bs(x, df={args.spline_df}, degree=3, include_intercept=False) - 1",
                {"x": model_df[args.prs_col].astype(float)},
                return_type="dataframe",
            )
            spline_basis.columns = [f"prs_spline_{i + 1}" for i in range(spline_basis.shape[1])]
            covar_design = prepare_design(model_df, args.prs_col, covar_cols).drop(columns=[args.prs_col])
            x_spline = pd.concat([covar_design.reset_index(drop=True), spline_basis.reset_index(drop=True)], axis=1)
            spline_fit = sm.Logit(y.reset_index(drop=True), x_spline.astype(float)).fit(disp=False, maxiter=100)

            lrt_stat = 2 * (spline_fit.llf - linear_fit.llf)
            lrt_df = spline_fit.df_model - linear_fit.df_model
            row["lrt_df"] = lrt_df
            row["lrt_stat"] = max(lrt_stat, 0)
            row["spline_lrt_p"] = st.chi2.sf(row["lrt_stat"], lrt_df)
            row["status"] = "ok"
        except Exception as exc:
            row["status"] = f"failed: {exc}"
        rows.append(row)

    out = pd.DataFrame(rows)
    tested = out["spline_lrt_p"].notna().sum()
    out["spline_bonferroni_p"] = np.where(
        out["spline_lrt_p"].notna(),
        np.minimum(out["spline_lrt_p"] * tested, 1),
        np.nan,
    )
    out["spline_fdr_bh"] = bh_fdr(out["spline_lrt_p"].to_numpy()) if len(out) else []
    out = out.sort_values(["spline_lrt_p", "phenotype_id"], na_position="last")
    out.to_csv(args.out, sep="\t", index=False)
    print(f"tested {tested} nonlinear PRS models")


if __name__ == "__main__":
    main()
