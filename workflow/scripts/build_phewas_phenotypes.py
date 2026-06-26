#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def parse_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def normalize_sex(value):
    value = str(value).strip().lower()
    if value.startswith("f"):
        return "female"
    if value.startswith("m"):
        return "male"
    return "unknown"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", required=True)
    parser.add_argument("--target-map", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--denominator", required=True)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--min-code-count", type=int, default=2)
    parser.add_argument("--force-include-all", action="store_true")
    parser.add_argument("--phenotypes-out", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    targets = pd.read_csv(args.targets, sep="\t", dtype=str)
    target_map = pd.read_csv(args.target_map, sep="\t", dtype=str)
    target_map = target_map[target_map["map_status"].eq("ok")]
    evidence = pd.read_csv(args.evidence, sep="\t", dtype=str)
    denominator = pd.read_csv(args.denominator, sep="\t", dtype={"IID": str})
    scores = pd.read_csv(args.scores, sep="\t", usecols=["IID"], dtype={"IID": str})
    covariates = pd.read_csv(args.covariates, sep="\t", dtype={"IID": str})

    covariates["sex_normalized"] = covariates["sex"].map(normalize_sex)
    base = (
        scores[["IID"]]
        .merge(covariates[["IID", "sex_normalized"]], on="IID", how="inner")
        .merge(denominator[["IID"]].drop_duplicates(), on="IID", how="inner")
        .drop_duplicates("IID")
    )

    mapped = evidence.merge(
        target_map[["phenotype_id", "vocabulary_id", "concept_code"]],
        on=["vocabulary_id", "concept_code"],
        how="inner",
    )
    if not mapped.empty:
        mapped = mapped[mapped["IID"].isin(set(base["IID"]))]
        counts = (
            mapped.dropna(subset=["code_date"])
            .drop_duplicates(["IID", "phenotype_id", "code_date"])
            .groupby(["IID", "phenotype_id"])
            .size()
            .reset_index(name="distinct_code_dates")
        )
        case_pairs = counts[counts["distinct_code_dates"] >= args.min_code_count][["IID", "phenotype_id"]]
    else:
        case_pairs = pd.DataFrame(columns=["IID", "phenotype_id"])

    phenotype_df = base[["IID"]].copy()
    summary_rows = []

    for _, target in targets.iterrows():
        phenotype_id = target["phenotype_id"]
        sex_restriction = str(target.get("sex_restriction", "both")).strip().lower()
        min_cases = int(target.get("min_cases", 0))
        allow_low_cases = parse_bool(target.get("allow_low_cases", "false"))

        eligible = base.copy()
        if sex_restriction in {"female", "male"}:
            eligible = eligible[eligible["sex_normalized"].eq(sex_restriction)]

        cases = set(case_pairs.loc[case_pairs["phenotype_id"].eq(phenotype_id), "IID"])
        eligible_iids = set(eligible["IID"])
        cases = cases & eligible_iids
        n_cases = len(cases)
        n_controls = len(eligible_iids) - n_cases
        include = args.force_include_all or n_cases >= min_cases or allow_low_cases
        status = "included" if include else "excluded_low_cases"
        if args.force_include_all and n_cases < min_cases and not allow_low_cases:
            status = "included_forced_low_cases"
        if target_map[target_map["phenotype_id"].eq(phenotype_id)].empty:
            include = False
            status = "excluded_no_icd_map"

        summary_rows.append({
            "phenotype_id": phenotype_id,
            "label": target["label"],
            "phecodes": target["phecodes"],
            "sex_restriction": sex_restriction,
            "category": target.get("category", ""),
            "min_cases": min_cases,
            "allow_low_cases": allow_low_cases,
            "n_eligible": len(eligible_iids),
            "cases": n_cases,
            "controls": n_controls,
            "include": include,
            "status": status,
        })

        if include:
            values = phenotype_df["IID"].map(lambda iid: 1 if iid in cases else (0 if iid in eligible_iids else pd.NA))
            phenotype_df[phenotype_id] = values.astype("Int64")

    Path(args.phenotypes_out).parent.mkdir(parents=True, exist_ok=True)
    phenotype_df.to_csv(args.phenotypes_out, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(args.summary_out, sep="\t", index=False)
    print(f"wrote phenotype matrix {phenotype_df.shape} to {args.phenotypes_out}")
    print(f"wrote phenotype summary to {args.summary_out}")


if __name__ == "__main__":
    main()
