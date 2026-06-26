#!/usr/bin/env python3
import argparse
import pandas as pd


def split_codes(value):
    return [x.strip() for x in str(value).split(",") if x.strip()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", required=True)
    parser.add_argument("--expanded-map", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    targets = pd.read_csv(args.targets, sep="\t", dtype=str)
    expanded = pd.read_csv(args.expanded_map, sep="\t", dtype=str)

    rows = []
    for _, target in targets.iterrows():
        target_phecodes = set(split_codes(target["phecodes"]))
        subset = expanded[expanded["phecode"].isin(target_phecodes)].copy()
        if subset.empty:
            rows.append({
                "phenotype_id": target["phenotype_id"],
                "label": target["label"],
                "phecode": "",
                "vocabulary_id": "",
                "concept_code": "",
                "map_status": "no_icd_map",
            })
            continue
        subset["phenotype_id"] = target["phenotype_id"]
        subset["label"] = target["label"]
        subset["map_status"] = "ok"
        rows.extend(subset[[
            "phenotype_id",
            "label",
            "phecode",
            "vocabulary_id",
            "concept_code",
            "map_status",
        ]].drop_duplicates().to_dict("records"))

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {len(out)} target ICD map rows to {args.out}")


if __name__ == "__main__":
    main()
