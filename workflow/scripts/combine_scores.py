#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict


def score_sum_column(fieldnames):
    candidates = [name for name in fieldnames if name.startswith("SCORE") and name.endswith("_SUM")]
    if len(candidates) != 1:
        raise ValueError(f"Expected exactly one *_SUM score column, found {candidates}")
    return candidates[0]


def normalize_row(row):
    if "#IID" in row and "IID" not in row:
        row["IID"] = row.pop("#IID")
    if "#FID" in row and "FID" not in row:
        row["FID"] = row.pop("#FID")
    return row


def numeric(value):
    if value in ("", "NA", "nan"):
        return 0.0
    return float(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-name", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("sscore", nargs="+")
    args = parser.parse_args()

    scores = defaultdict(float)
    denominators = defaultdict(float)
    nalleles = defaultdict(float)
    fids = {}

    for path in args.sscore:
        with open(path, newline="") as score_file:
            reader = csv.DictReader(score_file, delimiter="\t")
            sum_col = score_sum_column(reader.fieldnames)
            for row in reader:
                row = normalize_row(row)
                iid = row["IID"]
                if "FID" in row and row["FID"]:
                    fids[iid] = row["FID"]
                scores[iid] += numeric(row[sum_col])
                if "DENOM" in row:
                    denominators[iid] += numeric(row["DENOM"])
                if "ALLELE_CT" in row:
                    nalleles[iid] += numeric(row["ALLELE_CT"])

    values = list(scores.values())
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1) if len(values) > 1 else 0.0
    sd = variance ** 0.5

    with open(args.out, "w", newline="") as out_file:
        fieldnames = ["FID", "IID", f"{args.score_name}_SUM", f"{args.score_name}_Z", "DENOM", "ALLELE_CT"]
        writer = csv.DictWriter(out_file, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for iid in sorted(scores, key=lambda x: int(x) if x.isdigit() else x):
            score = scores[iid]
            writer.writerow({
                "FID": fids.get(iid, iid),
                "IID": iid,
                f"{args.score_name}_SUM": f"{score:.12g}",
                f"{args.score_name}_Z": f"{((score - mean) / sd):.12g}" if sd else "NA",
                "DENOM": f"{denominators[iid]:.12g}",
                "ALLELE_CT": f"{nalleles[iid]:.12g}",
            })

    print(f"wrote {len(scores)} samples to {args.out}")
    print(f"score_mean\t{mean:.12g}")
    print(f"score_sd\t{sd:.12g}")


if __name__ == "__main__":
    main()
