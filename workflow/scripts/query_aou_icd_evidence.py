#!/usr/bin/env python3
import argparse
import csv
import os
from pathlib import Path

import pandas as pd
from google.cloud import bigquery


DOMAIN_SQL = """
SELECT person_id, c.vocabulary_id, c.concept_code, condition_start_date AS code_date
FROM `{cdr}.condition_occurrence` co
JOIN `{cdr}.concept` c
  ON co.condition_source_concept_id = c.concept_id
WHERE c.vocabulary_id = @vocabulary_id
  AND c.concept_code IN UNNEST(@concept_codes)
UNION DISTINCT
SELECT person_id, c.vocabulary_id, c.concept_code, observation_date AS code_date
FROM `{cdr}.observation` o
JOIN `{cdr}.concept` c
  ON o.observation_source_concept_id = c.concept_id
WHERE c.vocabulary_id = @vocabulary_id
  AND c.concept_code IN UNNEST(@concept_codes)
UNION DISTINCT
SELECT person_id, c.vocabulary_id, c.concept_code, procedure_date AS code_date
FROM `{cdr}.procedure_occurrence` po
JOIN `{cdr}.concept` c
  ON po.procedure_source_concept_id = c.concept_id
WHERE c.vocabulary_id = @vocabulary_id
  AND c.concept_code IN UNNEST(@concept_codes)
UNION DISTINCT
SELECT person_id, c.vocabulary_id, c.concept_code, measurement_date AS code_date
FROM `{cdr}.measurement` m
JOIN `{cdr}.concept` c
  ON m.measurement_source_concept_id = c.concept_id
WHERE c.vocabulary_id = @vocabulary_id
  AND c.concept_code IN UNNEST(@concept_codes)
"""


def chunks(values, size):
    for i in range(0, len(values), size):
        yield values[i:i + size]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-map", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cdr", default=os.environ.get("WORKSPACE_CDR"))
    parser.add_argument("--billing-project", default=os.environ.get("GOOGLE_PROJECT") or os.environ.get("GOOGLE_CLOUD_PROJECT"))
    parser.add_argument("--batch-size", type=int, default=5000)
    args = parser.parse_args()

    if not args.cdr:
        raise ValueError("--cdr or WORKSPACE_CDR is required")
    if not args.billing_project:
        raise ValueError("--billing-project or GOOGLE_PROJECT is required")

    target_map = pd.read_csv(args.target_map, sep="\t", dtype=str)
    target_map = target_map[target_map["map_status"].eq("ok")]
    codes = (
        target_map[["vocabulary_id", "concept_code"]]
        .drop_duplicates()
        .sort_values(["vocabulary_id", "concept_code"])
    )

    client = bigquery.Client(project=args.billing_project)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    with out_path.open("w", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        writer.writerow(["IID", "vocabulary_id", "concept_code", "code_date"])

        for vocabulary_id, vocab_df in codes.groupby("vocabulary_id"):
            concept_codes = vocab_df["concept_code"].dropna().astype(str).tolist()
            for batch_no, batch in enumerate(chunks(concept_codes, args.batch_size), start=1):
                sql = f"SELECT DISTINCT CAST(person_id AS STRING) AS IID, vocabulary_id, concept_code, CAST(code_date AS STRING) AS code_date FROM ({DOMAIN_SQL.format(cdr=args.cdr)}) WHERE code_date IS NOT NULL"
                job_config = bigquery.QueryJobConfig(
                    query_parameters=[
                        bigquery.ScalarQueryParameter("vocabulary_id", "STRING", vocabulary_id),
                        bigquery.ArrayQueryParameter("concept_codes", "STRING", batch),
                    ]
                )
                print(f"querying {vocabulary_id} batch {batch_no} with {len(batch)} codes")
                rows = client.query(sql, job_config=job_config).result(page_size=10000)
                batch_rows = 0
                for row in rows:
                    writer.writerow([row.IID, row.vocabulary_id, row.concept_code, row.code_date])
                    batch_rows += 1
                total_rows += batch_rows
                print(f"wrote {batch_rows} rows for {vocabulary_id} batch {batch_no}")

    print(f"wrote {total_rows} ICD evidence rows to {args.out}")


if __name__ == "__main__":
    main()
