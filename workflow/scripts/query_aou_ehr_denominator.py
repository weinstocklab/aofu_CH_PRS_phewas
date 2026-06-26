#!/usr/bin/env python3
import argparse
import csv
import os
from pathlib import Path

from google.cloud import bigquery


DENOMINATOR_SQL = """
SELECT DISTINCT CAST(person_id AS STRING) AS IID
FROM (
  SELECT person_id
  FROM `{cdr}.condition_occurrence` co
  JOIN `{cdr}.concept` c
    ON co.condition_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM')
  UNION DISTINCT
  SELECT person_id
  FROM `{cdr}.observation` o
  JOIN `{cdr}.concept` c
    ON o.observation_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM')
  UNION DISTINCT
  SELECT person_id
  FROM `{cdr}.procedure_occurrence` po
  JOIN `{cdr}.concept` c
    ON po.procedure_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM')
  UNION DISTINCT
  SELECT person_id
  FROM `{cdr}.measurement` m
  JOIN `{cdr}.concept` c
    ON m.measurement_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM')
)
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--cdr", default=os.environ.get("WORKSPACE_CDR"))
    parser.add_argument("--billing-project", default=os.environ.get("GOOGLE_PROJECT") or os.environ.get("GOOGLE_CLOUD_PROJECT"))
    args = parser.parse_args()

    if not args.cdr:
        raise ValueError("--cdr or WORKSPACE_CDR is required")
    if not args.billing_project:
        raise ValueError("--billing-project or GOOGLE_PROJECT is required")

    client = bigquery.Client(project=args.billing_project)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = client.query(DENOMINATOR_SQL.format(cdr=args.cdr)).result(page_size=10000)
    n = 0
    with out_path.open("w", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        writer.writerow(["IID"])
        for row in rows:
            writer.writerow([row.IID])
            n += 1
    print(f"wrote {n} EHR denominator IIDs to {args.out}")


if __name__ == "__main__":
    main()
