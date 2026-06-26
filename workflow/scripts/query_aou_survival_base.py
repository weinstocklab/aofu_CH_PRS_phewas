#!/usr/bin/env python3
import argparse
import csv
import os
from pathlib import Path

from google.cloud import bigquery


SURVIVAL_SQL = """
WITH obs AS (
  SELECT
    person_id,
    MIN(observation_period_start_date) AS observation_start_date,
    MAX(observation_period_end_date) AS observation_end_date
  FROM `{cdr}.observation_period`
  GROUP BY person_id
)
SELECT
  CAST(p.person_id AS STRING) AS IID,
  CAST(DATE(p.birth_datetime) AS STRING) AS birth_date,
  CAST(obs.observation_start_date AS STRING) AS observation_start_date,
  CAST(obs.observation_end_date AS STRING) AS observation_end_date,
  CAST(d.death_date AS STRING) AS death_date
FROM `{cdr}.person` p
JOIN obs
  ON p.person_id = obs.person_id
LEFT JOIN `{cdr}.death` d
  ON p.person_id = d.person_id
WHERE p.birth_datetime IS NOT NULL
  AND obs.observation_start_date IS NOT NULL
  AND obs.observation_end_date IS NOT NULL
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
    rows = client.query(SURVIVAL_SQL.format(cdr=args.cdr)).result(page_size=10000)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out_path.open("w", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        writer.writerow(["IID", "birth_date", "observation_start_date", "observation_end_date", "death_date"])
        for row in rows:
            writer.writerow([row.IID, row.birth_date, row.observation_start_date, row.observation_end_date, row.death_date or ""])
            n += 1
    print(f"wrote {n} survival base rows to {args.out}")


if __name__ == "__main__":
    main()
