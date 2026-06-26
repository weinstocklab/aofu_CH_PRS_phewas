-- Template for a simple AoU OMOP condition-concept PheWAS phenotype table.
-- Replace PROJECT.DATASET with $WORKSPACE_CDR, and restrict person IDs to
-- your scored EUR cohort before materializing large outputs.

WITH eur AS (
  SELECT CAST(person_id AS STRING) AS IID
  FROM `PROJECT.DATASET.person`
),
conditions AS (
  SELECT
    CAST(co.person_id AS STRING) AS IID,
    co.condition_concept_id,
    c.concept_name,
    COUNT(*) AS n_events
  FROM `PROJECT.DATASET.condition_occurrence` co
  JOIN `PROJECT.DATASET.concept` c
    ON co.condition_concept_id = c.concept_id
  JOIN eur
    ON CAST(co.person_id AS STRING) = eur.IID
  WHERE co.condition_concept_id != 0
  GROUP BY IID, condition_concept_id, concept_name
)
SELECT *
FROM conditions;
