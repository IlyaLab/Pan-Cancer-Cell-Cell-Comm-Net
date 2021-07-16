#standardSQL
#
# here, the absolute value S1 (sxy) is compared to the
# permuted S1 values to call 'high value edges'

WITH
  nulltab AS(
  SELECT
    string_field_0 Study,
    int64_field_1 early_n,
    int64_field_2 late_n,
    MAX(ABS(double_field_3)) max_S1,
    MAX(ABS(double_field_4)) max_med_diff
  FROM
    `isb-cgc-02-0001.Cytokine_Network_Analysis_v2.study_PFI_S1_null`
  GROUP BY
    1,
    2,
    3 )
SELECT
  a.Study,
  EdgeID,
  med_diff,
  S1
FROM
  `isb-cgc-02-0001.Cytokine_Network_Analysis_v2.study_PFI_S1` AS a
JOIN
  nulltab AS b
ON
  a.Study = b.Study
  AND a.early_n = b.early_n
  AND a.late_n = b.late_n
WHERE
  ABS(a.S1) > b.max_S1
  AND ABS(a.med_diff) > b.max_med_diff
