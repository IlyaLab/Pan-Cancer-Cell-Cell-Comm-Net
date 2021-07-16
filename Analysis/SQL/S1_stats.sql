# saved query: S1 stat table
# note: uses the v2 tables.
SELECT
  a.Study,
  a.EdgeID,
  a.n AS early_n,
  a.median AS early_med,
  a.mad AS early_mad,
  b.n AS late_n,
  b.median AS late_med,
  b.mad AS late_mad,
  a.median - b.median AS med_diff,
  (a.median - b.median) / SQRT( 1.4826*a.mad  + 1.4826*b.mad ) AS S1
FROM
  `isb-cgc-02-0001.Cytokine_Network_Analysis_v2.summary_by_study` a
JOIN
  `isb-cgc-02-0001.Cytokine_Network_Analysis_v2.summary_by_NOT_study` b
ON
  a.EdgeID = b.EdgeID and a.Study = b.Study
WHERE
  (a.mad IS NOT NULL
    AND b.mad IS NOT NULL)
  AND (a.mad > 0 OR b.mad > 0)
GROUP BY
1,2,3,4,5,6,7,8
