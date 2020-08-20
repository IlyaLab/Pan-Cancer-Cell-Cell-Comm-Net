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
  (a.median - b.median) / SQRT(POW(a.mad,2) + POW(b.mad,2)) AS S1
FROM
  `isb-cgc-02-0001.Cytokine_Network_Analysis.net_wts_summary_early_stage_by_study`  a
JOIN
  `isb-cgc-02-0001.Cytokine_Network_Analysis.net_wts_summary_late_stage_by_study`  b
ON
  a.EdgeID = b.EdgeID
  AND a.Study = b.Study
WHERE
  (a.mad IS NOT NULL
    AND b.mad IS NOT NULL)
  AND a.mad > 0
  AND b.mad > 0
GROUP BY
1,2,3,4,5,6,7,8
ORDER BY
S1 DESC