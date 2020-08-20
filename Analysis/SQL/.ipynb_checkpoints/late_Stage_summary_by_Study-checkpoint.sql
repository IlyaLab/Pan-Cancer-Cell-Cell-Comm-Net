WITH
  cohort AS (
  SELECT
    b.Study,
    a.bcr_patient_barcode as PatientBarcode,
    b.AliquotBarcode as AliquotBarcode
  FROM
    `isb-cgc-02-0001.Cytokine_Network_Work.PanCancerSurvival` a
  JOIN
    `isb-cgc-02-0001.Cytokine_Network_Work.EBppAliquotBarcodes` b
  ON
     a.bcr_patient_barcode = SUBSTR(b.AliquotBarcode, 1, 12)
  WHERE
     bcr_patient_barcode IN (
        SELECT SUBSTR(SampleBarcode, 1, 12) FROM `isb-cgc-02-0001.Cytokine_Network_Work.xCell`)
    AND (ajcc_pathologic_tumor_stage = 'Stage III'  OR 
         ajcc_pathologic_tumor_stage = 'Stage IIIA' OR
         ajcc_pathologic_tumor_stage = 'Stage IIIB' OR
         ajcc_pathologic_tumor_stage = 'Stage IIIC' OR
         ajcc_pathologic_tumor_stage = 'Stage IV'   OR
         ajcc_pathologic_tumor_stage = 'Stage IVA'  OR
         ajcc_pathologic_tumor_stage = 'Stage IVB'  OR
         ajcc_pathologic_tumor_stage = 'Stage IVC')
    AND b.SampleTypeLetterCode = 'TP'
    AND a.type NOT IN ('DLBC','THYM','LAML')

GROUP BY
    1,2,3 ),
    --
    --
    --
  edges AS (
  SELECT
    EdgeID,
    EdgeWt,
    Study,
    adb.AliquotBarcode as AliquotBarcode,
    SUBSTR(adb.AliquotBarcode, 1, 12) as PatientBarcode
  FROM
    `isb-cgc-02-0001.Cytokine_Network_Work.weighted_cytokine_network` adb  
  JOIN
    cohort
  ON
    cohort.AliquotBarcode = adb.AliquotBarcode
  GROUP BY
    1,
    2,
    3,
    4,
    5),
    --
    --
    --
  medians AS (
  SELECT
    EdgeID,
    Study,
    COUNT(AliquotBarcode) as n,
    APPROX_QUANTILES(EdgeWt,1000)[
  OFFSET
    (500)] median,
    AVG(EdgeWt) as mean,
    VAR_SAMP(EdgeWt) as sampvar
  FROM
    edges
  GROUP BY
    1,2 ),
    ---
    ---
    --
  mads AS (
  SELECT
    edges.Study as Study,
    edges.EdgeID AS EdgeID,
    APPROX_QUANTILES(ABS(EdgeWt - median),1000)[
  OFFSET
    (500)] mad
  FROM
    edges
  JOIN
    medians
  ON
    edges.EdgeID = medians.EdgeID
    AND edges.Study = medians.Study
  GROUP BY
    EdgeID, Study )
--
--
--
SELECT
  medians.Study as Study,
  medians.EdgeID as EdgeID,
  medians.n as n,
  median,
  mad,
  mean,
  sampvar,
  SQRT(sampvar) as stddev
FROM
  medians
JOIN
  mads
ON
  medians.EdgeID = mads.EdgeID
  AND medians.Study = mads.Study
GROUP BY
  Study,EdgeID,n,median,mean,sampvar,mad  
  