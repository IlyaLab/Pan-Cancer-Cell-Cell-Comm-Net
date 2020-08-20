WITH
  cohort AS (
  SELECT
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
    AND (ajcc_pathologic_tumor_stage = 'Stage I'  OR 
         ajcc_pathologic_tumor_stage = 'Stage IA' OR
         ajcc_pathologic_tumor_stage = 'Stage IB' OR
         ajcc_pathologic_tumor_stage = 'Stage II' OR
         ajcc_pathologic_tumor_stage = 'Stage IIA'   OR
         ajcc_pathologic_tumor_stage = 'Stage IIB'  OR
         ajcc_pathologic_tumor_stage = 'Stage IIC')
    AND b.SampleTypeLetterCode = 'TP'
    AND a.type NOT IN ('DLBC','THYM','LAML')

GROUP BY
    1,2 ),
    --
    --
    --
  edges AS (
  SELECT
    EdgeID,
    EdgeWt,
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
    4),
    --
    --
    --
  medians AS (
  SELECT
    EdgeID,
    COUNT(AliquotBarcode) as n,
    APPROX_QUANTILES(EdgeWt,1000)[
  OFFSET
    (500)] median,
    AVG(EdgeWt) as mean,
    VAR_SAMP(EdgeWt) as sampvar
  FROM
    edges
  GROUP BY
    1 ),
    ---
    ---
    --
  mads AS (
  SELECT
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
  GROUP BY
    EdgeID )
--
--
--
SELECT
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
GROUP BY
  EdgeID,n,median,mean,sampvar,mad  
  