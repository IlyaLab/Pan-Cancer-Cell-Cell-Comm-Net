WITH
--
--
-- First we get the median PFI for each study
--
  PFIs AS (
  SELECT
    b.Study,
    APPROX_QUANTILES(a.PFI_time,1000) [
  OFFSET
    (500)] median_pfi
  FROM
    `isb-cgc-02-0001.Cytokine_Network_Work.PanCancerSurvival` a
  JOIN
    `isb-cgc-02-0001.Cytokine_Network_Work.EBppAliquotBarcodes` b
  ON
    a.bcr_patient_barcode = SUBSTR(b.AliquotBarcode, 1, 12)
  WHERE
    bcr_patient_barcode IN (
    SELECT
      SUBSTR(SampleBarcode, 1, 12)
    FROM
      `isb-cgc-02-0001.Cytokine_Network_Work.xCell`)
    AND a.PFI = 1
    AND b.SampleTypeLetterCode = 'TP'
    AND a.type NOT IN ('DLBC',
      'THYM',
      'LAML',
      'PCPG',
      'KICH')
  GROUP BY
    1 ),
    --
    --
    -- then get the samples with PFIs lower than study medians
    --
  cohort AS (
  SELECT
    b.Study,
    b.median_pfi,
    a.PFI_time,
    a.bcr_patient_barcode
  FROM
    `isb-cgc-02-0001.Cytokine_Network_Work.PanCancerSurvival` a
  JOIN
    PFIs b
  ON
    a.type = b.Study
  WHERE
    bcr_patient_barcode IN (
    SELECT
      SUBSTR(SampleBarcode, 1, 12)
    FROM
      `isb-cgc-02-0001.Cytokine_Network_Work.xCell`)
    AND a.PFI_time >= b.median_pfi
    AND a.type NOT IN ('DLBC',
      'THYM',
      'LAML',
      'PCPG',
      'KICH')
  GROUP BY
    1,
    2,
    3,
    4 ),    
  --
  --
  -- then get edges for samples with PFI as above
  -- AND in the XCell data, which means Primary solid tumor
  --
  edges AS (
  SELECT
    EdgeID,
    EdgeWt,
    Study AS Study,
    adb.AliquotBarcode AS AliquotBarcode,
    bcr_patient_barcode
  FROM
    `isb-cgc-02-0001.Cytokine_Network_Work.net_wts` adb
  JOIN
    cohort
  ON
    cohort.bcr_patient_barcode = SUBSTR(adb.AliquotBarcode, 1, 12)
  WHERE
    SUBSTR(adb.AliquotBarcode, 1, 15) IN (
    SELECT
      SampleBarcode
    FROM
      `isb-cgc-02-0001.Cytokine_Network_Work.xCell`)
  GROUP BY
    1,
    2,
    3,
    4,
    5),
  --
  --
  -- then calculate medians and such
  --
  medians AS (
  SELECT
    EdgeID,
    COUNT(AliquotBarcode) AS n,
    APPROX_QUANTILES(EdgeWt,1000)[
  OFFSET
    (500)] median,
    AVG(EdgeWt) AS mean,
    VAR_SAMP(EdgeWt) AS sampvar
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
  medians.EdgeID AS EdgeID,
  medians.n AS n,
  median,
  mad,
  mean,
  sampvar,
  SQRT(sampvar) AS stddev
FROM
  medians
JOIN
  mads
ON
  medians.EdgeID = mads.EdgeID
GROUP BY
  EdgeID,
  n,
  median,
  mean,
  sampvar,
  mad
