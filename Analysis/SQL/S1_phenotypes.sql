WITH
  cohort AS (
  SELECT
    AliquotBarcode,
    Study,
    SUBSTR(AliquotBarcode, 1, 12) as PatientBarcode
  FROM
    `isb-cgc-02-0001.Cytokine_Network_Work.EBppAliquotBarcodes`
  WHERE
    SampleType = 'Primary solid Tumor'
    AND SUBSTR(AliquotBarcode, 1, 15) IN (
       SELECT SampleBarcode FROM `isb-cgc-02-0001.Cytokine_Network_Work.xCell`)
  GROUP BY
    1,2 ),
--
--
--
stages AS (
  (SELECT
    a.bcr_patient_barcode as PatientBarcode,
    b.Study,
    '2' as Stage
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
    -- AND a.type NOT IN ('DLBC','THYM','LAML'))
    )
--  
  UNION ALL
--
(  SELECT
    a.bcr_patient_barcode as PatientBarcode,
    b.Study as Study,
    '1' as Stage
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
    --AND a.type NOT IN ('DLBC','THYM','LAML')

GROUP BY
    1,2,3) 

),

 PFIMeds AS (
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
    --AND a.type NOT IN ('DLBC',
    -- 'THYM',
    --  'LAML',
    --  'PCPG',
    --  'KICH')
  GROUP BY
    1 ),
  --
  --
  -- then get the samples with PFIs lower than study medians
  --
 pfitable AS (

  (SELECT
    b.Study,
    '2' as PFI,
    a.bcr_patient_barcode as PatientBarcode
  FROM
    `isb-cgc-02-0001.Cytokine_Network_Work.PanCancerSurvival` a
  JOIN
    PFIMeds b
  ON
    a.type = b.Study
  WHERE
    bcr_patient_barcode IN (
    SELECT
      SUBSTR(SampleBarcode, 1, 12)
    FROM
      `isb-cgc-02-0001.Cytokine_Network_Work.xCell`)
    AND a.PFI_time >= b.median_pfi
    --AND a.type NOT IN ('DLBC',
    --  'THYM',
    --  'LAML',
    --  'PCPG',
    --  'KICH')
  GROUP BY
    1,
    2,
    3 )
--
  UNION ALL
--
  (SELECT
    b.Study,
    '1' as PFI,
    a.bcr_patient_barcode as PatientBarcode
  FROM
    `isb-cgc-02-0001.Cytokine_Network_Work.PanCancerSurvival` a
  JOIN
    PFIMeds b
  ON
    a.type = b.Study
  WHERE
    bcr_patient_barcode IN (
    SELECT
      SUBSTR(SampleBarcode, 1, 12)
    FROM
      `isb-cgc-02-0001.Cytokine_Network_Work.xCell`)
    AND a.PFI_time < b.median_pfi
    --AND a.type NOT IN ('DLBC',
    --  'THYM',
    --  'LAML',
    --  'PCPG',
    --  'KICH')
  GROUP BY
    1,
    2,
    3 )
),  

j1 AS (
select 
  c.PatientBarcode,
  c.Study,
  Stage
FROM
cohort c
LEFT JOIN
stages s
ON
s.Study = c.Study
AND s.PatientBarcode = c.PatientBarcode
GROUP BY 
1,2,3
),

j2 AS (
select
    c.Study,
    pfi,
    c.PatientBarcode
FROM
cohort c
LEFT JOIN
pfitable p
ON
c.PatientBarcode = p.PatientBarcode
AND c.Study = p.Study
GROUP BY
1,2,3
)

select
  j2.PatientBarcode,
  j2.Study,
  Stage,
  PFI
FROM
j2 left join j1
ON j2.PatientBarcode = j1.PatientBarcode
AND j2.Study = j1.Study
GROUP BY
1,2,3,4
