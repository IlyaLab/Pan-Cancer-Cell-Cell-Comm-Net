# getting the gene probability distributions from the cloud bucket.
gsutil -m cp 'gs://pan-cancer-cell-cell-comm-net/GeneProbDist/prob_DC_*' .
gsutil -m cp 'gs://pan-cancer-cell-cell-comm-net/GeneProbDist/prob_CD8+ T-cells*' .