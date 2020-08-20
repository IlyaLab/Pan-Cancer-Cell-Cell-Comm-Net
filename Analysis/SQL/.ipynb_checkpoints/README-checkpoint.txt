
To generate the S1 statistics, first medians and MADs are computed.
This uses network weights as inputs, which are partitioned by Study.
The median-MAD tables are saved as new BigQuery tables, and are used
to compute the S1 stats; essentially a scaled difference of medians.

For example, for the PFI-study results:

In BigQuery:
1. run short_PFI_study.sql, save results as new BigQuery table.
2. run long_PFI_study.sql, save results as new BigQuery table.
3. use the S1_stats.sql to join the two results tables and compute stats.

Same procedure for early and late Stage sqls.