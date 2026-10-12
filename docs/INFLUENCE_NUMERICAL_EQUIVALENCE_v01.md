# Influence regression: exact ranked scores, explicit file-order limits

The two archived-affinity networks match the historical network files byte for byte. Native replay run network_replay_v09, using frozen source_v11, completed both influence directions (controller 30364743). The local compatibility copy changes only the DE index-column heading to index; all data-body bytes are retained.

Both differential GRNs have exactly the same edge identities and numeric weights as the saved results. All per-TF influence_score, target_score_scaled, G_score, G_score_scaled, direct_targets and factor_fc values match exactly in both directions.

Seven Old→Young TFs and three Young→Old TFs differ in raw target_score, and consequently raw influence_score_raw, by at most 2.220446049250313e-16. These differences do not alter any normalized ranking component. The strict all-component equality test deliberately remains failed; it has not been relaxed or presented as byte-equivalence evidence.

Row-order differences occur solely within tied influence scores: 88 positions in Old→Young and 35 in Young→Old. The score sequences and top-20 TF sets are identical. In particular, NFIC and SATB1 exchange presentation positions within an equal score. Display position inside a tie must not be interpreted as additional biological priority.

ANANSE traverses Python sets and network paths; the historical process hash seed and complete execution state were not recorded. The small raw-score/order differences are consistent with unordered traversal, but that cause is not authenticated. No rounding, numeric backfilling or fabricated component columns has been used to force byte identity.

The original wrapper's -f flag is incompatible with its two-column networks under the observed ANANSE code. Removing it and retaining the GTF option with the index-heading compatibility copy recovers the historical ranked scores. The archived wrapper stays unchanged and the discrepancies remain reviewable in ANANSE_FLAG_INVESTIGATION_v01.md.

Evidence: evidence/audit/influence_regression_v01/influence_regression.json and score_order_detail.json. This closes the normalized-ranking comparison for archived affinity inputs, not the upstream BAM/peak/contact reconstruction or report generation.

Successor saved-union run atac_saved_union_v12 reproduces every per-TF component, including the raw scores, exactly in both directions (qualification 30365485). Ranking bytes still differ; Old→Young top-20 presentation differs within ties, while Young→Old top-20 order matches in that run. This qualifies that execution and does not recover the missing historical hash seed or erase the earlier replay's machine-precision differences. The durable bundle receipt records its individual hashes and comparison scope.
