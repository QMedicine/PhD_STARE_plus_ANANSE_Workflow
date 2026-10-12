# ANANSE full-output flag incompatibility

The provided influence wrapper includes `-f`. In installed ANANSE 0.5.1 this is `--full-output`, which requires tf_expression, target_expression, weighted_binding and activity columns in addition to tf_target/prob. The custom hybrid networks have only tf_target/prob, as the user explicitly specified.

Native influence run 30364644 failed with the missing-column error. The saved historical differential-network headers are source, target and weight. This supports using ordinary probability-only output for the preserved two-column hybrid inputs. The historical wrapper remains unchanged in the archive; the adapter omits -f and records this discrepancy. No missing activity or expression columns are fabricated.

Both reconstructed archived-affinity networks already match their historical SYMBOLS.tsv files byte for byte. A native ranking/differential-network regression without -f is required to establish whether the compatible command recovers the preserved influence results. This evidence is still pending.
