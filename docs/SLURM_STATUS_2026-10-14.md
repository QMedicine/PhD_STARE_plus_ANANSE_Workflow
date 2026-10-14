# Slurm status checked on 14 October 2026

Read-only status check at **20:16 AEDT (UTC+11)**. This update supersedes the live-job observations in [development status v05](DEVELOPMENT_STATUS_v05.md); it does not change its pilot qualification or methods. No jobs were submitted, cancelled, restarted or modified during this check.

## Twelve original High/Low components, pooled into six libraries

| Step | Job | Observed state |
|---|---|---|
| Full raw Hi-C mapping | 30370773 | Completed, exit 0; ended 14 October at 10:22 |
| Raw BAM checks | 30372233, tasks 0–11 | All twelve completed, exit 0 |
| Pooling, counting, STARE and ANANSE | 30381432 | Completed, exit 0; ended at 13:37; all 40 Nextflow tasks completed |
| Result qualification/publication | 30372419 | Failed, exit 1, at 13:41: disk quota exceeded during project copying |
| Dependent reports | 30379823 | Cancelled after the publisher failed; never started |

The completed downstream trace includes all four influence tasks: both directions under both historical and corrected ATAC column ordering. The fixed-grid retained-pair contact route and the separate rounded `squareCounts()` route completed. These jobs used their frozen source snapshots and the existing pooled-TPM converter; they did **not** use the newly documented state-specific TPM correction.

The publisher traceback is an `OSError: [Errno 122] Disk quota exceeded` while creating:

```text
/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs/validated_results/hic_raw_to_stare_v29/diffhic_loess_fixed_grid__historical_column_order/networks/Young_HSC/raw
```

Scratch results remain at `/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/runs/hic_raw_to_stare_v29/results`. The destination may contain a partial copy. Completed computation does not establish completed qualification/publication; no complete final-result package is claimed. Project quota availability must be resolved before copying can succeed. No cleanup or retry was attempted.

## Six independent, already combined paired-end libraries

| Step | Job | Observed state |
|---|---|---|
| Raw mapping driver | 30381427 | Running for 28 h 33 min |
| Last mapping task: `Y_LTHSC_Rep3` | 30400776 | Running for 11 h 26 min; 32 CPUs, 64 GB |
| Full BAM checks | 30381428, tasks 0–5 | Pending dependency |
| Counting, STARE and ANANSE | 30381429 | Pending dependency |
| Qualification/publication | 30381430 | Pending dependency |
| Main and alternate-filter reports | 30381579 | Pending dependency |

Five of six mapping tasks have completed with exit 0 in the current Nextflow trace. Each library is mapped independently without High/Low pooling. The last mapping task is still active; downstream counting and qualification have not yet run for this full cohort. The launcher records END/FAIL email to the institutional address. The pending publisher uses the same durable project destination tree, so project quota also matters to that later copy.

## Evidence inspected

Existing scheduler records were read with `squeue`, `sacct` and `scontrol`. Existing files inspected under `/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow` were:

* `logs/raw_hic_publish_30372419.err` and `.out`.
* `logs/hic_single_library_full_v38/fraction_count_replacement_30381432.out` and `.err`.
* `logs/hic_single_library_full_v38/raw_30381427.out` and `.err`.
* `runs/hic_raw_to_stare_v29/results/provenance/trace_30381432.tsv`: 40 completed tasks, exit 0.
* `runs/hic_single_library_full_raw_v38/results/provenance/trace_30381427.tsv`: one completed validation task and five completed mapping tasks, exit 0; the running sixth mapping task is not yet recorded as a completed row.

This is a timestamped observation, not a continuing monitor or a fresh scientific rerun. Unrelated jobs in the same user's queue were left untouched.
