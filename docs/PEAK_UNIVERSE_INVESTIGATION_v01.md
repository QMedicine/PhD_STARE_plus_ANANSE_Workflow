# Peak-universe drift

The saved PhD union has 206,350 unmerged 501-bp summit windows. Rebuilding the union from the four currently stored summit files yields 201,829 windows, a difference of 4,521. Both subsequent ATAC column-order modes therefore used a different peak universe from the saved STARE run; neither comparison alone can resolve which ordering produced the historical results.

The historical all_summits.1bp.bed, all_summits.500bp.clean.bed, union_base.tsv, union_counts.tsv and union_summit500.named.bed remain intact in the project union_bed directory. They each describe the 206,350-window checkpoint. The current summit files may have been overwritten after an earlier peak-calling run; that provenance is not yet established.

`--union_windows` accepts the explicitly named preserved four-column BED. With it, the workflow recounts the pre-blacklist BAMs after the original blacklist filtering and evaluates both column assignments on the same historical windows. Counts/offsets and STARE affinities must be compared before selecting a reproduced ordering. Without this parameter, the workflow builds windows from the supplied peaks and records a distinct reconstruction.

The native reconstructed matrices had fewer gene/TSS records than the saved affinity matrices, while motif columns matched. This is evidence of a mismatch, not evidence that either historical ordering has been authenticated. The saved-union comparison remains pending.
