# Repository instructions

This repository reconstructs the PhD method, whose scientific authority is the user request, accepted corrections and `docs/RECONSTRUCTION_PLAN_v01.md`.

Preserve historical scripts verbatim. Put adapters and corrections in separate files. Maintain both explicitly labelled ATAC ordering modes; use H13CORE DE tables for rankings and plots. Do not infer that current cached file contents authenticate historical bytes.

Keep all HPC compute, caches, images, intermediates and development logs beneath `/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow`. The user-authorized durable launcher/manifest and final-result destination is `/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs`. Historical inputs and named environments are read-only.

Use Nextflow for task dependencies and resume, Python/R for calculations, and Bash for cluster setup/launch. No silent host-software fallback. Validate native container operations, sample identities, contrast signs, scores, output integrity and exact resume; scheduler success alone is insufficient.

Retain failures, commands, software/image identities and qualified/deferred routes in compact evidence. Inspect logs before retries. Slurm jobs use END/FAIL mail to the recorded institutional email.

Git contains source, historical scripts, documentation, recipes and compact provenance. Exclude reads, BAMs, contact matrices, images, large scientific outputs, credentials and authentication material. Stage an explicit reviewed allowlist before publication.
