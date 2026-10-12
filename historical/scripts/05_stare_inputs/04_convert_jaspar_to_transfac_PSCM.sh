#!/bin/bash
#SBATCH --job-name=jaspar2transfac
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:60:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/04_convert_jaspar_to_transfac_PSCM_%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/04_convert_jaspar_to_transfac_PSCM_%x_%j.err

set -euo pipefail

# --------------------------
# Inputs / outputs (edit as needed)
# --------------------------
JASPAR_IN="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/motifs.jaspar"
OUTDIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/motifs"
OUTFILE="${OUTDIR}/motifs.transfac"     # pass this to STARE with: -s motifs.transfac

# Optional: set to "true" if you want to replace spaces and odd chars in names
SANITIZE="false"

mkdir -p "$OUTDIR"

# Export so Rscript can see them
export JASPAR_IN OUTDIR OUTFILE SANITIZE

# --------------------------
# Modules
# --------------------------
module load R/4.5.1

# --------------------------
# R: JASPAR (PFM) -> TRANSFAC (PSCM)
# --------------------------
Rscript - <<'RSCRIPT'
infile   <- Sys.getenv("JASPAR_IN", "")
outdir   <- Sys.getenv("OUTDIR", "")
outfile  <- Sys.getenv("OUTFILE", "")
sanitize <- tolower(Sys.getenv("SANITIZE", "false")) %in% c("1","true","yes")

if (infile == "" || outfile == "") {
  stop("Environment variables JASPAR_IN and OUTFILE must be set.")
}
if (!file.exists(infile)) stop("Input JASPAR file not found: ", infile)
if (outdir != "" && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("universalmotif", quietly = TRUE)) {
  stop("The 'universalmotif' package is not available in this R environment.")
}
library(universalmotif)

message("[*] Reading JASPAR: ", infile)
motifs <- read_jaspar(infile)  # PFM counts -> list of universalmotif objects
if (length(motifs) == 0) stop("No motifs were read from: ", infile)

# Make names unique (only adds suffixes if needed). No changes otherwise.
orig_names <- vapply(motifs, function(m) m@name, character(1))
new_names  <- orig_names

# Optional light sanitizing (if you choose to enable it)
if (sanitize) {
  clean_name <- function(x) {
    x <- gsub("\\s+", "_", x)            # spaces -> underscores
    x <- gsub("[^A-Za-z0-9_:+.-]", "", x) # drop odd chars that can annoy CLIs
    x
  }
  new_names <- vapply(motifs, function(m) clean_name(m@name), character(1))
}

# Always ensure uniqueness (TRANSFAC prefers unique IDs; downstream matrices too)
uniq_names <- make.unique(new_names, sep = "_")

if (!all(uniq_names == orig_names)) {
  changed <- sum(uniq_names != orig_names)
  message(sprintf("[i] Adjusted %d motif name(s) (%s sanitizing).",
                  changed, if (sanitize) "with" else "without"))
}

for (i in seq_along(motifs)) motifs[[i]]@name <- uniq_names[i]

message("[*] Writing TRANSFAC PSCM: ", outfile)
write_transfac(motifs, file = outfile, overwrite = TRUE)
message("[✓] Wrote ", length(motifs), " motifs to TRANSFAC.")
RSCRIPT
