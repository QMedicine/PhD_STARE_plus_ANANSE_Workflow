# Third-party code and attribution

Historical scripts supplied by the user are archived verbatim with HPC paths and SHA-256 manifests. They remain evidence rather than generic workflow code.

The public diffHic mapper is archived unchanged at historical/upstream_diffhic_standard with its source commit, checksums, package DESCRIPTION and GPL-3 text. Its separately maintained task adapter documents compatibility changes; see HIC_STANDARD_PRESPLIT_AND_ARIMA_v01.md. No claim is made that this recovers the missing custom _t32 helper.

The tfactivity source templates and support modules are from nf-core/tfactivity revision a6d7fec1ebf00ddabaaa347ed9f6025626da2697. Its copied MIT licence is retained with the archive. The TPM/filter and differential-expression adapters document their historical source and preserve the recovered scientific defaults.

The installed ANANSE 0.5.1 Python source snapshot is from the user's named environment. Its installed MIT licence is preserved at historical/software/ananse/LICENSE. The official project is https://github.com/vanheeringen-lab/ANANSE and its release record is https://github.com/vanheeringen-lab/ANANSE/releases/tag/v0.5.1. The source snapshot serves as an audit of actual influence scoring, and the workflow executes the pinned container instead of importing the archived snapshot.

STARE modified_STARE.sh is the supplied locally modified historical shell script. The original STARE 1.0.5 image remains an HPC input with its recorded digest; no compiled software or container image is committed to Git.

References: Hecker et al. (2023), https://doi.org/10.1093/bioinformatics/btad062; Xu et al. (2021), https://doi.org/10.1093/nar/gkab598. The customized STARE-to-ANANSE converter is the user's PhD hybrid step, rather than ANANSE's binding/network implementation.
