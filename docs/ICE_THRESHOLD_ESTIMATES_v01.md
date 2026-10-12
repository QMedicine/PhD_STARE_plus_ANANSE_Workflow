# Exploratory ICE threshold estimates

The original threshold TSV was in expired scratch and has not been recovered. The user authorized diagnostic-based estimates on 12 October 2026. These estimates do not replace the authenticated corrected-H5 checkpoint route when testing exact STARE reproduction.

Native diagnostic job 30364457 converted the two preserved merged 100-kb coolers, normalized together with hicNormalize --normalize smallest, and generated HiCExplorer 3.7.6 modified-z-score diagnostic plots. The original blacklist step before the supplied conversion scripts remains unresolved; these plots therefore describe a pre-blacklist exploratory reconstruction.

Visual inspection shows a small near-zero-count population separated from the principal broad peak, centred around modified z=0 to 1. As an initial exploratory filter, use lower=-3 and upper=3 for both states. This cuts the low-coverage shoulder and the sparse far-right tail while retaining the main peak. It is a visual estimate, not a recovered historical setting or an optimized threshold.

The software's automatic lower-line values were -5.326644732771318 and -5.247589020217729; these are diagnostic outputs and do not prove the user originally selected them. Estimated correction outputs must carry a reconstructed/estimated status and remain separate from historical corrected matrices.

Scratch diagnostic files: runs/hic_diagnostics_v01_30364457/Old_HSC.diagnostic.png and Young_HSC.diagnostic.png. Image hashes and input manifest are recorded in that run's checksums.sha256. No diagnostic image is committed to Git.
