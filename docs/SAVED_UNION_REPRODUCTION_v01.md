# Reproduction with the preserved summit universe

Native run atac_saved_union_v12 completed 27 Nextflow tasks (controller 30365019), including pre-blacklist BAM filtering/recounting, both explicit column-order modes, corrected-H5 contact export, STARE gABC, the custom converters and four influence comparisons. Reports were disabled explicitly because strict original plotting failures remain unresolved.

The preserved 206,350-window BED restores the old peak universe. Regression job 30365372 compares all four historical-order affinity matrices: identical gene identities, identical 759 numeric columns and zero numeric differences. Both resulting SYMBOLS GRNs also reproduce the original network SHA-256 values byte for byte.

The corrected-column mode produces different gene/TSS records and network hashes. Under the recovered references, contacts and saved peak universe, only the historical column assignment reproduces the saved PhD affinity/GRN outputs. This is stronger evidence than the original script paths or current cached filenames alone. The user may still review their original analysis history; the new pipeline keeps both modes separate.

The first current-summit reconstruction used a smaller, 201,829-window universe and therefore did not reproduce the saved affinities in either mode. It remains a separately labelled development run, rather than being overwritten or conflated with the saved-union reproduction.

Historical-order network hashes: Old_HSC 04b7d275e65551fb93a7d52f157e6a3f62039302513b827fac85302053d8ee94; Young_HSC 31c66fe294329cdc825aa0a6d9fc6c01af185654649f9cbf4c2c79cf82da95eb.

This qualifies BAM-to-STARE/network reconstruction with preserved corrected H5, saved windows and cached RNA TPM/DE checkpoints. Cooler-to-ICE reconstruction, raw Hi-C mapping and fully integrated raw FASTQ routes remain distinct, incomplete checks. Influence byte-level limits are described separately in INFLUENCE_NUMERICAL_EQUIVALENCE_v01.md.
