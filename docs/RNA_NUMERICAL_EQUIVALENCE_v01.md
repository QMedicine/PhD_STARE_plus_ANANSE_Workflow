# RNA numerical equivalence

The BAM recount plus recovered same-chromosome/strand duplicate cleanup reproduces the selected historical CSV byte for byte (SHA-256 7566126cabd67c2fe55a7275e0654b1005d8f5caecff608f8d344b65ea2d8363). Native Nextflow RNA run rna_recount_v06 completed eight tasks. The selected-ID checkpoint remains explicit; de novo representative selection is not yet qualified.

The archived TPM process reopens counts.clean.tsv. Restoring this boundary removes the observed summation-layout rounding differences. Native boundary regression 30364484 reproduces unfiltered TPM byte for byte (b51e5c1e08368f7e613e6940dc47498e92293d52870ebd267212e01c854a82b3), and filtered raw counts byte for byte (29d771000313d2fd78b90f116cbdeb9c18562c8f1704b4e9c71ff55d0387099c). The intermediate clean counts have identical numeric values, with a different index-column heading.

The archived DE template uses a zero-intercept cell-means design, `~0 + condition`, rather than an intercept model. Native regression 30364483, using the pinned R 4.1.3/DESeq2 1.34.0/ashr 2.2.54 image and Old-first factor levels, reproduces the primary Old→Young DE table byte for byte (b6b36094c8a084355a05b2b334fdfeb78c7607766eeda2b3df9b2e11e2a91a30).

The historical reverse DE file was fitted independently with Young-first factor levels. It is reproduced byte for byte when fitting those levels (job 30364451). It differs from the sign-negated primary fit by at most 2.4110035701241372e-6 across components. Therefore exact byte identity to both independent historical tables and mandatory sign reversal of one fitted table are different contracts.

The user requested sign reversal for the new workflow. Its primary fit uses explicit state ordering, and its reverse output negates only log2FoldChange while retaining all other components. The historical replay route uses the two original directional DE files. Neither comparison should be described as proving that independently refitting both directions has identical numerical bytes.

The corrections are implemented in bin/deseq_historical.R and bin/tpm_filter.py. Their combined updated Nextflow execution and exact resume remain to be qualified separately; native helper regressions do not establish that integration by themselves.
