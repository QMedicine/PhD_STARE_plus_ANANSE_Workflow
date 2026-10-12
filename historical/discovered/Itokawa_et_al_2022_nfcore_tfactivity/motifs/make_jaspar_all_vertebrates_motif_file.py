#!/usr/bin/env python3
"""
Build a JASPAR-style motif file for All vertebrates (all TF classes, CORE
collection, all motif versions) from the JASPAR2024 release.

Output: /vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/motifs/JASPAR2024_all_vertebrates_motifs.jaspar
"""

import os
from pyjaspar import jaspardb

OUTPUT = ("/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/motifs/JASPAR2024_all_vertebrates_motifs.jaspar")
# TAX_ID = 10090            # Use this alternatively for just Mus musculus
COLLECTION = "CORE"
RELEASE = "JASPAR2024"

def main() -> None:
    # Connect to the desired JASPAR release
    jdb = jaspardb(release=RELEASE)

    # Retrieve all CORE motifs for all vertebrates, all TF classes, all versions
    motifs = jdb.fetch_motifs(
        collection=COLLECTION,
        tax_group = ['vertebrates'],
        tf_class=None,       # keep all classes
        all_versions=True
    )

    # Ensure destination directory exists
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)

    # Write motifs in classic JASPAR PFM format
    with open(OUTPUT, "w", encoding="utf-8") as fh:
        for mot in motifs:
            fh.write(f">{mot.matrix_id} {mot.name.upper()}\n")
            for base in ("A", "C", "G", "T"):
                counts = " ".join(str(int(x)) for x in mot.counts[base])
                fh.write(f"{base} [ {counts} ]\n")
            fh.write("\n")  # blank line between motifs

if __name__ == "__main__":
    main()