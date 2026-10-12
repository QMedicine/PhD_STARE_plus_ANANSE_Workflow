# ATAC ordering evidence and two-mode decision v0.1

On 12 October 2026 the user requested preservation of historical behaviour, a separately labelled corrected mode, the snippets for both, and execution of both to determine which reproduces the original results. Both ranking and plotting use H13CORE DESeq2 files.

The current cached `union_counts.tsv` header is Young_rep1, Young_rep2, Old_rep1, Old_rep2. Native audit 30363816 counted the same three saved windows in the four BAMs:

| BAM age/replicate | peak_1 | peak_2 | peak_3 |
|---|---:|---:|---:|
| Young_rep1 | 24 | 26 | 25 |
| Young_rep2 | 26 | 26 | 20 |
| Old_rep1 | 27 | 28 | 27 |
| Old_rep2 | 23 | 23 | 23 |

The cached data rows begin `27,23,24,26`, `28,23,26,26`, `27,23,25,20`. Thus the cached matrix is reversed relative to its header. The saved STARE metadata records 23 January 2026 and consumption of activity columns 7,8 for Old and 5,6 for Young, using this named activity path. This links the recorded STARE invocation to the activity path, but file paths and current content alone do not authenticate historical file bytes. Both modes must be compared with saved affinities and rankings.

Historical snippet (preserved verbatim in the archive):

```bash
LABELS=("Young_rep1" "Young_rep2" "Old_rep1" "Old_rep2")
# The preceding loop writes Young_rep1.counts etc.
paste <(cut -f4 "$BEDWIN") $(find "$TMP_COUNTS_DIR" -name "*.counts" | sort)
```

Alphabetical data order is Old_rep1, Old_rep2, Young_rep1, Young_rep2.

Corrected snippet (the intended file list follows header order):

```bash
count_files=()
for label in "${LABELS[@]}"; do
    count_files+=("$TMP_COUNTS_DIR/${label}.counts")
done
paste <(cut -f4 "$BEDWIN") "${count_files[@]}"
```

The executable adapter is `bin/assemble_atac.py`. It joins filenames to explicit samplesheet identities and writes `column_assignment.json`. In `historical_column_order` the Young-first header and Old-first data are deliberately preserved, only for the two-Old/two-Young fixture. In `corrected_column_order`, each header has its matching sample's data. Intermediate activities, networks, comparisons and reports must retain mode labels; never overwrite one with the other.
