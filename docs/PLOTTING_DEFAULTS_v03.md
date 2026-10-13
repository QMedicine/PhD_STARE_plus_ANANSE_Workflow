# Selected plotting defaults

The 13 October 2026 user decision supersedes the earlier strict empty-plot policy. Historical notebooks and earlier evidence remain unchanged.

The default focal plots derive from `GRN_plotting_Young2Old.Rmd` lines 1997–2554 and `GRN_plotting_Old2Young.Rmd` lines 14–571. Both use `layout_nicely`, seed 1, an outward two-hop graph, path2 edge widths, square TFs, circular genes, a directional logFC gradient and the scaled-score/unfiltered-neighbourhood panel. Hop distances in the information panel are unweighted.

| Transition | Main edge cutoff | Main adjusted P cutoff | Main absolute logFC cutoff | Alternate edge cutoff | Alternate absolute logFC cutoff |
|---|---:|---:|---:|---:|---:|
| Young_HSC → Old_HSC | 0.01 | 0.05 | 0.25 | 0.001 | 0.1 |
| Old_HSC → Young_HSC | 0.0001 | 0.05 | 0.2 | 0.00001 | 0.1 |

The alternate thresholds are an explicit development choice: tenfold lower edge cutoffs, absolute logFC 0.1, and the same adjusted P cutoff. They change display only. Both variants reuse exactly the same ranking, component scores and directional H13CORE DE table. They are published into separate folders; an alternate folder's name includes its edge, adjusted P and absolute logFC cutoffs. The main report contains no alternate panels.

Known no-edge conditions produce a centred TF square coloured by its measured directional logFC and the exact requested warning: `Warning, <TF> edge values fall below the minimum cutoff value of <min_edge_weight_plot>`. The actual reason, including a DE or outward-edge filter removing the last edge, is retained in `plot_status.tsv`. Unexpected errors remain fatal. Missing DE values are grey and are not converted to zero.

Each report embeds exactly 20 focal PDF panels and one top-10 overview PDF. The overview derives from the first Young2Old notebook section, lines 13–471, retaining its top-10 influence selection and up to 15 non-TF and 15 TF targets per focal TF. Its missing `width_raw` attribute is repaired. Added point layers make every retained node visible, including unlabelled genes. Edges are black. TFs and genes use their actual directional logFC regardless of DE significance; labels are translucent and displaced. Overview edge cutoffs follow the report variant. The overview does not apply the focal plot's DE node-exclusion rule.

Generic workflow defaults are edge 0.01, adjusted P 0.05, absolute logFC 0.25 and logFC limits −2 to 2. Study-specific differences reside in the `plot_settings` map keyed by exact `source->target` states. Engine scripts contain no study labels. `make_relaxed_reports`, `relaxed_plot_edge_factor` and `relaxed_plot_abs_lfc_cutoff` control alternate reports.

`report_replay.nf` renders existing explicitly named ranking/differential-network/DE tables without rerunning scientific inference. Its TSV manifest has `mode`, `source`, `target`, `ranking`, `diffnetwork`, and `de` columns. This allows report updates without disturbing active Hi-C mapping or downstream inference jobs. Plot scripts, node/edge exports, settings and statuses accompany the PDFs; HTML packages contain the PDF bytes and require no remote viewer assets.

Template identities and all adapters are recorded in `PLOT_TEMPLATE_PROVENANCE_v03.json`.
