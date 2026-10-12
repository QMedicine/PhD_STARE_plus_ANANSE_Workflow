#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(scales)
  library(units)
  library(Cairo)
  library(grid)
  library(patchwork)
})

# ============================================================
# 2-step branching subnet plot for ANANSE influence intuition
#
# - Uses UNFILTERED diffnetwork for TABLE counts (true ANANSE diffnetwork)
# - Uses FILTERED diffnetwork for PLOT (min_edge_weight)
# - Builds outward ego graph (radius=max_steps) around tf_of_interest
# - Optional outward-edge filter for plotting graph (layer_u < layer_v)
# - TF nodes squares; gene nodes circles
# - Node fill is directional log2FC (blue-white-red)
# - Edge width is either direct edge or 2-step proxy (path2)
# - Adds:
#     (1) edge-width legend (via scale_edge_width)
#     (2) right-side info panel (outside plot border) BELOW the log2FC legend
#         showing scaled ANANSE scores + true 2-step neighbourhood counts
#
# Notes on scores:
# - ANANSE influence_score_raw = target_score + G_score (raw)
# - ANANSE influence_score (scaled 0-1) comes from rank/minmax on
#   target_score_scaled + G_score_scaled.
#   In this script, we display the scaled values as requested.
# ============================================================

# ----------------------------
# Input files
# ----------------------------
influence_file <- "" # supplied by the workflow adapter

diffnet_file <- "" # supplied by the workflow adapter

deseq_file <- "" # supplied by the workflow adapter


# ----------------------------
# User knobs
# ----------------------------
tf_of_interest <- "" # supplied by the workflow adapter
max_steps <- 2

# PLOT edge filter only (table uses the unfiltered diffnetwork)
min_edge_weight <- 0.01

# DE definition (used for node colouring; also used for "DE genes reachable (strict)" count)
padj_cutoff    <- 0.05
abs_lfc_cutoff <- 0.25

# Densify plot (optional)
add_top_nonDE_genes <- FALSE
n_top_nonDE_genes   <- 40

include_unmatched_genes <- FALSE
n_unmatched_genes       <- 25

# Branching structure
keep_outward_edges <- TRUE

# Edge width definition:
#  "edge"  = direct edge weight
#  "path2" = for non-TF1 edges, width_raw ~ w(TF1->source) * w(source->target) if TF1->source exists
edge_width_mode <- "path2"

# Aesthetics
tf_shape    <- 22
gene_shape  <- 21
node_border <- "grey85"

edge_col   <- "grey10"
edge_alpha <- 1
edge_width_range <- c(0.4, 2.2)
arrow_mm   <- 2.3

nonDE_fill     <- "grey70"
unmatched_fill <- "grey85"

fc_low  <- "blue3"
fc_mid  <- "white"
fc_high <- "red3"

# Labels
label_tfs <- TRUE
label_top_de_genes <- TRUE
n_label_de_genes <- 25

label_nudge_x <- 0.00
label_nudge_y <- 0.1
label_hjust   <- 0.5
label_vjust   <- -1
label_box_padding_lines <- 0.12
label_box_radius_lines  <- 0.12
label_fill_col   <- "white"
label_fill_alpha <- 0.50
label_text_col   <- "black"
label_border_size <- 0.0

# Output
out_dir <- "" # supplied by the workflow adapter
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- "" # supplied by the workflow adapter
out_width <- 12
out_height <- 9


# ----------------------------
# Helpers
# ----------------------------
stop_if_missing <- function(path) if (!file.exists(path)) stop("File not found: ", path)

msg_counts <- function(label, g) {
  message(label, ": nodes=", vcount(g), " edges=", ecount(g))
}

# ----------------------------
# Load
# ----------------------------
stop_if_missing(influence_file)
stop_if_missing(diffnet_file)
stop_if_missing(deseq_file)

inf <- read_tsv(influence_file, show_col_types = FALSE, na = c("", "NA"))
if (!("factor" %in% colnames(inf))) stop("Influence file missing column 'factor'.")

diffnet_full <- read_tsv(diffnet_file, show_col_types = FALSE, na = c("", "NA")) %>%
  select(source, target, weight) %>%
  filter(!is.na(source), !is.na(target), !is.na(weight))

if (nrow(diffnet_full) == 0) stop("diffnet_full is empty after reading.")

# filtered copy for plotting only
diffnet_plot <- diffnet_full %>% filter(weight >= min_edge_weight)
if (nrow(diffnet_plot) == 0) stop("No edges after min_edge_weight filtering; lower min_edge_weight.")

# TF universe: all scored TFs + all sources in the diffnetwork
tf_universe <- union(unique(inf$factor), unique(diffnet_full$source))
if (!(tf_of_interest %in% tf_universe)) stop("tf_of_interest not found in TF universe: ", tf_of_interest)

de <- read_tsv(deseq_file, show_col_types = FALSE, na = c("", "NA")) %>%
  select(gene_id, log2FoldChange, padj)

de2 <- de %>%
  mutate(
    abs_lfc = abs(log2FoldChange),
    is_DE_strict = !is.na(padj) & !is.na(log2FoldChange) & (padj <= padj_cutoff) & (abs_lfc >= abs_lfc_cutoff),
    is_DE_padj   = !is.na(padj) & (padj <= padj_cutoff)
  )

# ----------------------------
# TABLE metrics from TRUE ANANSE diffnetwork (unfiltered)
# ----------------------------
g_all_table <- graph_from_data_frame(diffnet_full, directed = TRUE)
E(g_all_table)$weight <- diffnet_full$weight

if (!(tf_of_interest %in% V(g_all_table)$name)) {
  stop("TF not found in unfiltered diffnetwork graph (table graph): ", tf_of_interest)
}

# Ego graph within max_steps on the true network
ego_tbl_vids <- ego(g_all_table, order = max_steps, nodes = tf_of_interest, mode = "out")[[1]]
g_tbl0 <- induced_subgraph(g_all_table, vids = ego_tbl_vids)

# hop distances (unweighted) for layer counts
dmat_tbl <- distances(g_tbl0, v = tf_of_interest, to = V(g_tbl0), mode = "out")
layer_tbl <- as.numeric(dmat_tbl[1, ])
names(layer_tbl) <- colnames(dmat_tbl)

# total primary/secondary nodes (all node types)
tbl_direct_targets_total    <- sum(layer_tbl == 1, na.rm = TRUE)
tbl_secondary_targets_total <- sum(layer_tbl == 2, na.rm = TRUE)
tbl_targets_within2_total   <- sum(is.finite(layer_tbl) & layer_tbl >= 1 & layer_tbl <= max_steps, na.rm = TRUE)

# DE genes reachable (strict and padj-only), based on DE table overlap with ego nodes
ego_tbl_nodes <- names(layer_tbl)[is.finite(layer_tbl) & layer_tbl >= 1 & layer_tbl <= max_steps]

tbl_de_strict <- de2 %>% filter(gene_id %in% ego_tbl_nodes, is_DE_strict) %>% nrow()
tbl_de_padj   <- de2 %>% filter(gene_id %in% ego_tbl_nodes, is_DE_padj)   %>% nrow()

# Pull ANANSE scores for the TF (scaled)
tf_inf <- inf %>% filter(factor == tf_of_interest) %>% slice(1)
if (nrow(tf_inf) == 0) stop("TF not found in influence table: ", tf_of_interest)

# ----------------------------
# PLOT graph from FILTERED diffnetwork
# ----------------------------
g_all <- graph_from_data_frame(diffnet_plot, directed = TRUE)
E(g_all)$weight <- diffnet_plot$weight
msg_counts("g_all (plot graph, filtered)", g_all)

if (!(tf_of_interest %in% V(g_all)$name)) {
  stop("TF not found as a node in filtered plot graph. Lower min_edge_weight or use unfiltered for plotting.")
}

# Ego graph (outward) within max_steps on plot graph
ego_vids <- ego(g_all, order = max_steps, nodes = tf_of_interest, mode = "out")[[1]]
g_sub0 <- induced_subgraph(g_all, vids = ego_vids)
msg_counts("g_sub0 (plot ego induced)", g_sub0)

if (ecount(g_sub0) == 0) stop("g_sub0 has 0 edges; lower min_edge_weight.")

# Layers with safe name mapping (plot graph)
dmat <- distances(g_sub0, v = tf_of_interest, to = V(g_sub0), mode = "out")
layer_vec <- as.numeric(dmat[1, ])
names(layer_vec) <- colnames(dmat)
V(g_sub0)$layer <- layer_vec[V(g_sub0)$name]

keep_nodes <- V(g_sub0)$name[is.finite(V(g_sub0)$layer) & V(g_sub0)$layer <= max_steps]
g_sub <- induced_subgraph(g_sub0, vids = keep_nodes)
V(g_sub)$layer <- layer_vec[V(g_sub)$name]

message("Layer counts in g_sub (plot graph):")
print(table(V(g_sub)$layer, useNA = "ifany"))

# Outward-edge filter for plot structure
if (isTRUE(keep_outward_edges)) {
  edf <- as_data_frame(g_sub, what = "edges") %>%
    rename(source = from, target = to) %>%
    mutate(
      layer_u = layer_vec[source],
      layer_v = layer_vec[target]
    )

  edf2 <- edf %>%
    filter(is.finite(layer_u), is.finite(layer_v),
           layer_u < layer_v,
           layer_v <= max_steps) %>%
    select(source, target, weight)

  if (nrow(edf2) == 0) stop("Outward-edge filter removed all edges. Set keep_outward_edges=FALSE or lower min_edge_weight.")

  g_sub <- graph_from_data_frame(
    edf2,
    directed = TRUE,
    vertices = data.frame(name = unique(c(edf2$source, edf2$target)))
  )
  V(g_sub)$layer <- layer_vec[V(g_sub)$name]
  E(g_sub)$weight <- edf2$weight
  msg_counts("g_sub (after outward-edge filter)", g_sub)
}

# ----------------------------
# Node annotations (plot graph)
# ----------------------------
nodes_df <- tibble(name = V(g_sub)$name) %>%
  mutate(
    layer = V(g_sub)$layer,
    is_center = name == tf_of_interest,
    is_tf = name %in% tf_universe,
    node_type = if_else(is_tf, "TF", "Gene")
  ) %>%
  left_join(de2, by = c("name" = "gene_id")) %>%
  mutate(
    has_deseq = !is.na(log2FoldChange) & !is.na(padj),
    is_DE = dplyr::coalesce(is_DE_strict, FALSE),
    log2FoldChange = if_else(is_center & !has_deseq, 0, log2FoldChange),
    has_deseq = if_else(is_center, TRUE, has_deseq)
  )

edges_df <- as_data_frame(g_sub, what = "edges") %>%
  rename(source = from, target = to) %>%
  mutate(abs_weight = abs(weight))

incoming_gene <- edges_df %>%
  group_by(target) %>%
  summarise(max_incoming = max(abs_weight, na.rm = TRUE), .groups = "drop") %>%
  mutate(max_incoming = ifelse(is.infinite(max_incoming), 0, max_incoming))

nodes_df <- nodes_df %>%
  left_join(incoming_gene, by = c("name" = "target")) %>%
  mutate(max_incoming = if_else(is.na(max_incoming), 0, max_incoming))

extra_nonDE_genes <- character(0)
if (isTRUE(add_top_nonDE_genes)) {
  extra_nonDE_genes <- nodes_df %>%
    filter(node_type == "Gene", has_deseq, !is_DE) %>%
    arrange(desc(max_incoming)) %>%
    slice_head(n = n_top_nonDE_genes) %>%
    pull(name)
}

extra_unmatched_genes <- character(0)
if (isTRUE(include_unmatched_genes)) {
  extra_unmatched_genes <- nodes_df %>%
    filter(node_type == "Gene", !has_deseq) %>%
    arrange(desc(max_incoming)) %>%
    slice_head(n = n_unmatched_genes) %>%
    pull(name)
}

keep_final_nodes <- nodes_df %>%
  filter(
    is_tf |
      (node_type == "Gene" & is_DE) |
      (node_type == "Gene" & name %in% extra_nonDE_genes) |
      (node_type == "Gene" & name %in% extra_unmatched_genes)
  ) %>%
  pull(name) %>%
  unique()

g_final <- induced_subgraph(g_sub, vids = keep_final_nodes)
msg_counts("g_final (plot, after keeping TFs + DE + extras)", g_final)
if (ecount(g_final) == 0) stop("g_final has 0 edges; relax filters or disable outward edges.")

# Attach vertex attrs
nodes_df_final <- nodes_df %>% filter(name %in% V(g_final)$name)

V(g_final)$layer          <- nodes_df_final$layer[match(V(g_final)$name, nodes_df_final$name)]
V(g_final)$is_center      <- nodes_df_final$is_center[match(V(g_final)$name, nodes_df_final$name)]
V(g_final)$node_type      <- nodes_df_final$node_type[match(V(g_final)$name, nodes_df_final$name)]
V(g_final)$has_deseq      <- nodes_df_final$has_deseq[match(V(g_final)$name, nodes_df_final$name)]
V(g_final)$is_DE          <- nodes_df_final$is_DE[match(V(g_final)$name, nodes_df_final$name)]
V(g_final)$log2FoldChange <- nodes_df_final$log2FoldChange[match(V(g_final)$name, nodes_df_final$name)]
V(g_final)$abs_lfc        <- nodes_df_final$abs_lfc[match(V(g_final)$name, nodes_df_final$name)]

V(g_final)$shape <- ifelse(V(g_final)$node_type == "TF", tf_shape, gene_shape)
V(g_final)$size  <- ifelse(V(g_final)$is_center, 13,
                           ifelse(V(g_final)$node_type == "TF", 9, 7))

# Labels: TFs + top DE genes by abs_lfc
plot_label <- rep(NA_character_, vcount(g_final))
names(plot_label) <- V(g_final)$name

if (isTRUE(label_tfs)) {
  tf_names <- V(g_final)$name[V(g_final)$node_type == "TF"]
  plot_label[tf_names] <- tf_names
}
if (isTRUE(label_top_de_genes)) {
  de_gene_names <- V(g_final)$name[V(g_final)$node_type == "Gene" & V(g_final)$is_DE %in% TRUE]
  abs_lfc_vec <- V(g_final)$abs_lfc
  abs_lfc_vec[is.na(abs_lfc_vec)] <- 0
  ord <- order(abs_lfc_vec[match(de_gene_names, V(g_final)$name)], decreasing = TRUE)
  top_genes <- head(de_gene_names[ord], n_label_de_genes)
  plot_label[top_genes] <- top_genes
}
V(g_final)$plot_label <- unname(plot_label[V(g_final)$name])

# Edge widths
edf_final <- as_data_frame(g_final, what = "edges") %>%
  rename(source = from, target = to) %>%
  mutate(abs_weight = abs(weight))

# TF1->X lookup
tf1_out <- edf_final %>%
  filter(source == tf_of_interest) %>%
  select(target, w_tf1 = abs_weight) %>%
  distinct()

if (tolower(edge_width_mode) == "path2") {
  edf_final <- edf_final %>%
    left_join(tf1_out, by = c("source" = "target")) %>%
    mutate(
      width_raw = case_when(
        source == tf_of_interest ~ abs_weight,
        !is.na(w_tf1)            ~ abs_weight * w_tf1,
        TRUE                     ~ abs_weight
      ),
      width_raw = ifelse(is.na(width_raw), abs_weight, width_raw)
    )
} else {
  edf_final <- edf_final %>% mutate(width_raw = abs_weight)
}

edf_final <- edf_final %>% mutate(edge_width = rescale(width_raw, to = edge_width_range))

# Rebuild final plot graph with edge attributes + vertex attributes
vdf <- tibble(
  name           = V(g_final)$name,
  layer          = V(g_final)$layer,
  is_center      = V(g_final)$is_center,
  node_type      = V(g_final)$node_type,
  has_deseq      = V(g_final)$has_deseq,
  is_DE          = V(g_final)$is_DE,
  log2FoldChange = V(g_final)$log2FoldChange,
  abs_lfc        = V(g_final)$abs_lfc,
  shape          = V(g_final)$shape,
  size           = V(g_final)$size,
  plot_label     = V(g_final)$plot_label
)

g_final2 <- graph_from_data_frame(
  d = edf_final %>% transmute(from = source, to = target, weight, width_raw),
  directed = TRUE,
  vertices = vdf
)

# Layout
layers_vec <- V(g_final2)$layer
names(layers_vec) <- V(g_final2)$name

coords <- tryCatch({
  sug <- layout_with_sugiyama(g_final2, layers = layers_vec)
  sug$layout
}, error = function(e) {
  message("layout_with_sugiyama failed (", conditionMessage(e), "). Falling back to layout_with_fr().")
  layout_with_fr(g_final2)
})

center_idx <- which(V(g_final2)$name == tf_of_interest)
coords <- sweep(coords, 2, coords[center_idx, ], FUN = "-")

# Drop isolated nodes
iso <- V(g_final2)[degree(g_final2, mode = "all") == 0]
if (length(iso) > 0) g_final2 <- delete_vertices(g_final2, iso)

# ----------------------------
# Build RIGHT-SIDE TABLE panel (outside plot border)
# (Use scaled values as requested)
# ----------------------------
summary_lines <- c(
  paste0("TF: ", tf_of_interest),
  paste0("Direct targets (1-step, total): ", tbl_direct_targets_total),
  paste0("Secondary targets (2-step, total): ", tbl_secondary_targets_total),
  paste0("Targets within ", max_steps, " steps (total): ", tbl_targets_within2_total),
  paste0("DE genes reachable (padj & absLFC≥", abs_lfc_cutoff, "): ", tbl_de_strict),
  paste0("DE genes reachable (padj only): ", tbl_de_padj),
  "",
  paste0("target_score_scaled: ", signif(tf_inf$target_score_scaled, 4)),
  paste0("G_score_scaled: ", signif(tf_inf$G_score_scaled, 4)),
  paste0("influence_score (scaled): ", signif(tf_inf$influence_score, 4))
)

table_grob <- grobTree(
  rectGrob(gp = gpar(fill = alpha("white", 0.0), col = NA)),  # borderless
  textGrob(
    label = paste(summary_lines, collapse = "\n"),
    x = unit(0, "npc"), y = unit(1, "npc"),
    just = c("left", "top"),
    gp = gpar(fontsize = 9)
  )
)

p_table <- ggplot() +
  theme_void() +
  annotation_custom(table_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

# ----------------------------
# Plot (main)
# ----------------------------
p_main <- ggraph(g_final2, layout = "manual", x = coords[,1], y = coords[,2]) +

  # Edges (width legend uses raw widths)
  geom_edge_link(
    aes(edge_width = width_raw),
    colour = edge_col,
    alpha = edge_alpha,
    arrow = arrow(length = unit(arrow_mm, "mm"), type = "closed"),
    end_cap = circle(2.4, "mm"),
    show.legend = TRUE
  ) +
  scale_edge_width(
    range = edge_width_range,
    breaks = function(x) pretty(x, n = 5),
    name = "Edge weight"
  ) +

  # non-DE matched genes (grey)
  geom_node_point(
    data = function(x) dplyr::filter(x, node_type == "Gene", has_deseq, !is_DE),
    aes(size = size, shape = shape),
    fill = nonDE_fill,
    colour = node_border,
    stroke = 0.6
  ) +

  # unmatched genes (grey85)
  geom_node_point(
    data = function(x) dplyr::filter(x, node_type == "Gene", !has_deseq),
    aes(size = size, shape = shape),
    fill = unmatched_fill,
    colour = node_border,
    stroke = 0.6
  ) +

  # TFs + DE genes (directional log2FC)
  geom_node_point(
    data = function(x) dplyr::filter(x, node_type == "TF" | is_DE),
    aes(fill = log2FoldChange, size = size, shape = shape),
    colour = node_border,
    stroke = 0.6
  ) +
  scale_size_identity() +
  scale_shape_identity() +
  scale_fill_gradient2(
    low = fc_low, mid = fc_mid, high = fc_high,
    midpoint = 0,
    limits = fc_limits,
    breaks = fc_breaks,                # optional
    oob = scales::squish,              # keeps out-of-range values at the ends
    na.value = "grey80",
    name = paste0("log2FC\n(", target_label, " vs ", source_label, ")")
  ) +

  # Labels
  geom_node_label(
    aes(label = plot_label),
    data = function(x) dplyr::filter(x, !is.na(plot_label)),
    nudge_x = label_nudge_x,
    nudge_y = label_nudge_y,
    hjust = label_hjust,
    vjust = label_vjust,
    fill = scales::alpha(label_fill_col, label_fill_alpha),
    colour = label_text_col,
    label.size = label_border_size,
    label.padding = unit(label_box_padding_lines, "lines"),
    label.r = unit(label_box_radius_lines, "lines"),
    show.legend = FALSE
  ) +

  theme_void() +
  ggtitle(paste0(
    "2-step branching subnet around ", tf_of_interest,
    " | plot minEdge≥", min_edge_weight,
    " | DE: padj≤", padj_cutoff, ", absLFC≥", abs_lfc_cutoff,
    " | outwardEdges=", keep_outward_edges,
    " | edgeWidth=", edge_width_mode
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  )

# Stack the table under the legends by placing it in a right-side column.
# The legends live on p_main (right side). We'll give the table its own right column
# and align heights using patchwork.
final_plot <- p_main | p_table +
  plot_layout(widths = c(1, 0.33))

# ----------------------------
# Save
# ----------------------------
ggsave(out_pdf, final_plot, width = out_width, height = out_height, device = cairo_pdf)
message("Saved plot to: ", out_pdf)
