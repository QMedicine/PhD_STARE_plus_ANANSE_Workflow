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

# - Uses FILTERED diffnetwork for PLOT (min_edge_weight)
# Selected final force-directed notebook section; adapters are recorded in provenance.
# ----------------------------
# Input files
# ----------------------------
influence_file <- "influence.tsv"

diffnet_file <- "diffnetwork.tsv"

deseq_file <- "directional_de.tsv"


# ----------------------------
# User knobs
# ----------------------------
tf_of_interest <- "FOCAL_TF"
max_steps <- 2

# Plot edges only (table uses unfiltered diffnetwork)
min_edge_weight_plot <- 0.01

# DE definition
padj_cutoff    <- 0.05
abs_lfc_cutoff <- 0.25

# Densify plot (optional)
add_top_nonDE_genes <- FALSE
n_top_nonDE_genes   <- 40

include_unmatched_genes <- FALSE
n_unmatched_genes       <- 25

# Branching structure
keep_outward_edges <- TRUE

# Edge width definition
edge_width_mode <- "path2"  # "edge" or "path2"

# Plot layout (THIS is what fixes the 3-row issue)
plot_layout <- "nicely"
plot_seed <- 1              # reproducible force layouts
tf_y_offset_frac <- 0.20    # move TF a bit upward after centering (0 = centered)

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

# Fixed colour scale across plots
fc_limits <- c(-2, 2)     # pick your global min/max (e.g. c(-3, 3))
fc_breaks <- c(-2, -1, 0, 1, 2)   # optional: consistent legend ticks

# Labels
label_tfs <- TRUE
label_top_de_genes <- TRUE
n_label_de_genes <- 25

label_nudge_x <- 0.00
label_nudge_y <- 0.05
label_hjust   <- 0.5
label_vjust   <- -0.25
label_box_padding_lines <- 0.12
label_box_radius_lines  <- 0.12
label_fill_col   <- "white"
label_fill_alpha <- 0.90
label_text_col   <- "black"
label_border_size <- 0.0

# Right-hand table style
table_fontsize <- 10
table_bg_alpha <- 0.0

# Output
out_dir <- "plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_pdf <- file.path(out_dir, paste0("ANANSE_2step_branch_", tf_of_interest, "_", source_label, "_to_", target_label, ".pdf"))
out_width <- 15
out_height <- 6


# ----------------------------
# Helpers
# ----------------------------
stop_if_missing <- function(path) if (!file.exists(path)) stop("File not found: ", path)

msg_counts <- function(label, g) {
  message(label, ": nodes=", vcount(g), " edges=", ecount(g))
}

# Convert ALL-CAPS symbols to mouse-style Titlecase for DISPLAY ONLY
mouse_case <- function(x) {
  x2 <- as.character(x)
  is_all_caps <- (x2 == toupper(x2)) & grepl("^[A-Z0-9_]+$", x2) & nchar(x2) > 1
  y <- x2
  y[is_all_caps] <- {
    z <- tolower(x2[is_all_caps])
    paste0(toupper(substr(z, 1, 1)), substr(z, 2, nchar(z)))
  }
  y
}

extract_legend_grob <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

wrap_grob_plot <- function(grob) {
  ggplot() + theme_void() + annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}

make_table_grob <- function(lines, fontsize = 9, bg_alpha = 0.0) {
  grobTree(
    rectGrob(gp = gpar(fill = alpha("white", bg_alpha), col = NA)),
    textGrob(
      label = paste(lines, collapse = "\n"),
      x = unit(0, "npc"), y = unit(1, "npc"),
      just = c("left", "top"),
      gp = gpar(fontsize = fontsize)
    )
  )
}

get_plot_coords <- function(g, layout_name, seed = 1) {
  set.seed(seed)
  layout_name <- tolower(layout_name)
  if (layout_name == "nicely") return(layout_nicely(g))
  if (layout_name == "fr")     return(layout_with_fr(g))
  if (layout_name == "kk")     return(layout_with_kk(g))
  if (layout_name == "circle") return(layout_in_circle(g))
  if (layout_name == "random") return(layout_randomly(g))
  # fallback
  layout_nicely(g)
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

diffnet_plot <- diffnet_full %>% filter(weight >= min_edge_weight_plot)
if (nrow(diffnet_plot) == 0) stop("No plot edges after min_edge_weight_plot filtering; lower it.")

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

tf_inf <- inf %>% filter(factor == tf_of_interest) %>% slice(1)
if (nrow(tf_inf) == 0) stop("TF not found in influence table: ", tf_of_interest)

# ============================================================
# TABLE metrics from TRUE ANANSE diffnetwork (UNFILTERED)
# ============================================================
g_all_table <- graph_from_data_frame(diffnet_full, directed = TRUE)
v_tbl <- match(tf_of_interest, V(g_all_table)$name)
if (is.na(v_tbl)) stop("TF not found in unfiltered diffnetwork graph: ", tf_of_interest)

ego_tbl_vids <- ego(g_all_table, order = max_steps, nodes = v_tbl, mode = "out")[[1]]
g_tbl0 <- induced_subgraph(g_all_table, vids = ego_tbl_vids)

# IMPORTANT: hop distances must be UNWEIGHTED
dmat_tbl <- distances(
  g_tbl0,
  v = match(tf_of_interest, V(g_tbl0)$name),
  to = V(g_tbl0),
  mode = "out",
  weights = NA
)
dist_tbl <- as.numeric(dmat_tbl[1, ])

tbl_direct_targets_total    <- sum(dist_tbl == 1, na.rm = TRUE)
tbl_secondary_targets_total <- sum(dist_tbl == 2, na.rm = TRUE)
tbl_targets_within2_total   <- sum(dist_tbl >= 1 & dist_tbl <= max_steps, na.rm = TRUE)

reachable_tbl_nodes <- V(g_tbl0)$name[dist_tbl >= 1 & dist_tbl <= max_steps]

tbl_de_strict_nodes <- de2 %>% filter(gene_id %in% reachable_tbl_nodes, is_DE_strict) %>% pull(gene_id)
tbl_de_padj_nodes   <- de2 %>% filter(gene_id %in% reachable_tbl_nodes, is_DE_padj)   %>% pull(gene_id)

tbl_de_strict <- length(tbl_de_strict_nodes)
tbl_de_padj   <- length(tbl_de_padj_nodes)

tbl_de_strict_genes_only <- sum(!(tbl_de_strict_nodes %in% tf_universe))
tbl_de_padj_genes_only   <- sum(!(tbl_de_padj_nodes %in% tf_universe))

message("TABLE sanity: out-degree=", degree(g_all_table, v = v_tbl, mode = "out"),
        " | ego nodes=", vcount(g_tbl0), " | ego edges=", ecount(g_tbl0),
        " | direct=", tbl_direct_targets_total, " | secondary=", tbl_secondary_targets_total)

# ============================================================
# PLOT graph from FILTERED diffnetwork
# ============================================================
g_all <- graph_from_data_frame(diffnet_plot, directed = TRUE)
msg_counts("g_all (plot graph, filtered)", g_all)

v_plot <- match(tf_of_interest, V(g_all)$name)
if (is.na(v_plot)) stop("TF not found in filtered plot graph. Lower min_edge_weight_plot.")

ego_vids <- ego(g_all, order = max_steps, nodes = v_plot, mode = "out")[[1]]
g_sub0 <- induced_subgraph(g_all, vids = ego_vids)
msg_counts("g_sub0 (plot ego induced)", g_sub0)
if (ecount(g_sub0) == 0) stop("g_sub0 has 0 edges; lower min_edge_weight_plot.")

# Hop layers (UNWEIGHTED)
dmat <- distances(
  g_sub0,
  v = match(tf_of_interest, V(g_sub0)$name),
  to = V(g_sub0),
  mode = "out",
  weights = NA
)
layer_vec <- as.integer(dmat[1, ])
names(layer_vec) <- V(g_sub0)$name
V(g_sub0)$layer <- layer_vec[V(g_sub0)$name]

keep_nodes <- V(g_sub0)$name[is.finite(V(g_sub0)$layer) & V(g_sub0)$layer <= max_steps]
g_sub <- induced_subgraph(g_sub0, vids = keep_nodes)
V(g_sub)$layer <- layer_vec[V(g_sub)$name]

# Optional outward-edge filter (keeps branching, reduces back-edges)
if (isTRUE(keep_outward_edges)) {
  edf <- as_data_frame(g_sub, what = "edges") %>%
    rename(source = from, target = to) %>%
    mutate(layer_u = layer_vec[source],
           layer_v = layer_vec[target])

  edf2 <- edf %>%
    filter(is.finite(layer_u), is.finite(layer_v),
           layer_u < layer_v,
           layer_v <= max_steps) %>%
    select(source, target, weight)

  if (nrow(edf2) == 0) stop("Outward-edge filter removed all plot edges. Disable keep_outward_edges or lower min_edge_weight_plot.")

  g_sub <- graph_from_data_frame(
    edf2,
    directed = TRUE,
    vertices = data.frame(name = unique(c(edf2$source, edf2$target)))
  )
  V(g_sub)$layer <- layer_vec[V(g_sub)$name]
  E(g_sub)$weight <- edf2$weight
  msg_counts("g_sub (after outward-edge filter)", g_sub)
}

# Node annotations (plot graph)
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
    log2FoldChange = log2FoldChange,
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
if (ecount(g_final) == 0) stop("g_final has 0 plot edges; relax plot filters.")

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

# Labels: store DISPLAY labels using mouse_case()
plot_label <- rep(NA_character_, vcount(g_final))
names(plot_label) <- V(g_final)$name

if (isTRUE(label_tfs)) {
  tf_names <- V(g_final)$name[V(g_final)$node_type == "TF"]
  plot_label[tf_names] <- mouse_case(tf_names)
}
if (isTRUE(label_top_de_genes)) {
  de_gene_names <- V(g_final)$name[V(g_final)$node_type == "Gene" & V(g_final)$is_DE %in% TRUE]
  abs_lfc_vec <- V(g_final)$abs_lfc
  abs_lfc_vec[is.na(abs_lfc_vec)] <- 0
  ord <- order(abs_lfc_vec[match(de_gene_names, V(g_final)$name)], decreasing = TRUE)
  top_genes <- head(de_gene_names[ord], n_label_de_genes)
  plot_label[top_genes] <- mouse_case(top_genes)
}
V(g_final)$plot_label <- unname(plot_label[V(g_final)$name])

# Edge widths (legend uses width_raw)
edf_final <- as_data_frame(g_final, what = "edges") %>%
  rename(source = from, target = to) %>%
  mutate(abs_weight = abs(weight))

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

# Rebuild plot graph with width_raw + vertex attrs
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

# Drop isolated nodes before layout
iso <- V(g_final2)[degree(g_final2, mode = "all") == 0]
if (length(iso) > 0) g_final2 <- delete_vertices(g_final2, iso)

# ---- Force-directed plot layout (fixes the "3 rows" problem) ----
coords <- get_plot_coords(g_final2, plot_layout, seed = plot_seed)

# center TF, then move it slightly upward if requested
center_idx <- match(tf_of_interest, V(g_final2)$name)
coords <- sweep(coords, 2, coords[center_idx, ], FUN = "-")
span_y <- diff(range(coords[,2]))
coords[,2] <- coords[,2] + span_y * tf_y_offset_frac

# ----------------------------
# Build plot with legends (for extraction)
# ----------------------------
p_with_legends <- ggraph(g_final2, layout = "manual", x = coords[,1], y = coords[,2]) +
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
  geom_node_point(
    data = function(x) dplyr::filter(x, node_type == "Gene", has_deseq, !is_DE),
    aes(size = size, shape = shape),
    fill = nonDE_fill,
    colour = node_border,
    stroke = 0.6
  ) +
  geom_node_point(
    data = function(x) dplyr::filter(x, node_type == "Gene", !has_deseq),
    aes(size = size, shape = shape),
    fill = unmatched_fill,
    colour = node_border,
    stroke = 0.6
  ) +
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
  ggtitle(stringr::str_wrap(paste0(
    "2-step branching subnet around ", mouse_case(tf_of_interest),
    " | plot minEdge≥", min_edge_weight_plot,
    " | DE: padj≤", padj_cutoff, ", absLFC≥", abs_lfc_cutoff,
    " | outwardEdges=", keep_outward_edges,
    " | edgeWidth=", edge_width_mode,
    " | layout=", plot_layout
  ), width = 90)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  )

legend_grob <- extract_legend_grob(p_with_legends)
if (is.null(legend_grob)) stop("Could not extract legend grob; unexpected.")

p_main <- p_with_legends + theme(legend.position = "none")

# ----------------------------
# Right-side table (scaled values requested)
# ----------------------------
summary_lines <- c(
  paste0("TF: ", mouse_case(tf_of_interest)),
  paste0("Direct targets (1-step, total): ", tbl_direct_targets_total),
  paste0("Secondary targets (2-step, total): ", tbl_secondary_targets_total),
  paste0("Targets within ", max_steps, " steps (total): ", tbl_targets_within2_total),
  paste0("DE reachable (padj & absLFC≥", abs_lfc_cutoff, "): ", tbl_de_strict,
         "  (genes-only: ", tbl_de_strict_genes_only, ")"),
  paste0("DE reachable (padj only): ", tbl_de_padj,
         "  (genes-only: ", tbl_de_padj_genes_only, ")"),
  "",
  paste0("target_score_scaled: ", signif(tf_inf$target_score_scaled, 4)),
  paste0("G_score_scaled: ", signif(tf_inf$G_score_scaled, 4)),
  paste0("influence_score (scaled): ", signif(tf_inf$influence_score, 4))
)

table_grob <- make_table_grob(summary_lines, fontsize = table_fontsize, bg_alpha = table_bg_alpha)

p_legend <- wrap_grob_plot(legend_grob)
p_table  <- wrap_grob_plot(table_grob)

right_col <- p_legend / p_table + plot_layout(heights = c(0.62, 0.38))
final_plot <- p_main | right_col + plot_layout(widths = c(1, 0.36))

# ----------------------------
# Save
# ----------------------------
ggsave(out_pdf, final_plot, width = out_width, height = out_height, device = cairo_pdf, bg = "white")
message("Saved plot to: ", out_pdf)
