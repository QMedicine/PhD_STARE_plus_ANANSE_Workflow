set.seed(1) # Reproducible layout
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
})

# ----------------------------
# User-configurable parameters
# ----------------------------

# Input files (your paths)
influence_file <- "influence.tsv"
diffnet_file <- "diffnetwork.tsv"
deseq_file <- "directional_de.tsv"

# Core selection knobs
n_tfs             <- 10        # top TFs from influence file
rank_col          <- "influence_score" # or other column
min_abs_edge_wt   <- 0.01     # edge included if abs(weight) >= this
target_mode       <- "per_tf" # "per_tf" or "overall"
n_targets         <- 15       # if per_tf: per TF; if overall: total across TFs
# --- TF-target handling ---
tf_target_fill <- "plum3"      # light purple for TFs that are targets but not in top_tfs
keep_tf_targets <- TRUE        # include TF targets beyond "top N genes"
n_tf_targets_per_tf <- 15      # cap TF-target edges per source TF (prevents hairballs)
always_keep_top_tf_to_top_tf_edges <- TRUE

# DE target highlighting thresholds
padj_cutoff       <- 0.05
lfc_cutoff        <- 0.0      # set to e.g. 0.25 if you want a minimum effect size

# Layout options
layout_name <- "nicely"
center_top_tf     <- FALSE     # if TRUE, tries to place highest-ranked TF centrally (best with layout_name="star")
circle_top_tf_12oclock <- FALSE # for circle layout: place highest-ranked TF at the top position

# Aesthetics
tf_fill           <- "purple"
gene_fill_default <- "grey30"
up_fill           <- "red3"
down_fill         <- "blue3"
node_border       <- "grey10"

tf_size_range     <- c(10, 20)   # min/max size for TF nodes (scaled by influence)
gene_size         <- 10

edge_col <- "black"
edge_alpha        <- 1
edge_width_range  <- c(0.2, 1.6) # scaled by abs(weight)
arrow_mm          <- 2.5

# Labels
label_tfs         <- TRUE
label_genes       <- TRUE
max_gene_labels   <- 50          # avoid unreadable hairballs
gene_label_by     <- "abs_lfc"    # "abs_lfc" or "abs_edge"

# Label frame (background box behind the text)
label_fill_col   <- "white"   # frame fill colour
label_fill_alpha <- 0.30      # 0 (transparent) -> 1 (opaque)

label_text_col   <- "black"   # text colour

label_border_col <- "white"   # frame border colour (see note below)
label_border_size <- 0.00     # 0 removes border

label_padding_lines <- 0.15   # how much padding around the text
label_radius_lines  <- 0.15   # rounded corners

# Output
out_pdf <- "top10.pdf"
out_width         <- 10
out_height        <- 8

# ----------------------------
# Helper functions
# ----------------------------

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

read_influence <- function(path) {
  # Reads ANANSE influence output. Keeps all cols but requires at least: factor + rank_col
  read_tsv(path, na = c("", "NA"), show_col_types = FALSE)
}

read_diffnet <- function(path) {
  # Reads differential network edges (source -> target with weight)
  read_tsv(path, na = c("", "NA"), show_col_types = FALSE) %>%
    select(source, target, weight)
}

read_deseq <- function(path) {
  # Reads DESeq2 results (gene_id is gene symbol in your case)
  read_tsv(path, na = c("", "NA"), show_col_types = FALSE) %>%
    select(gene_id, log2FoldChange, padj)
}

make_circle_layout <- function(g, top_tf = NULL, top_at_12 = TRUE) {
  # Manual circle coordinates so we can put a specific node at 12 o'clock.
  vnames <- V(g)$name
  n <- length(vnames)

  # Default order: keep TFs first, then genes, so circle is a bit more interpretable.
  is_tf <- V(g)$type == "TF"
  ord <- c(which(is_tf), which(!is_tf))

  # If requested, rotate so top_tf is first in ord, which we map to 12 o'clock.
  if (!is.null(top_tf) && top_at_12 && top_tf %in% vnames) {
    idx <- which(vnames[ord] == top_tf)
    if (length(idx) == 1) ord <- c(ord[idx:length(ord)], ord[1:(idx - 1)])
  }

  angles <- seq(from = pi/2, length.out = n, by = -2*pi/n) # start at 12 o'clock
  x <- cos(angles)
  y <- sin(angles)

  coords <- matrix(NA_real_, nrow = n, ncol = 2)
  coords[ord, 1] <- x
  coords[ord, 2] <- y
  coords
}

pick_layout <- function(g, layout_name, top_tf = NULL) {
  # Returns either:
  # - list(layout="manual", x=..., y=...) for manual layouts, or
  # - list(layout="<ggraph layout name>") for built-in layouts.
  layout_name <- tolower(layout_name)

  if (layout_name == "circle") {
    coords <- make_circle_layout(g, top_tf = top_tf, top_at_12 = circle_top_tf_12oclock)
    return(list(layout = "manual", x = coords[,1], y = coords[,2]))
  }

  if (layout_name == "star") {
    # Radial star layout centred at top_tf (best interpretation of "place top TF in centre").
    if (!is.null(top_tf) && top_tf %in% V(g)$name) {
      center_idx <- which(V(g)$name == top_tf)
      coords <- layout_as_star(g, center = center_idx, mode = "out")
      return(list(layout = "manual", x = coords[,1], y = coords[,2]))
    } else {
      coords <- layout_as_star(g, mode = "out")
      return(list(layout = "manual", x = coords[,1], y = coords[,2]))
    }
  }

  # Map some friendly names to ggraph layouts
  if (layout_name %in% c("fr", "force", "force-directed")) return(list(layout = "fr"))
  if (layout_name %in% c("kk", "kamada-kawai")) return(list(layout = "kk"))
  if (layout_name %in% c("random", "randomly")) return(list(layout = "randomly"))
  if (layout_name %in% c("nicely", "nice")) return(list(layout = "nicely"))
  if (layout_name %in% c("hierarchical", "sugiyama")) return(list(layout = "sugiyama"))

  # Fallback
  list(layout = "fr")
}

# ----------------------------
# Main: build and plot network
# ----------------------------

stop_if_missing(influence_file)
stop_if_missing(diffnet_file)
stop_if_missing(deseq_file)

inf <- read_influence(influence_file)
diffnet <- read_diffnet(diffnet_file)
de <- read_deseq(deseq_file)

# Basic input validation
required_inf_cols <- c("factor", rank_col)
missing_inf_cols <- setdiff(required_inf_cols, colnames(inf))
if (length(missing_inf_cols) > 0) {
  stop("Influence file is missing required columns: ", paste(missing_inf_cols, collapse = ", "))
}

# Any factor appearing in ANANSE influence output is treated as a TF name universe
tf_universe <- unique(inf$factor)

# Select top TFs
top_tfs_tbl <- inf %>%
  arrange(desc(.data[[rank_col]])) %>%
  slice_head(n = n_tfs) %>%
  mutate(rank = row_number())

top_tfs <- top_tfs_tbl$factor
top_tf1 <- top_tfs[1]

message("Top TFs selected (", n_tfs, "): ", paste(top_tfs, collapse = ", "))

# Filter edges to those outgoing from top TFs
edges0 <- diffnet %>%
  filter(source %in% top_tfs)

message("Edges from top TFs (pre-weight filter): ", nrow(edges0))

# Filter by abs(weight)
edges1 <- edges0 %>%
  mutate(abs_weight = abs(weight)) %>%
  filter(abs_weight >= min_abs_edge_wt)

message("Edges after abs(weight) >= ", min_abs_edge_wt, ": ", nrow(edges1))

if (nrow(edges1) == 0) {
  stop("No edges remain after filtering. Lower min_abs_edge_wt or increase n_tfs.")
}

# Split edges into:
#  (A) TF targets (target is a TF name from tf_universe)
#  (B) non-TF gene targets
edges1_tf_targets   <- edges1 %>% filter(target %in% tf_universe)
edges1_gene_targets <- edges1 %>% filter(!(target %in% tf_universe))

# A) Always keep edges between top TFs (prevents JUND<->FOS being dropped by "top N targets")
edges_keep_top_tf_pairs <- edges1 %>%
  filter(source %in% top_tfs, target %in% top_tfs)

# B) Select gene targets as before (top N per TF or overall), but only for NON-TF targets
if (tolower(target_mode) == "per_tf") {
  edges_gene_sel <- edges1_gene_targets %>%
    group_by(source) %>%
    arrange(desc(abs_weight)) %>%
    slice_head(n = n_targets) %>%
    ungroup()
} else if (tolower(target_mode) == "overall") {
  edges_gene_sel <- edges1_gene_targets %>%
    arrange(desc(abs_weight)) %>%
    slice_head(n = n_targets)
} else {
  stop("target_mode must be 'per_tf' or 'overall'")
}

# C) Optionally include TF targets (capped per TF so it doesn't explode)
edges_tf_sel <- edges_keep_top_tf_pairs[0, ]  # empty default with same cols
if (isTRUE(keep_tf_targets)) {
  edges_tf_sel <- edges1_tf_targets %>%
    group_by(source) %>%
    arrange(desc(abs_weight)) %>%
    slice_head(n = n_tf_targets_per_tf) %>%
    ungroup()
}

# Combine: gene targets + TF targets + (optionally) guaranteed topTF->topTF edges
edges2 <- bind_rows(
  edges_gene_sel,
  edges_tf_sel,
  if (isTRUE(always_keep_top_tf_to_top_tf_edges)) edges_keep_top_tf_pairs else NULL
) %>%
  distinct(source, target, .keep_all = TRUE)

message("Edges kept total after TF-aware selection: ", nrow(edges2))

# Nodes: TFs + targets
targets <- unique(edges2$target)
nodes <- tibble(name = unique(c(top_tfs, unique(edges2$target)))) %>%
  mutate(
    # Any node in tf_universe is considered a TF (even if it's only a target)
    is_tf = name %in% tf_universe,
    is_top_tf = name %in% top_tfs,
    type = if_else(is_tf, "TF", "Gene")
  )

# Join TF influence for sizing (only for TF nodes)
nodes <- nodes %>%
  left_join(top_tfs_tbl %>% select(factor, !!rank_col, rank), by = c("name" = "factor")) %>%
  rename(rank_value = !!rank_col)

# Join DE info for gene targets (match by gene symbol in gene_id; Ensembl-like targets won't match and will stay NA)
de2 <- de %>%
  mutate(
    de_status = case_when(
      is.na(padj) | is.na(log2FoldChange) ~ NA_character_,
      padj <= padj_cutoff & log2FoldChange >=  lfc_cutoff ~ "up",
      padj <= padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "down",
      TRUE ~ "ns"
    ),
    abs_lfc = abs(log2FoldChange)
  )

nodes <- nodes %>%
  left_join(de2, by = c("name" = "gene_id"))

nodes <- nodes %>%
  mutate(
    fill = case_when(
      is_top_tf ~ tf_fill,            # strong purple
      is_tf     ~ tf_target_fill,     # light purple for TF targets
      de_status == "up"   ~ up_fill,
      de_status == "down" ~ down_fill,
      TRUE ~ gene_fill_default
    )
  )

# Sanity-check: how many targets matched DESeq2?
n_targets_total <- length(targets)
n_targets_matched <- nodes %>%
  filter(type == "Gene") %>%
  summarise(n = sum(!is.na(log2FoldChange))) %>%
  pull(n)

message("Targets total: ", n_targets_total,
        " | Targets matched in DESeq2: ", n_targets_matched,
        " | Unmatched targets (kept, labelled as-is): ", n_targets_total - n_targets_matched)

# Build igraph object
g <- graph_from_data_frame(
  d = edges2 %>% select(source, target, weight, abs_weight),
  vertices = nodes %>% select(name, type, rank_value, rank, log2FoldChange, padj, de_status, abs_lfc, fill),
  directed = TRUE
)

# Compute node sizes
# TF size scales with rank_value (influence); genes fixed.
tf_vals <- V(g)$rank_value
tf_vals[is.na(tf_vals)] <- 0

# Rescale only among TFs to avoid genes contaminating range
tf_only <- V(g)$type == "TF"
tf_scaled <- rep(NA_real_, vcount(g))
if (sum(tf_only) >= 2) {
  tf_scaled[tf_only] <- rescale(tf_vals[tf_only], to = tf_size_range)
} else {
  tf_scaled[tf_only] <- mean(tf_size_range)
}
node_size <- ifelse(V(g)$type == "TF", tf_scaled, gene_size)

V(g)$node_size <- node_size

# Compute edge width scaling
E(g)$edge_width <- rescale(E(g)$abs_weight, to = edge_width_range)

E(g)$width_raw <- E(g)$abs_weight # Audit repair: restore the raw edge attribute referenced by the original scale

# Layout selection
# If user requests centering of top TF, "star" is the most reliable interpretation.
if (center_top_tf && layout_name != "star") {
  message("center_top_tf=TRUE: switching layout to 'star' to reliably centre top TF.")
  layout_name_use <- "star"
} else {
  layout_name_use <- layout_name
}

lay <- pick_layout(g, layout_name_use, top_tf = top_tf1)

# Decide labels
# Always label TFs (optional), and label up to max_gene_labels genes by chosen criterion.
label_tbl <- tibble(name = V(g)$name, type = V(g)$type) %>%
  mutate(label = NA_character_)

if (label_tfs) {
  label_tbl <- label_tbl %>%
    mutate(label = if_else(type == "TF", name, label))
}

if (label_genes) {
  gene_names <- V(g)$name[V(g)$type == "Gene"]
  gene_df <- nodes %>%
    filter(type == "Gene", name %in% gene_names) %>%
    mutate(score_for_label = case_when(
      gene_label_by == "abs_lfc" ~ if_else(is.na(abs_lfc), -Inf, abs_lfc),
      TRUE ~ -Inf
    ))

  if (gene_label_by == "abs_edge") {
    # label genes with strongest incoming abs edge weight among the selected edges
    incoming_max <- edges2 %>%
      group_by(target) %>%
      summarise(max_abs_edge = max(abs_weight), .groups = "drop")
    gene_df <- gene_df %>%
      left_join(incoming_max, by = c("name" = "target")) %>%
      mutate(score_for_label = if_else(is.na(max_abs_edge), -Inf, max_abs_edge))
  }

  top_gene_labels <- gene_df %>%
    arrange(desc(score_for_label)) %>%
    slice_head(n = max_gene_labels) %>%
    pull(name)

  label_tbl <- label_tbl %>%
    mutate(label = if_else(type == "Gene" & name %in% top_gene_labels, name, label))
}

# Map labels onto graph vertices so ggraph's node data includes them
label_map <- setNames(label_tbl$label, label_tbl$name)
V(g)$plot_label <- unname(label_map[V(g)$name])

write_tsv(nodes, file.path(dirname(out_pdf), "top10_nodes.tsv"))
write_tsv(edges2, file.path(dirname(out_pdf), "top10_edges.tsv"))

# ----------------
# Plot
# ----------------

# 1) Build the base ggraph object WITHOUT passing x/y unless manual
if (identical(lay$layout, "manual")) {
  p <- ggraph(g, layout = "manual", x = lay$x, y = lay$y)
} else {
  p <- ggraph(g, layout = lay$layout)
}

# 2) Add layers (same as before)
p <- p +
  geom_edge_link(
    aes(edge_width = width_raw),   # IMPORTANT: use the *raw* values, not the rescaled ones
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
guides(
    edge_width = guide_legend(
      override.aes = list(
       # ggraph draws edges as lines; legend will show line widths
        alpha = 1
     )
    )
  ) +
  geom_node_point(aes(size = node_size, fill = log2FoldChange, shape = type),
                  colour = node_border, stroke = 0.6) +
  scale_shape_manual(values = c(TF = 22, Gene = 21), name = "Node type") +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "red3", midpoint = 0,
                      limits = fc_limits, breaks = fc_breaks, oob = scales::squish,
                      na.value = "grey80", name = paste0("log2FC\n(", target_label, " vs ", source_label, ")")) +
  scale_size_identity() +
  geom_node_label(
    aes(label = plot_label),
    repel = TRUE, max.overlaps = Inf, nudge_y = 0.15,
    data = function(x) dplyr::filter(x, !is.na(plot_label)),  # keeps x/y from  ggraph layout
    fill = scales::alpha(label_fill_col, label_fill_alpha),
    colour = label_text_col,
    label.size = label_border_size,
    label.padding = unit(label_padding_lines, "lines"),
    label.r = unit(label_radius_lines, "lines"),
    show.legend = FALSE
  ) +
  theme_void() +
  ggtitle(paste0("ANANSE: Top ", n_tfs, " TFs and targets (", layout_name_use, " layout)")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# ----------------
# Save (Cairo)
# ----------------
library(Cairo)

# Option A (recommended): just use ggsave with cairo_pdf
ggsave(out_pdf, p, width = out_width, height = out_height, device = cairo_pdf, bg = "white")
message("Saved plot to: ", out_pdf)

# Option B (if you want explicit CairoPDF control):
# CairoPDF(file = out_pdf, width = out_width, height = out_height,
#          pointsize = 12, bg = "transparent", onefile = TRUE)
# print(p)
# dev.off()

# Print to screen if interactive
if (interactive()) print(p)
