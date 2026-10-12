# install.packages("universalmotif")  # if needed
library(universalmotif)
library(stringr)
library(dplyr)

# INPUTS
in_jaspar <- "/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/motifs/motifs.jaspar"    # the original file
out_jaspar <- "/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/motifs/jaspar_expanded_motif.jaspar"  # feed this to tfactivity via --motifs
map_csv    <- "/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/motifs/jaspar_expanded_map.csv" # sidecar metadata for later joins

# --- Helpers ---
ensure_extrainfo_char <- function(m) {
  # extrainfo must be a named character vector (or length-0 character())
  if (is.null(m@extrainfo)) m@extrainfo <- character()
  if (!is.character(m@extrainfo)) m@extrainfo <- as.character(m@extrainfo)
  m
}
safe_label <- function(m) {
  # prefer TF label from altname; fallback to name if needed
  lbl <- if (!is.null(m@altname) && nzchar(m@altname)) m@altname else m@name
  str_trim(ifelse(is.null(lbl), "", lbl))
}
append_extrainfo <- function(m, key, value) {
  m <- ensure_extrainfo_char(m)
  m@extrainfo[key] <- value
  m
}

# --- Read original motifs (list of universalmotif) ---
motifs <- read_jaspar(in_jaspar)

# --- Expand dimers A::B into two copies (altname=A and altname=B) ---
expand_dimer <- function(m) {
  label <- safe_label(m)    # TF label (may be dimer)
  id    <- m@name           # JASPAR ID (e.g., MA0006.1)

  if (str_detect(label, "[:]{2}")) {
    parts <- str_split(label, "[:]{2}", n = 2, simplify = TRUE)
    A <- str_trim(parts[1]); B <- str_trim(parts[2])
    dimer_lab <- paste0(A, "::", B)

    # Copy 1 (monomer A): keep PWM, new unique ID, altname = A
    mA <- m
    mA@name    <- paste0(id, ".A")  # ensure unique ID
    mA@altname <- A                 # used by tfactivity for expression filtering
    mA <- append_extrainfo(mA, "dimer", dimer_lab)
    mA <- append_extrainfo(mA, "monomer", A)
    mA <- append_extrainfo(mA, "expanded_from", id)

    # Copy 2 (monomer B)
    mB <- m
    mB@name    <- paste0(id, ".B")
    mB@altname <- B
    mB <- append_extrainfo(mB, "dimer", dimer_lab)
    mB <- append_extrainfo(mB, "monomer", B)
    mB <- append_extrainfo(mB, "expanded_from", id)

    list(mA, mB)
  } else {
    # Monomer: ensure altname is usable for filtering
    if (is.null(m@altname) || !nzchar(m@altname)) m@altname <- label
    m <- append_extrainfo(m, "dimer", NA_character_)
    m <- append_extrainfo(m, "monomer", m@altname)
    m <- append_extrainfo(m, "expanded_from", id)
    list(m)
  }
}

expanded <- unlist(lapply(motifs, expand_dimer), recursive = FALSE)

# --- Write expanded motifs as JASPAR (ID = @name, TF label = @altname) ---
write_jaspar(expanded, file = out_jaspar, overwrite = TRUE)

# --- Sidecar mapping for later joins (e.g., to display dimer labels) ---
map_rows <- lapply(expanded, function(m) {
  m <- ensure_extrainfo_char(m)
  data.frame(
    id      = m@name,                        # JASPAR ID written in file
    tf      = m@altname,                     # monomer used for filtering
    dimer   = ifelse("dimer"   %in% names(m@extrainfo),   m@extrainfo[["dimer"]],   NA),
    monomer = ifelse("monomer" %in% names(m@extrainfo),   m@extrainfo[["monomer"]], NA),
    from_id = ifelse("expanded_from" %in% names(m@extrainfo), m@extrainfo[["expanded_from"]], NA),
    width   = ncol(m@motif),
    stringsAsFactors = FALSE
  )
})
map <- do.call(rbind, map_rows)
write.csv(map, map_csv, row.names = FALSE)

cat(sprintf(
  "Expanded %d → %d motifs. Wrote: %s\nMapping: %s\n",
  length(motifs), length(expanded), out_jaspar, map_csv
))