suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(apeglm)

  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(cluster)
  library(grid)
})

# Silence ComplexHeatmap rasterization suggestion
ht_opt$message <- FALSE

# ============================================================
# Paths
# ============================================================
# Default: current working directory (recommended to run from repo root)
root <- Sys.getenv("ROOT_DIR", unset = getwd())
root <- normalizePath(root, mustWork = TRUE)

samples_file <- file.path(root, "samples.tsv")
counts_file  <- file.path(root, "counts", "gene_counts.matrix.tsv")

out_dir <- file.path(root, "de")
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Liftoff mapping (IOPOBE -> LFTS + gene_name symbol)
# Header:
# new_gene_id old_gene_id gene_name overlap_bp frac_old frac_new
gene_map_tsv <- file.path(
  root, "L_ferri_bakta", "output_liftoff_gbk", "pypolca_corrected_to_LFTS.tsv"
)

# ============================================================
# Tunables
# ============================================================
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1
MAX_UNION_GENES <- 2000

# Union heatmap gene labels (right side):
# Goal = not all labels; just enough to cover the full height ("head to toe")
LABEL_TARGET_TOTAL    <- NA_integer_  # NA = auto; set 60-90 if you want a fixed number
LABEL_MAX_TOTAL       <- 80           # hard cap
LABEL_MIN_PER_CLUSTER <- 8            # ensure each cluster slice has some labels
LABEL_FONTSIZE        <- 6

# ============================================================
# Helper functions (robust input + DESeq2)
# ============================================================
clean_names <- function(x) {
  x %>%
    stringr::str_remove("^.*/") %>%
    stringr::str_remove("\\.bam$") %>%
    stringr::str_remove("\\.sorted$")
}

find_level <- function(target, levs) {
  if (target %in% levs) return(target)

  m1 <- levs[tolower(levs) == tolower(target)]
  if (length(m1) == 1) return(m1[1])

  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  m2 <- levs[norm(levs) == norm(target)]
  if (length(m2) == 1) return(m2[1])

  m3 <- levs[grepl(norm(target), norm(levs))]
  if (length(m3) == 1) return(m3[1])

  stop(
    "Could not uniquely match condition '", target, "'.\n",
    "Available levels:\n- ", paste(levs, collapse = "\n- ")
  )
}

make_dds <- function(cts, coldata, ref_level) {
  cd <- coldata
  cd$condition <- relevel(factor(cd$condition), ref = ref_level)

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(cts)),
    colData   = cd,
    design    = ~ condition
  )

  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  dds
}

get_shrunken <- function(dds, numerator_level, ref_level) {
  coef_name <- paste0(
    "condition_", make.names(numerator_level),
    "_vs_", make.names(ref_level)
  )

  if (!(coef_name %in% resultsNames(dds))) {
    stop(
      "Coefficient not found: ", coef_name, "\n\n",
      "Available coefficients:\n",
      paste(resultsNames(dds), collapse = "\n")
    )
  }

  res <- results(dds, name = coef_name)
  res_shr <- lfcShrink(dds, coef = coef_name, res = res, type = "apeglm")

  df <- as.data.frame(res_shr) %>%
    tibble::rownames_to_column("gene_id") %>%
    arrange(is.na(padj), padj)

  list(coef_name = coef_name, res = res, res_shr = res_shr, df = df)
}

write_sig_list <- function(df, out_file, padj_cutoff = 0.05, lfc_cutoff = 1) {
  sig <- df %>%
    filter(!is.na(padj), padj < padj_cutoff,
           !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_cutoff)
  readr::write_tsv(sig, out_file)
  nrow(sig)
}

# ============================================================
# QC plot helpers
# ============================================================
plot_ma_pdf <- function(res_shr, out_pdf, main = NULL) {
  pdf(out_pdf, width = 7, height = 6)
  plotMA(res_shr, main = main, ylim = c(-5, 5))
  abline(h = c(-1, 1), lty = 2)
  dev.off()
}

plot_volcano_pdf <- function(df, out_pdf, title,
                             padj_cutoff = 0.05, lfc_cutoff = 1) {
  dd <- df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      padj_plot = pmax(padj, 1e-300),
      neglog10  = -log10(padj_plot),
      sig = (padj < padj_cutoff) & (abs(log2FoldChange) >= lfc_cutoff)
    )

  p <- ggplot(dd, aes(x = log2FoldChange, y = neglog10)) +
    geom_point(aes(color = sig), alpha = 0.75, size = 1) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = title, x = "log2 fold-change (apeglm)", y = "-log10(padj)") +
    theme_bw()

  ggsave(out_pdf, p, width = 7, height = 6)
}

plot_pca_pdf <- function(vsd, coldata, out_pdf) {
  pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    labs(
      title = "PCA (VST)",
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance")
    ) +
    theme_bw()

  ggsave(out_pdf, p, width = 7, height = 6)
}

# ============================================================
# Gene symbol mapping: IOPOBE -> gene_name (symbol) from mapping TSV
# ============================================================
standardize_names <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "_", x))
}

load_gene_symbol_map <- function(map_tsv) {
  if (!file.exists(map_tsv)) {
    message("[INFO] Gene map TSV not found, symbols disabled: ", map_tsv)
    return(list(symbol = setNames(character(0), character(0)),
                label  = setNames(character(0), character(0)),
                table  = NULL))
  }

  gm <- readr::read_tsv(
    map_tsv,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  colnames(gm) <- standardize_names(colnames(gm))

  # Your header: new_gene_id old_gene_id gene_name ...
  new_col   <- intersect(colnames(gm), c("new_gene_id","new_geneid","gene_id","geneid"))[1]
  gname_col <- intersect(colnames(gm), c("gene_name","genename","symbol","gene"))[1]

  if (is.na(new_col)) {
    stop("Could not find new_gene_id/gene_id column in: ", map_tsv,
         "\nColumns:\n- ", paste(colnames(gm), collapse = "\n- "))
  }
  if (is.na(gname_col)) {
    message("[WARN] gene_name column not found in mapping TSV; labels will be gene_id only.")
  }

  gm2 <- gm %>%
    transmute(
      gene_id     = .data[[new_col]],
      gene_symbol = if (!is.na(gname_col)) .data[[gname_col]] else NA_character_
    ) %>%
    mutate(
      gene_id = na_if(gene_id, "NA"),
      gene_id = na_if(gene_id, ""),
      gene_symbol = na_if(gene_symbol, "NA"),
      gene_symbol = na_if(gene_symbol, "")
    ) %>%
    filter(!is.na(gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    mutate(
      gene_label = if_else(!is.na(gene_symbol),
                           paste0(gene_id, " (", gene_symbol, ")"),
                           gene_id)
    )

  sym_map <- setNames(gm2$gene_symbol, gm2$gene_id)
  lab_map <- setNames(gm2$gene_label,  gm2$gene_id)

  message("[INFO] Loaded gene symbols for ", sum(!is.na(gm2$gene_symbol)), " genes from mapping TSV.")
  list(symbol = sym_map, label = lab_map, table = gm2)
}

gene_maps <- load_gene_symbol_map(gene_map_tsv)
gene_symbol_map <- gene_maps$symbol
gene_label_map  <- gene_maps$label

symbol_for_gene <- function(gids) {
  sym <- gene_symbol_map[gids]
  sym[is.na(sym)] <- NA_character_
  unname(sym)
}
label_for_gene <- function(gids) {
  lab <- gene_label_map[gids]
  lab[is.na(lab)] <- gids[is.na(lab)]
  unname(lab)
}

# ============================================================
# ComplexHeatmap utilities
# ============================================================
.height_from_rows <- function(n) max(7, min(42, 0.012 * n + 7))

zscore_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

make_condition_colors <- function(cond_factor) {
  levs <- levels(cond_factor)
  pal <- RColorBrewer::brewer.pal(max(3, length(levs)), "Set2")[seq_along(levs)]
  setNames(pal, levs)
}

make_col_annotation <- function(meta_df) {
  meta_df$condition <- factor(meta_df$condition)
  cond_cols <- make_condition_colors(meta_df$condition)

  HeatmapAnnotation(
    df = data.frame(Condition = meta_df$condition),
    which = "col",
    col = list(Condition = cond_cols),
    na_col = "white",
    annotation_height = unit(4, "mm"),
    gap = unit(1, "mm")
  )
}

choose_pam_k <- function(n_rows) {
  if (n_rows < 10) return(NA_integer_)
  if (n_rows < 30) return(2L)
  if (n_rows < 80) return(3L)
  return(4L)
}

# SAFE correlation distance between ROWS of the input matrix
safe_cor_dist <- function(mat) {
  cmat <- suppressWarnings(cor(t(mat), method = "pearson", use = "pairwise.complete.obs"))
  cmat <- pmax(pmin(cmat, 1), -1)
  cmat[is.na(cmat)] <- 0
  diag(cmat) <- 1
  as.dist(1 - cmat)
}

# Label selection that "fills head to toe" (even spacing by cluster slice)
select_genes_to_label_fill_height <- function(heat_z,
                                             cluster_vec,
                                             row_split_levels,
                                             target_total = NA_integer_,
                                             max_total = 80,
                                             min_per_cluster = 8) {
  stopifnot(all(rownames(heat_z) %in% names(cluster_vec)))

  n_all <- nrow(heat_z)
  if (is.na(target_total) || target_total <= 0) {
    target_total <- max(24, round(n_all / 8))
  }
  target_total <- min(target_total, max_total, n_all)

  sizes <- sapply(row_split_levels, function(cl) sum(cluster_vec[rownames(heat_z)] == cl))
  sizes <- pmax(sizes, 0)
  names(sizes) <- row_split_levels

  k <- round(target_total * sizes / sum(sizes))
  k <- pmax(k, min_per_cluster)
  k <- pmin(k, sizes)

  # adjust to hit target_total
  cur <- sum(k)
  if (cur > target_total) {
    extra <- cur - target_total
    ord <- order(k, decreasing = TRUE)
    i <- 1
    while (extra > 0) {
      cl <- names(k)[ord[i]]
      if (k[cl] > min_per_cluster) {
        k[cl] <- k[cl] - 1
        extra <- extra - 1
      }
      i <- i + 1
      if (i > length(ord)) i <- 1
      if (all(k <= min_per_cluster)) break
    }
  } else if (cur < target_total) {
    missing <- target_total - cur
    ord <- order(sizes - k, decreasing = TRUE)
    i <- 1
    while (missing > 0) {
      cl <- names(k)[ord[i]]
      if (k[cl] < sizes[cl]) {
        k[cl] <- k[cl] + 1
        missing <- missing - 1
      }
      i <- i + 1
      if (i > length(ord)) i <- 1
      if (all(k >= sizes)) break
    }
  }

  selected <- character(0)

  for (cl in row_split_levels) {
    kk <- k[cl]
    if (is.na(kk) || kk <= 0) next

    genes <- rownames(heat_z)[cluster_vec[rownames(heat_z)] == cl]
    if (length(genes) == 0) next

    # approximate displayed order using the same dist+linkage as heatmap
    if (length(genes) >= 2) {
      submat <- heat_z[genes, , drop = FALSE]
      hc <- hclust(safe_cor_dist(submat), method = "ward.D2")
      genes_ord <- genes[hc$order]
    } else {
      genes_ord <- genes
    }

    # evenly spaced selection
    pos <- unique(pmax(1, pmin(length(genes_ord),
                              round(seq(1, length(genes_ord), length.out = kk)))))
    selected <- c(selected, genes_ord[pos])
  }

  selected <- unique(selected)
  at <- match(selected, rownames(heat_z))
  labels <- label_for_gene(selected)

  o <- order(at)
  list(at = at[o], labels = labels[o], genes = selected[o])
}

# ============================================================
# Read inputs
# ============================================================
samples <- read.delim(samples_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

required_cols <- c("sample_id", "condition")
missing_cols <- setdiff(required_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop("samples.tsv is missing columns: ", paste(missing_cols, collapse = ", "))
}

cts <- read.delim(counts_file, sep = "\t", header = TRUE, check.names = FALSE)
if (!("GeneID" %in% colnames(cts))) stop("Counts matrix must have column 'GeneID'.")

# Fix count column names if needed
count_cols <- colnames(cts)[colnames(cts) != "GeneID"]
if (length(setdiff(samples$sample_id, count_cols)) > 0) {
  cleaned <- clean_names(count_cols)
  if (length(setdiff(samples$sample_id, cleaned)) == 0) {
    message("[INFO] Applying clean_names() to count matrix columns.")
    colnames(cts)[colnames(cts) != "GeneID"] <- cleaned
  }
}

rownames(cts) <- cts$GeneID
cts <- cts[, setdiff(colnames(cts), "GeneID"), drop = FALSE]

missing_in_counts <- setdiff(samples$sample_id, colnames(cts))
if (length(missing_in_counts) > 0) {
  stop("These sample_id are in samples.tsv but NOT in counts matrix columns:\n",
       paste(missing_in_counts, collapse = "\n"))
}

# Align columns to samples.tsv order
cts <- cts[, samples$sample_id, drop = FALSE]

coldata <- samples %>%
  dplyr::select(sample_id, condition) %>%
  tibble::column_to_rownames("sample_id")
coldata$condition <- factor(coldata$condition)

# ============================================================
# Map condition names robustly + set order
# ============================================================
levs <- levels(coldata$condition)
cat("Condition levels detected:\n")
print(levs)

cond_control <- find_level("control", levs)
cond_diamide <- find_level("diamide", levs)
cond_noFe    <- find_level("noFe", levs)
cond_h2o2    <- find_level("h2o2", levs)

# fixed order for plotting
coldata$condition <- factor(
  as.character(coldata$condition),
  levels = c(cond_control, cond_diamide, cond_h2o2, cond_noFe)
)

# ============================================================
# Fit DESeq2 models
# ============================================================
cat("\nFitting DESeq2 with reference =", cond_control, "\n")
dds_controlRef <- make_dds(cts, coldata, cond_control)

cat("\nFitting DESeq2 with reference =", cond_noFe, "\n")
dds_noFeRef <- make_dds(cts, coldata, cond_noFe)

# ============================================================
# Normalized counts + VST
# ============================================================
norm_counts <- counts(dds_controlRef, normalized = TRUE)
readr::write_tsv(tibble::as_tibble(norm_counts, rownames = "gene_id"),
                 file.path(out_dir, "normalized_counts.tsv"))
message("Normalized counts written: ", file.path(out_dir, "normalized_counts.tsv"))

vsd <- vst(dds_controlRef, blind = FALSE)
vsd_mat <- assay(vsd)

plot_pca_pdf(vsd, coldata, file.path(fig_dir, "PCA_vst.pdf"))
message("Wrote: ", file.path(fig_dir, "PCA_vst.pdf"))

pdf(file.path(fig_dir, "Dispersion.pdf"), width = 7, height = 6)
plotDispEsts(dds_controlRef, main = "Dispersion estimates")
dev.off()
message("Wrote: ", file.path(fig_dir, "Dispersion.pdf"))

# ============================================================
# Contrasts
# ============================================================
comparisons <- list(
  diamide_vs_control = list(dds = dds_controlRef, num = cond_diamide, ref = cond_control,
                            out = file.path(out_dir, "DE_diamide_vs_control.tsv")),
  noFe_vs_control    = list(dds = dds_controlRef, num = cond_noFe,    ref = cond_control,
                            out = file.path(out_dir, "DE_noFe_vs_control.tsv")),
  h2o2_vs_noFe       = list(dds = dds_noFeRef,     num = cond_h2o2,    ref = cond_noFe,
                            out = file.path(out_dir, "DE_h2o2_vs_noFe.tsv"))
)

sig_counts <- tibble(contrast = character(), n_sig = integer())
res_store <- list()

for (nm in names(comparisons)) {
  spec <- comparisons[[nm]]

  cat("\nRunning: ", nm, " (", spec$num, " vs ", spec$ref, ")\n", sep = "")
  rr <- get_shrunken(spec$dds, spec$num, spec$ref)

  readr::write_tsv(rr$df, spec$out)
  cat("Wrote: ", spec$out, "\n", sep = "")

  sig_file <- file.path(out_dir, paste0("DE_", nm, ".sig.tsv"))
  n_sig <- write_sig_list(rr$df, sig_file, padj_cutoff = PADJ_CUTOFF, lfc_cutoff = LFC_CUTOFF)
  message("Wrote: ", sig_file, " (n_sig = ", n_sig, ")")
  sig_counts <- bind_rows(sig_counts, tibble(contrast = nm, n_sig = n_sig))

  ma_pdf <- file.path(fig_dir, paste0("MA_", nm, ".pdf"))
  plot_ma_pdf(rr$res_shr, ma_pdf, main = paste0(nm, " (apeglm)"))
  message("Wrote: ", ma_pdf)

  vol_pdf <- file.path(fig_dir, paste0("Volcano_", nm, ".pdf"))
  plot_volcano_pdf(rr$df, vol_pdf, title = nm, padj_cutoff = PADJ_CUTOFF, lfc_cutoff = LFC_CUTOFF)
  message("Wrote: ", vol_pdf)

  # Top50 ComplexHeatmap (row labels include gene symbol when present)
  top50 <- rr$df %>% filter(!is.na(padj)) %>% slice_head(n = 50) %>% pull(gene_id)
  top50 <- top50[top50 %in% rownames(vsd_mat)]
  if (length(top50) >= 2) {
    z <- zscore_rows(vsd_mat[top50, , drop = FALSE])

    meta2 <- coldata[colnames(z), , drop = FALSE]
    ord <- order(meta2$condition, rownames(meta2))
    z <- z[, ord, drop = FALSE]
    meta2 <- meta2[ord, , drop = FALSE]

    colAnn <- make_col_annotation(meta2)

    myCol <- colorRampPalette(c("royalblue", "white", "red3"))(100)
    myBreaks <- seq(-2, 2, length.out = 100)

    ht <- Heatmap(
      z,
      name = "Z-score",
      col = colorRamp2(myBreaks, myCol),
      top_annotation = colAnn,

      cluster_rows = TRUE,
      clustering_distance_rows = safe_cor_dist,
      clustering_method_rows = "ward.D2",

      show_row_names = TRUE,
      row_labels = label_for_gene(rownames(z)),
      row_names_gp = gpar(fontsize = 7),

      column_split = meta2$condition,
      cluster_columns = FALSE,
      show_column_names = TRUE,
      column_names_rot = 90,
      column_names_gp = gpar(fontsize = 8),

      use_raster = FALSE
    )

    out_pdf <- file.path(fig_dir, paste0("ComplexHeatmap_Top50_", nm, ".pdf"))
    pdf(out_pdf, width = 10, height = max(7, 0.18 * nrow(z) + 5))
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
    dev.off()
    message("Wrote: ", out_pdf)
  }

  res_store[[nm]] <- rr
}

readr::write_tsv(sig_counts, file.path(out_dir, "DEG_summary_counts.tsv"))
message("Wrote: ", file.path(out_dir, "DEG_summary_counts.tsv"))

# ============================================================
# UNION genes across contrasts
# ============================================================
sig_union <- unique(unlist(lapply(names(res_store), function(nm) {
  df <- res_store[[nm]]$df
  df %>%
    filter(!is.na(padj), padj < PADJ_CUTOFF,
           !is.na(log2FoldChange), abs(log2FoldChange) >= LFC_CUTOFF) %>%
    pull(gene_id)
})))
message(length(sig_union), " genes are DE in any contrast (union).")

# Cap union if huge by min padj across contrasts
if (length(sig_union) > MAX_UNION_GENES) {
  minpadj <- rep(1, length(sig_union)); names(minpadj) <- sig_union
  for (nm in names(res_store)) {
    df <- res_store[[nm]]$df %>% select(gene_id, padj)
    m <- setNames(df$padj, df$gene_id)
    common <- intersect(names(minpadj), names(m))
    minpadj[common] <- pmin(minpadj[common], m[common], na.rm = TRUE)
  }
  sig_union <- names(sort(minpadj))[1:MAX_UNION_GENES]
  message("Union capped to ", MAX_UNION_GENES, " genes by min(padj).")
}

sig_union <- sig_union[sig_union %in% rownames(vsd_mat)]
if (length(sig_union) < 2) stop("Union DEG list has <2 genes present in VST matrix.")

# Matrices
vst_union <- vsd_mat[sig_union, , drop = FALSE]
z_union   <- zscore_rows(vst_union)

# Order columns by condition then sample
meta2 <- coldata[colnames(z_union), , drop = FALSE]
ord <- order(meta2$condition, rownames(meta2))
z_union <- z_union[, ord, drop = FALSE]
vst_union <- vst_union[, ord, drop = FALSE]
meta2 <- meta2[ord, , drop = FALSE]

# PAM clusters for row split
k <- choose_pam_k(nrow(z_union))
if (!is.na(k) && k >= 2 && nrow(z_union) > k) {
  pamClusters <- cluster::pam(z_union, k = k)
  cluster_vec <- paste0("Cluster ", pamClusters$clustering)
} else {
  cluster_vec <- rep("Cluster 1", nrow(z_union))
}
names(cluster_vec) <- rownames(z_union)
row_split <- factor(cluster_vec, levels = sort(unique(cluster_vec)))

# ============================================================
# EXPORT CSVs (gene + values across conditions + cluster label)
# ============================================================
cond_levels <- levels(meta2$condition)

mean_by_condition <- function(mat, meta) {
  sapply(cond_levels, function(lv) {
    idx <- which(meta$condition == lv)
    if (length(idx) == 0) return(rep(NA_real_, nrow(mat)))
    rowMeans(mat[, idx, drop = FALSE])
  })
}

vst_means <- mean_by_condition(vst_union, meta2)

norm_union <- norm_counts[rownames(vst_union), colnames(vst_union), drop = FALSE]
norm_means <- mean_by_condition(norm_union, meta2)

log2norm_union <- log2(norm_union + 1)
log2norm_means <- mean_by_condition(log2norm_union, meta2)

z_means <- mean_by_condition(z_union, meta2)

out_means <- tibble(
  gene_id = rownames(vst_union),
  gene_symbol = symbol_for_gene(rownames(vst_union)),
  gene_label  = label_for_gene(rownames(vst_union)),
  cluster = unname(cluster_vec[rownames(vst_union)])
) %>%
  bind_cols(as_tibble(vst_means)      %>% rename_with(~paste0("VST_mean_", .x))) %>%
  bind_cols(as_tibble(norm_means)     %>% rename_with(~paste0("norm_mean_", .x))) %>%
  bind_cols(as_tibble(log2norm_means) %>% rename_with(~paste0("log2norm_mean_", .x))) %>%
  bind_cols(as_tibble(z_means)        %>% rename_with(~paste0("Z_mean_", .x)))

out_means_csv <- file.path(out_dir, "DEG_union_cluster_condition_means.csv")
readr::write_csv(out_means, out_means_csv)
message("Wrote: ", out_means_csv)

# Per-sample VST (union genes)
out_vst <- tibble(
  gene_id = rownames(vst_union),
  gene_symbol = symbol_for_gene(rownames(vst_union)),
  gene_label  = label_for_gene(rownames(vst_union)),
  cluster = unname(cluster_vec[rownames(vst_union)])
) %>%
  bind_cols(as_tibble(vst_union, rownames = NULL))

out_vst_csv <- file.path(out_dir, "DEG_union_clusters_vst_per_sample.csv")
readr::write_csv(out_vst, out_vst_csv)
message("Wrote: ", out_vst_csv)

# Per-sample Z-score (exact heatmap values)
out_z <- tibble(
  gene_id = rownames(z_union),
  gene_symbol = symbol_for_gene(rownames(z_union)),
  gene_label  = label_for_gene(rownames(z_union)),
  cluster = unname(cluster_vec[rownames(z_union)])
) %>%
  bind_cols(as_tibble(z_union, rownames = NULL))

out_z_csv <- file.path(out_dir, "DEG_union_clusters_zscore_per_sample.csv")
readr::write_csv(out_z, out_z_csv)
message("Wrote: ", out_z_csv)

# ============================================================
# UNION ComplexHeatmap (all conditions) with "head to toe" labels
# ============================================================
colAnn <- make_col_annotation(meta2)

boxplotCol <- HeatmapAnnotation(
  boxplot = anno_boxplot(
    z_union,
    border = FALSE,
    gp = gpar(fill = "#CCCCCC"),
    pch = ".",
    size = unit(1.5, "mm"),
    axis = TRUE,
    axis_param = list(gp = gpar(fontsize = 9), side = "left")
  ),
  which = "col",
  annotation_height = unit(1.2, "cm")
)

label_sel <- select_genes_to_label_fill_height(
  heat_z = z_union,
  cluster_vec = cluster_vec,
  row_split_levels = levels(row_split),
  target_total = LABEL_TARGET_TOTAL,
  max_total = LABEL_MAX_TOTAL,
  min_per_cluster = LABEL_MIN_PER_CLUSTER
)

genelabels <- rowAnnotation(
  Genes = anno_mark(
    at = label_sel$at,
    labels = label_sel$labels,
    labels_gp = gpar(fontsize = LABEL_FONTSIZE, fontface = "plain"),
    padding = unit(1, "mm")
  ),
  width = ComplexHeatmap::max_text_width(label_sel$labels, gp = gpar(fontsize = LABEL_FONTSIZE)) +
    unit(6, "mm")
)

myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
myBreaks <- seq(-4, 4, length.out = 100)

# Safety: if any condition has <2 replicates, disable column clustering
min_reps <- min(as.integer(table(meta2$condition)))
cluster_cols_flag <- is.finite(min_reps) && (min_reps >= 2)

ht_union <- Heatmap(
  z_union,
  name = "Gene\nZ-score",
  col = colorRamp2(myBreaks, myCol),

  split = row_split,
  row_gap = unit(2.5, "mm"),
  cluster_row_slices = FALSE,

  top_annotation = colAnn,
  bottom_annotation = boxplotCol,

  cluster_rows = TRUE,
  clustering_distance_rows = safe_cor_dist,
  clustering_method_rows = "ward.D2",
  show_row_dend = TRUE,
  row_dend_width = unit(25, "mm"),
  show_row_names = FALSE,

  column_split = meta2$condition,
  column_gap = unit(2.5, "mm"),
  cluster_columns = cluster_cols_flag,
  cluster_column_slices = FALSE,
  show_column_dend = cluster_cols_flag,
  show_column_names = TRUE,
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 8),

  # FIX: safe distance (prevents NA/NaN/Inf crash)
  clustering_distance_columns = safe_cor_dist,
  clustering_method_columns = "ward.D2",

  use_raster = nrow(z_union) > 300
)

union_pdf <- file.path(fig_dir, "ComplexHeatmap_DEG_union_all_conditions.pdf")
pdf(union_pdf, width = 14, height = .height_from_rows(nrow(z_union)))
draw(
  ht_union + genelabels,
  heatmap_legend_side = "left",
  annotation_legend_side = "right",
  row_sub_title_side = "left"
)
dev.off()
message("Wrote: ", union_pdf)

# Session info
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Wrote: ", file.path(out_dir, "sessionInfo.txt"))

message("\nDONE\n- Tables: ", out_dir,
        "\n- Figures: ", fig_dir, "\n")
