# Shared publication aesthetics for comparative genomics / longevity plots.
# Aligned with birds_bats_and_mammals/scripts/shared/plot_aesthetics.R and
# comparative_reprogramming/scripts/_project.R.
#
#   source("aesthetics.R")
#
# Signed heatmaps always use ColorBrewer RdBu ends (never purple/green).

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 required for aesthetics.R")
}
if (!requireNamespace("stringr", quietly = TRUE)) {
  stop("stringr required for aesthetics.R")
}

# -----------------------------------------------------------------------------
# Core signed / categorical colors
# -----------------------------------------------------------------------------

# Up / down fills for bars, forests, volcanoes (iPSC / fibro style)
col_up <- "#9DBA91"
col_down <- "#A391BA"
col_concordant <- "#009E73"
col_delta <- "#E69F00"

# ColorBrewer RdBu (signed heatmaps)
gradient_low <- "#2166AC"
gradient_high <- "#B2182B"

PUB_BASE <- 16L
PUB_DPI <- 300L

sig_caption <- "^ FDR<0.1, * FDR<0.05, ** FDR<0.01, *** FDR<0.001"

# -----------------------------------------------------------------------------
# Optional named palettes (use when the project has that axis)
# -----------------------------------------------------------------------------

order_colors <- c(
  Primates = "#9b383a",
  Rodentia = "#31a6ad",
  Artiodactyla = "#D99F6A",
  Chiroptera = "#6755A3",
  Perissodactyla = "#E69F00",
  Carnivora = "#C48B9F",
  Didelphimorphia = "#808080"
)

celltype_colors <- c(
  Fibroblasts = "#A391BA",
  iPSCs = "#9DBA91",
  ESCs = "#009E73"
)

tissue_colors <- c(
  Liver = "#9b383a",
  Brain = "#31a6ad",
  Kidney = "#D99F6A",
  Heart = "#CC5454",
  Blood = "#8B0000",
  Unknown = "#7A7A7A"
)

# Bird / bat / mammal clade colors (birds_bats_and_mammals)
bird_color <- "#6B8E9F"
mammal_color <- "#7CB577"
all_mammal_color <- "#336600"
bat_color <- "#B68AB9"
flightless_bird_color <- "#3F545E"

clade_group_colors <- c(
  "Flying birds" = bird_color,
  "Flightless birds" = flightless_bird_color,
  "Bats" = bat_color,
  "Non-flying mammals" = mammal_color,
  "All mammals" = all_mammal_color
)

# Longevity trait accents
ml_color <- "#C76F84"      # logML
mlres_color <- "#7671A3"   # MLres
both_color <- "#111111"

# Compact jewel palette for arbitrary categorical axes
jewel_palette <- c(
  "#2c6e5f", "#9b383a", "#4878a0", "#64a590", "#945a87",
  "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
  "#20854c", "#9c534a", "#31a6ad", "#f0a0a3", "#c57251",
  "#856890", "#7e9960", "#7d9a9e", "#B65B32", "#c07a7a",
  "#8a7c6d", "#AA4839", "#3D7B68"
)

get_jewel_palette <- function(n) {
  n <- as.integer(n)[1]
  if (is.na(n) || n < 1) stop("n must be a positive integer")
  if (n <= length(jewel_palette)) {
    jewel_palette[seq_len(n)]
  } else {
    rep(jewel_palette, length.out = n)
  }
}

# -----------------------------------------------------------------------------
# Theme + scales + save
# -----------------------------------------------------------------------------

theme_pub <- function(base_size = PUB_BASE, title_color_override = NULL) {
  tc <- if (is.null(title_color_override)) "grey10" else title_color_override
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = base_size + 4, hjust = 0.5, color = tc
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 2, hjust = 0.5, color = "grey30"
      ),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 2),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.title = ggplot2::element_text(size = base_size - 1, face = "bold"),
      plot.caption = ggplot2::element_text(
        size = base_size - 4, hjust = 1, color = "grey40"
      ),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 1),
      strip.background = ggplot2::element_rect(fill = "white"),
      plot.margin = ggplot2::margin(10, 15, 10, 10)
    )
}

# Alias used in comparative_reprogramming
theme_reprog <- theme_pub

scale_fill_heatmap <- function(...) {
  ggplot2::scale_fill_gradient2(
    low = gradient_low, mid = "white", high = gradient_high,
    midpoint = 0, ...
  )
}

scale_color_heatmap <- function(...) {
  ggplot2::scale_color_gradient2(
    low = gradient_low, mid = "white", high = gradient_high,
    midpoint = 0, ...
  )
}

heatmap_palette <- function(n = 100L) {
  grDevices::colorRampPalette(c(gradient_low, "white", gradient_high))(n)
}

save_plot <- function(p, stub, width = 10, height = 8, dir = ".", dpi = PUB_DPI) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  pdf_path <- file.path(dir, paste0(stub, ".pdf"))
  png_path <- file.path(dir, paste0(stub, ".png"))
  ggplot2::ggsave(pdf_path, p, width = width, height = height)
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = dpi, bg = "white")
  invisible(list(pdf = pdf_path, png = png_path))
}

# -----------------------------------------------------------------------------
# Pathway labels + significance
# -----------------------------------------------------------------------------

.pathway_prefix_re <- paste0(
  "^(REACTOME_|KEGG_|WP_|BIOCARTA_|PID_|GO_|GOBP_|GOCC_|GOMF_|",
  "HP_|NABA_|HALLMARK_)"
)

clean_pathway <- function(name, wrap_width = 40, aliases = NULL) {
  sapply(as.character(name), function(raw) {
    used_alias <- FALSE
    if (!is.null(aliases) && raw %in% names(aliases)) {
      cleaned <- unname(aliases[[raw]])
      used_alias <- TRUE
    } else {
      cleaned <- gsub(.pathway_prefix_re, "", raw)
      cleaned <- tolower(gsub("_", " ", cleaned))
      cleaned <- trimws(cleaned)
      if (!is.null(aliases) && cleaned %in% names(aliases)) {
        cleaned <- unname(aliases[[cleaned]])
        used_alias <- TRUE
      } else if (nzchar(cleaned)) {
        cleaned <- paste0(
          toupper(substr(cleaned, 1, 1)),
          substr(cleaned, 2, nchar(cleaned))
        )
      }
    }
    ww <- if (used_alias) max(wrap_width, nchar(cleaned)) else wrap_width
    stringr::str_wrap(cleaned, width = ww)
  }, USE.NAMES = FALSE)
}

# Alias used in comparative_reprogramming
clean_pathway_name <- function(name, wrap_width = 42) {
  clean_pathway(name, wrap_width = wrap_width)
}

pathway_dedup_key <- function(name) {
  toupper(gsub(.pathway_prefix_re, "", as.character(name)))
}

get_sig_star <- function(padj) {
  out <- rep("", length(padj))
  out[!is.na(padj) & padj < 0.1] <- "^"
  out[!is.na(padj) & padj < 0.05] <- "*"
  out[!is.na(padj) & padj < 0.01] <- "**"
  out[!is.na(padj) & padj < 0.001] <- "***"
  out
}

# Alias
fgsea_sig_star <- get_sig_star

################################################################################
# Example usage (does not run on source)
################################################################################
if (FALSE) {
  library(ggplot2)

  df <- data.frame(
    x = c(-1.2, -0.4, 0.3, 1.1),
    y = c("A", "B", "C", "D"),
    dir = c("down", "down", "up", "up")
  )

  p <- ggplot(df, aes(x, y, fill = dir)) +
    geom_col() +
    scale_fill_manual(values = c(up = col_up, down = col_down)) +
    labs(title = "Example", fill = NULL) +
    theme_pub()

  save_plot(p, "example_bars", width = 6, height = 4, dir = "path/to/plots")

  clean_pathway(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_CELL_CYCLE"))
  get_sig_star(c(0.0001, 0.02, 0.08, 0.2))
}
