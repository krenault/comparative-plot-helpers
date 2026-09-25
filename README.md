# comparative-plot-helpers

Reusable R plotting helpers for comparative genomics / longevity projects.

| File | Function / exports | Purpose |
|------|--------------------|---------|
| `aesthetics.R` | `theme_pub` / `theme_reprog`, colors, `save_plot`, `clean_pathway`, RdBu scales | Shared publication skin used across projects |
| `Tree.R` | `plot_phylogenetic_tree()` | ggtree plots with tip colors, shapes, and bar annotations |
| `Enrichment_plot.R` | `plot_go_enrichment()` | GO / pathway enrichment bubble plots |
| `Heatmap.R` | `plot_custom_heatmap()` | Multi-dataset heatmaps (ComplexHeatmap or ggplot2) |

## Usage

```r
source("aesthetics.R")   # theme, colors, save_plot — source first
source("Tree.R")
source("Enrichment_plot.R")
source("Heatmap.R")
```

Each plot file ends with an `if (FALSE) { ... }` example block. Sourcing does **not** run examples.

### Aesthetics quick reference

```r
source("aesthetics.R")

# theme_bw + centered bold titles (alias: theme_reprog)
p + theme_pub()

# Signed heatmap fill (ColorBrewer RdBu)
p + scale_fill_heatmap()

# PDF + PNG (300 dpi, white background)
save_plot(p, "fig1_volcano", width = 8, height = 6, dir = "plots")

# Pathway label cleanup
clean_pathway("HALLMARK_OXIDATIVE_PHOSPHORYLATION")
```

Shared constants: `col_up`, `col_down`, `gradient_low`, `gradient_high`, `order_colors`, `celltype_colors`, `tissue_colors`, `clade_group_colors`.

## Dependencies

**aesthetics.R:** `ggplot2`, `stringr`

**Tree.R:** `ape`, `ggtree`, `ggtreeExtra`, `ggplot2`, `dplyr`, `ggnewscale`

**Enrichment_plot.R:** `ggplot2`, `dplyr`, `tidyr`, `scales`, `stringr`

**Heatmap.R:** `ggplot2`, `dplyr`, `pheatmap`, `ComplexHeatmap`, `circlize`, `RColorBrewer`, `viridis`, `reshape2`, `scales`, `tidyr`, `tibble`

## License

MIT
