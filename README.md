# comparative-plot-helpers

Reusable R plotting helpers for comparative genomics / longevity projects.

| File | Function | Purpose |
|------|----------|---------|
| `Tree.R` | `plot_phylogenetic_tree()` | ggtree plots with optional tip colors, shapes, and bar annotations |
| `Enrichment_plot.R` | `plot_go_enrichment()` | GO / pathway enrichment bubble plots (single or multi-dataset) |
| `Heatmap.R` | `plot_custom_heatmap()` | Multi-dataset heatmaps via ComplexHeatmap or ggplot2 |

## Usage

```r
source("Tree.R")
source("Enrichment_plot.R")
source("Heatmap.R")
```

Each file ends with an `if (FALSE) { ... }` example block — edit the placeholder paths and run interactively. Sourcing a file does **not** execute the examples.

## Dependencies

**Shared:** `ggplot2`, `dplyr`

**Tree.R:** `ape`, `ggtree`, `ggtreeExtra`, `ggnewscale`

**Enrichment_plot.R:** `tidyr`, `scales`, `stringr`

**Heatmap.R:** `pheatmap`, `ComplexHeatmap`, `circlize`, `RColorBrewer`, `viridis`, `reshape2`, `scales`, `tidyr`, `tibble`

## License

MIT
