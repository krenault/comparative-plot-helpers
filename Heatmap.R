# Multi-dataset enrichment / value heatmaps (ComplexHeatmap or ggplot2).
# Source this file, then call plot_custom_heatmap(...).

library(ggplot2)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(scales)
library(tidyr)
library(tibble)

plot_custom_heatmap <- function(
    data_file = NULL, data_dir = NULL, data_pattern = NULL, data = NULL, file_list = NULL, dataset_names = NULL,
    dataset_groups = NULL, group_colors = NULL, show_dataset_annotation = TRUE, show_group_colors = TRUE, 
    show_dataset_colors = FALSE,  # Changed default to FALSE
    annotation_height = 0.5, row_column = NULL, value_column = NULL, pvalue_column = NULL,
    min_datasets = 1, max_elements = NULL, selection_criteria = "value", selection_direction = "high", remove_na_rows = TRUE,
    sig_levels = c(0.001, 0.01, 0.05), sig_symbols = c("***", "**", "*"), show_significance = TRUE, sig_size = 3, sig_color = "white",
    color_palette = "viridis", custom_colors = NULL, na_color = "grey90", cluster_rows = TRUE, cluster_cols = TRUE, 
    scale_data = "none", cell_border_color = "white", cell_border_width = 0.5, show_cell_borders = TRUE, 
    family_groups = NULL, family_colors = NULL, show_cell_values = FALSE, cell_value_size = 2.5,
    cell_value_color = "black", cell_width = NULL, cell_height = NULL, main_title = NULL, subtitle = NULL,
    x_label = NULL, y_label = NULL, legend_title = NULL, legend_position = NULL, row_annotation = NULL,
    col_annotation = NULL, row_ann_colors = NULL, col_ann_colors = NULL, figure_width = 12, figure_height = 10,
    heatmap_width = NULL, heatmap_height = NULL, replace_underscores = TRUE, text_columns = c("row", "column"),
    row_label_size = 10, col_label_size = 10, title_size = 14, legend_text_size = 10, rotate_col_labels = 45,
    rotate_row_labels = 0, col_label_hjust = NULL, col_label_vjust = NULL, row_label_hjust = NULL,
    row_label_vjust = NULL, output_dir = NULL, output_prefix = NULL, save_formats = c("pdf", "png"),
    dpi = 300, heatmap_engine = "ComplexHeatmap", return_plot_data = FALSE, show_row_names = TRUE, verbose = TRUE) {
  
  if (verbose) cat("=== STARTING MULTI-DATASET HEATMAP GENERATION ===\n")
  
  get_jewel_palette <- function(n) {
    jewel_colors <- c("#2c6e5f", "#9b383a", "#4878a0", "#64a590", "#945a87", "#b68ab9", "#9c6e5a", "#9fbd8b",
                      "#a85472", "#86acb9", "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3", "#c57251",
                      "#856890", "#7e9960", "#7d9a9e", "#B65B32", "#c07a7a", "#8a7c6d", "#AA4839", "#4878a0", "#3D7B68")
    if (n <= length(jewel_colors)) return(jewel_colors[1:n]) else return(rep(jewel_colors, ceiling(n/length(jewel_colors)))[1:n])  # Use jewel palette instead of rainbow
  }
  
  load_multi_dataset <- function() {
    datasets <- list()
    if (!is.null(data)) {
      if (verbose) cat("Using provided data object\n")
      datasets[["provided_data"]] <- data
    } else if (!is.null(file_list)) {
      if (verbose) cat("Loading data from", length(file_list), "files\n")
      for (i in seq_along(file_list)) {
        file_path <- file_list[i]
        dataset_name <- if (!is.null(dataset_names) && length(dataset_names) >= i) {
          dataset_names[i]
        } else {
          tools::file_path_sans_ext(basename(file_path))
        }
        if (tools::file_ext(file_path) == "csv") {
          datasets[[dataset_name]] <- read.csv(file_path, stringsAsFactors = FALSE)
        } else {
          datasets[[dataset_name]] <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        }
        if (verbose) cat("  Loaded", dataset_name, "with", nrow(datasets[[dataset_name]]), "rows\n")
      }
    } else if (!is.null(data_dir) && !is.null(data_pattern)) {
      if (verbose) cat("Loading data from directory:", data_dir, "with pattern:", data_pattern, "\n")
      files <- list.files(data_dir, pattern = data_pattern, full.names = TRUE)
      if (length(files) == 0) stop("No files found matching pattern")
      for (i in seq_along(files)) {
        file_path <- files[i]
        dataset_name <- if (!is.null(dataset_names) && length(dataset_names) >= i) {
          dataset_names[i]
        } else {
          tools::file_path_sans_ext(basename(file_path))
        }
        if (tools::file_ext(file_path) == "csv") {
          datasets[[dataset_name]] <- read.csv(file_path, stringsAsFactors = FALSE)
        } else {
          datasets[[dataset_name]] <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        }
        if (verbose) cat("  Loaded", dataset_name, "with", nrow(datasets[[dataset_name]]), "rows\n")
      }
    } else if (!is.null(data_file)) {
      if (verbose) cat("Loading data from single file:", data_file, "\n")
      if (tools::file_ext(data_file) == "csv") {
        datasets[["single_file"]] <- read.csv(data_file, stringsAsFactors = FALSE)
      } else {
        datasets[["single_file"]] <- read.table(data_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
      }
    } else {
      stop("Must provide either 'data', 'data_file', 'file_list', or both 'data_dir' and 'data_pattern'")
    }
    return(datasets)
  }
  
  datasets <- load_multi_dataset()
  
  if (verbose) cat("Combining", length(datasets), "datasets...\n")
  combined_data <- NULL
  for (dataset_name in names(datasets)) {
    dataset <- datasets[[dataset_name]]
    dataset$dataset_source <- dataset_name
    if (is.null(combined_data)) {
      combined_data <- dataset
    } else {
      combined_data <- rbind(combined_data, dataset)
    }
  }
  
  if (verbose) cat("Combined data has", nrow(combined_data), "rows\n")
  
  required_cols <- c(row_column, value_column)
  missing_cols <- setdiff(required_cols, colnames(combined_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (verbose) cat("Processing data and handling duplicates...\n")
  heatmap_data_clean <- combined_data %>%
    group_by(!!sym(row_column), dataset_source) %>%
    summarise(
      value = mean(!!sym(value_column), na.rm = TRUE),
      pvalue = if (!is.null(pvalue_column) && pvalue_column %in% colnames(combined_data)) {
        mean(!!sym(pvalue_column), na.rm = TRUE)
      } else {
        NA
      },
      .groups = "drop"
    )
  
  value_matrix <- heatmap_data_clean %>%
    dplyr::select(!!sym(row_column), dataset_source, value) %>%
    pivot_wider(names_from = dataset_source, values_from = value) %>%
    column_to_rownames(row_column) %>%
    as.matrix()
  
  pvalue_matrix <- NULL
  if (!is.null(pvalue_column) && pvalue_column %in% colnames(combined_data)) {
    pvalue_matrix <- heatmap_data_clean %>%
      dplyr::select(!!sym(row_column), dataset_source, pvalue) %>%
      pivot_wider(names_from = dataset_source, values_from = pvalue) %>%
      column_to_rownames(row_column) %>%
      as.matrix()
  }
  
  if (verbose) cat("Initial matrix created with", nrow(value_matrix), "rows and", ncol(value_matrix), "columns\n")
  
  if (verbose) cat("Applying filtering and selection criteria...\n")
  
  if (min_datasets > 1) {
    datasets_per_row <- rowSums(!is.na(value_matrix))
    keep_rows <- datasets_per_row >= min_datasets
    value_matrix <- value_matrix[keep_rows, , drop = FALSE]
    if (!is.null(pvalue_matrix)) {
      pvalue_matrix <- pvalue_matrix[keep_rows, , drop = FALSE]
    }
    if (verbose) cat("After min_datasets filter (>= ", min_datasets, "):", nrow(value_matrix), "rows remaining\n")
  }
  
  if (remove_na_rows) {
    all_na_rows <- rowSums(!is.na(value_matrix)) == 0
    if (any(all_na_rows)) {
      value_matrix <- value_matrix[!all_na_rows, , drop = FALSE]
      if (!is.null(pvalue_matrix)) {
        pvalue_matrix <- pvalue_matrix[!all_na_rows, , drop = FALSE]
      }
      if (verbose) cat("Removed", sum(all_na_rows), "rows with all NA values\n")
    }
  }
  
  if (!is.null(max_elements) && nrow(value_matrix) > max_elements) {
    if (verbose) cat("Selecting top", max_elements, "elements based on", selection_criteria, "...\n")
    if (selection_criteria == "frequency") {
      frequency_scores <- rowSums(!is.na(value_matrix))
      selected_indices <- order(frequency_scores, decreasing = TRUE)[1:max_elements]
    } else if (selection_criteria == "significance" && !is.null(pvalue_matrix)) {
      sig_counts <- rowSums(pvalue_matrix <= 0.05, na.rm = TRUE)
      selected_indices <- order(sig_counts, decreasing = TRUE)[1:max_elements]
    } else if (selection_criteria == "pvalue" && !is.null(pvalue_matrix)) {
      min_pvalues <- apply(pvalue_matrix, 1, min, na.rm = TRUE)
      if (selection_direction == "low") {
        selected_indices <- order(min_pvalues, decreasing = FALSE)[1:max_elements]
      } else {
        selected_indices <- order(min_pvalues, decreasing = TRUE)[1:max_elements]
      }
    } else {
      if (selection_direction == "high") {
        max_values <- apply(abs(value_matrix), 1, max, na.rm = TRUE)
        selected_indices <- order(max_values, decreasing = TRUE)[1:max_elements]
      } else {
        min_values <- apply(abs(value_matrix), 1, min, na.rm = TRUE)
        selected_indices <- order(min_values, decreasing = FALSE)[1:max_elements]
      }
    }
    value_matrix <- value_matrix[selected_indices, , drop = FALSE]
    if (!is.null(pvalue_matrix)) {
      pvalue_matrix <- pvalue_matrix[selected_indices, , drop = FALSE]
    }
    if (verbose) cat("Selected", nrow(value_matrix), "elements for plotting\n")
  }
  
  if (replace_underscores) {
    if ("row" %in% text_columns) {
      rownames(value_matrix) <- gsub("_", " ", rownames(value_matrix))
      if (!is.null(pvalue_matrix)) {
        rownames(pvalue_matrix) <- gsub("_", " ", rownames(pvalue_matrix))
      }
      if (verbose) cat("Replaced underscores with spaces in row names\n")
    }
    if ("column" %in% text_columns) {
      colnames(value_matrix) <- gsub("_", " ", colnames(value_matrix))
      if (!is.null(pvalue_matrix)) {
        colnames(pvalue_matrix) <- gsub("_", " ", colnames(pvalue_matrix))
      }
      if (verbose) cat("Replaced underscores with spaces in column names\n")
    }
  }
  
  if (scale_data == "row") {
    value_matrix <- t(scale(t(value_matrix)))
  } else if (scale_data == "column") {
    value_matrix <- scale(value_matrix)
  } else if (scale_data == "both") {
    value_matrix <- scale(scale(value_matrix))
  }
  
  if (verbose) cat("Matrix created with", nrow(value_matrix), "rows and", ncol(value_matrix), "columns\n")
  
  if (!is.null(custom_colors)) {
    colors <- custom_colors
  } else if (color_palette == "viridis") {
    colors <- viridis::viridis(100)
  } else if (color_palette == "plasma") {
    colors <- viridis::plasma(100)
  } else if (color_palette == "inferno") {
    colors <- viridis::inferno(100)
  } else if (color_palette == "magma") {
    colors <- viridis::magma(100)
  } else if (color_palette == "RdYlBu") {
    colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  } else if (color_palette == "RdBu") {
    colors <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  } else if (color_palette == "jewel") {
    colors <- get_jewel_palette(100)
  } else {
    colors <- viridis::viridis(100)
  }
  
  dataset_annotation <- NULL
  dataset_ann_colors <- NULL
  
  if (show_dataset_annotation) {
    if (verbose) cat("Creating dataset annotations...\n")
    column_info <- data.frame(
      column_name = colnames(value_matrix),
      dataset = colnames(value_matrix),
      stringsAsFactors = FALSE
    )
    
    dataset_annotation <- data.frame(
      row.names = column_info$column_name,
      stringsAsFactors = FALSE
    )
    
    if (!is.null(dataset_groups)) {
      dataset_annotation$Group <- dataset_groups[column_info$dataset]
      dataset_annotation$Group[is.na(dataset_annotation$Group)] <- "Ungrouped"
    }
    
    if (!is.null(family_groups)) {
      dataset_annotation$Family <- family_groups[column_info$dataset]
      dataset_annotation$Family[is.na(dataset_annotation$Family)] <- "Other"
    }
    
    unique_groups <- if (!is.null(dataset_groups)) unique(dataset_annotation$Group) else NULL
    unique_families <- if (!is.null(family_groups)) unique(dataset_annotation$Family) else NULL
    
    dataset_ann_colors <- list()
    
    # Removed dataset color section
    
    if (show_group_colors && !is.null(unique_groups)) {
      if (!is.null(group_colors)) {
        group_color_palette <- group_colors
        missing_groups <- setdiff(unique_groups, names(group_colors))
        if (length(missing_groups) > 0) {
          additional_colors <- get_jewel_palette(length(missing_groups))
          names(additional_colors) <- missing_groups
          group_color_palette <- c(group_color_palette, additional_colors)
        }
      } else {
        group_color_palette <- get_jewel_palette(length(unique_groups))
        names(group_color_palette) <- unique_groups
      }
      dataset_ann_colors$Group <- group_color_palette
    }
    
    if (!is.null(unique_families)) {
      if (!is.null(family_colors)) {
        family_color_palette <- family_colors
        missing_families <- setdiff(unique_families, names(family_colors))
        if (length(missing_families) > 0) {
          additional_colors <- get_jewel_palette(length(missing_families))
          names(additional_colors) <- missing_families
          family_color_palette <- c(family_color_palette, additional_colors)
        }
      } else {
        family_color_palette <- get_jewel_palette(length(unique_families))
        names(family_color_palette) <- unique_families
      }
      dataset_ann_colors$Family <- family_color_palette
    }
  }
  
  sig_matrix <- NULL
  if (show_significance && !is.null(pvalue_matrix)) {
    if (verbose) cat("Creating significance annotations...\n")
    sig_matrix <- matrix("", nrow = nrow(pvalue_matrix), ncol = ncol(pvalue_matrix))
    rownames(sig_matrix) <- rownames(pvalue_matrix)
    colnames(sig_matrix) <- colnames(pvalue_matrix)
    for (i in seq_along(sig_levels)) {
      sig_matrix[pvalue_matrix <= sig_levels[i] & !is.na(pvalue_matrix)] <- sig_symbols[i]
    }
    if (verbose) {
      for (i in seq_along(sig_levels)) {
        count <- sum(pvalue_matrix <= sig_levels[i] & !is.na(pvalue_matrix))
        cat("Cells with p <=", sig_levels[i], "(", sig_symbols[i], "):", count, "\n")
      }
    }
  }
  
  if (heatmap_engine == "ComplexHeatmap") {
    if (verbose) cat("Creating ComplexHeatmap...\n")
    col_fun <- circlize::colorRamp2(
      seq(min(value_matrix, na.rm = TRUE), max(value_matrix, na.rm = TRUE), length = length(colors)),
      colors
    )
    
    cell_fun <- NULL
    if (show_significance && !is.null(sig_matrix) && nrow(value_matrix) <= 500) {
      cell_fun <- function(j, i, x, y, width, height, fill) {
        if (!is.na(sig_matrix[i, j]) && sig_matrix[i, j] != "") {
          grid::grid.text(sig_matrix[i, j], x, y,
                          gp = grid::gpar(col = sig_color, fontsize = sig_size * 3, fontface = "bold"))
        }
        if (show_cell_values) {
          grid::grid.text(round(value_matrix[i, j], 2), x, y - grid::unit(0.15, "npc"),
                          gp = grid::gpar(col = cell_value_color, fontsize = cell_value_size * 2.5))
        }
      }
    } else if (show_cell_values && nrow(value_matrix) <= 500) {
      cell_fun <- function(j, i, x, y, width, height, fill) {
        grid::grid.text(round(value_matrix[i, j], 2), x, y,
                        gp = grid::gpar(col = cell_value_color, fontsize = cell_value_size * 3))
      }
    } else if (nrow(value_matrix) > 500) {
      if (verbose) cat("Large matrix detected (", nrow(value_matrix), " rows). Disabling cell annotations for performance.\n")
      cell_fun <- NULL
    }
    
    top_annotation <- NULL
    if (!is.null(dataset_annotation)) {
      if ("Family" %in% colnames(dataset_annotation)) {
        unique_families <- unique(dataset_annotation$Family)
        family_color_palette <- get_jewel_palette(length(unique_families))
        names(family_color_palette) <- unique_families
        dataset_ann_colors$Family <- family_color_palette
      }
    }
    
    if (!is.null(dataset_annotation) && ncol(dataset_annotation) > 0) {
      top_annotation <- ComplexHeatmap::HeatmapAnnotation(
        df = dataset_annotation,
        col = dataset_ann_colors,
        annotation_height = grid::unit(annotation_height, "cm"),
        show_annotation_name = TRUE,
        show_legend = legend_position != "none",
        annotation_name_gp = grid::gpar(fontsize = legend_text_size),
        annotation_legend_param = list(
          title_gp = grid::gpar(fontsize = legend_text_size + 1),
          labels_gp = grid::gpar(fontsize = legend_text_size),
          ncol = 4
        )
      )
    }
    
    final_plot <- ComplexHeatmap::Heatmap(
      value_matrix,
      col = col_fun,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_cols,
      show_row_names = if (is.logical(show_row_names)) show_row_names else (nrow(value_matrix) <= 100),
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = row_label_size),
      column_names_gp = grid::gpar(fontsize = col_label_size),
      row_names_rot = rotate_row_labels,
      column_names_rot = rotate_col_labels,
      column_title = main_title,
      column_title_gp = grid::gpar(fontsize = title_size, fontface = "bold"),
      top_annotation = top_annotation,
      heatmap_legend_param = list(
        title = ifelse(is.null(legend_title), value_column, legend_title),
        title_gp = grid::gpar(fontsize = legend_text_size + 1),
        labels_gp = grid::gpar(fontsize = legend_text_size),
        ncol = 4,
        legend_direction = ifelse(legend_position %in% c("top", "bottom"), "horizontal", "vertical")
      ),
      cell_fun = cell_fun,
      na_col = na_color,
      border = ifelse(show_cell_borders && nrow(value_matrix) <= 200, cell_border_color, FALSE),
      rect_gp = if (show_cell_borders && nrow(value_matrix) <= 200) {
        grid::gpar(col = cell_border_color, lwd = cell_border_width)
      } else {
        grid::gpar(col = NA)
      },
      width = if (!is.null(heatmap_width)) {
        grid::unit(heatmap_width, "cm")
      } else if (!is.null(cell_width)) {
        grid::unit(cell_width * ncol(value_matrix), "mm")
      } else {
        NULL
      },
      height = if (!is.null(heatmap_height)) {
        grid::unit(heatmap_height, "cm")
      } else if (!is.null(cell_height)) {
        grid::unit(cell_height * nrow(value_matrix), "mm")
      } else {
        NULL
      },
      show_heatmap_legend = legend_position != "none"
    )
  }
  
  if (heatmap_engine == "ComplexHeatmap") {
    draw_options <- list()
    if (legend_position == "none") {
      draw_options$heatmap_legend_list <- list()
      draw_options$annotation_legend_list <- list()
    }
    if (legend_position == "left") {
      draw_options$heatmap_legend_side = "left"
      draw_options$annotation_legend_side = "left"
    } else if (legend_position == "top") {
      draw_options$heatmap_legend_side = "top"
      draw_options$annotation_legend_side = "top"
    } else if (legend_position == "bottom") {
      draw_options$heatmap_legend_side = "bottom"
      draw_options$annotation_legend_side = "bottom"
    } else if (legend_position == "right") {
      draw_options$heatmap_legend_side = "right"
      draw_options$annotation_legend_side = "right"
    }
    do.call(ComplexHeatmap::draw, c(list(final_plot), draw_options))
  } else {
    print(final_plot)
  }
  
  if (!is.null(output_dir) && !is.null(output_prefix)) {
    if (verbose) cat("Saving outputs...\n")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    for (format in save_formats) {
      output_file <- file.path(output_dir, paste0(output_prefix, "_multi_heatmap.", format))
      if (format == "pdf") {
        pdf(output_file, width = figure_width, height = figure_height)
      } else if (format == "png") {
        png(output_file, width = figure_width * dpi/72, height = figure_height * dpi/72, res = dpi)
      }
      if (heatmap_engine == "ComplexHeatmap") {
        ComplexHeatmap::draw(final_plot)
      } else {
        print(final_plot)
      }
      dev.off()
      if (verbose) cat("Saved:", output_file, "\n")
    }
  }
  
  if (verbose) {
    cat("\n=== MULTI-DATASET HEATMAP SUMMARY ===\n")
    cat("Datasets combined:", length(datasets), "\n")
    cat("Dataset names:", paste(names(datasets), collapse = ", "), "\n")
    cat("Final matrix size:", nrow(value_matrix), "x", ncol(value_matrix), "\n")
    cat("Value range:", round(min(value_matrix, na.rm = TRUE), 3), "to", round(max(value_matrix, na.rm = TRUE), 3), "\n")
    cat("Missing values:", sum(is.na(value_matrix)), "\n")
    if (min_datasets > 1) {
      cat("Filtered to elements appearing in >=", min_datasets, "datasets\n")
    }
    if (!is.null(max_elements)) {
      cat("Selected top", max_elements, "elements based on", selection_criteria, "(", selection_direction, ")\n")
    }
    if (replace_underscores) {
      cat("Text formatting: replaced underscores in", paste(text_columns, collapse = ", "), "names\n")
    }
    if (show_dataset_annotation) {
      annotation_info <- c()
      if (show_group_colors && !is.null(dataset_groups)) annotation_info <- c(annotation_info, "group colors")
      if (length(annotation_info) > 0) {
        cat("Annotations displayed:", paste(annotation_info, collapse = ", "), "\n")
      }
    }
    if (!is.null(dataset_groups)) {
      cat("Dataset groups:\n")
      for (group in unique(dataset_groups)) {
        datasets_in_group <- names(dataset_groups)[dataset_groups == group]
        cat("  ", group, ":", paste(datasets_in_group, collapse = ", "), "\n")
      }
    }
    if (!is.null(pvalue_matrix)) {
      cat("Significant cells (p < 0.05):", sum(pvalue_matrix < 0.05, na.rm = TRUE), "\n")
      for (i in seq_along(sig_levels)) {
        cat("Cells with", sig_symbols[i], "(p <", sig_levels[i], "):",
            sum(pvalue_matrix <= sig_levels[i], na.rm = TRUE), "\n")
      }
    }
  }
  
  result_list <- list(
    plot = final_plot,
    value_matrix = value_matrix,
    pvalue_matrix = pvalue_matrix,
    significance_matrix = sig_matrix,
    dataset_annotation = dataset_annotation,
    combined_data = heatmap_data_clean
  )
  
  if (return_plot_data && heatmap_engine == "ggplot2") {
    result_list$plot_data <- plot_data
  }
  
  return(result_list)
}

################################################################################
# Example usage (edit paths; does not run on source)
################################################################################
if (FALSE) {
  result <- plot_custom_heatmap(
    file_list = c(
      "path/to/species_a_go_enrichment.csv",
      "path/to/species_b_go_enrichment.csv"
    ),
    dataset_names = c("Species A", "Species B"),
    dataset_groups = c(
      "Species A" = "Group1",
      "Species B" = "Group2"
    ),
    group_colors = c("Group1" = "#9b383a", "Group2" = "#4878a0"),
    row_column = "Description",
    value_column = "GeneRatio",
    pvalue_column = "p.adjust",
    scale_data = "none",
    min_datasets = 1,
    max_elements = 40,
    selection_criteria = "frequency",
    main_title = "Pathway enrichment across species",
    show_dataset_annotation = TRUE,
    show_significance = TRUE,
    heatmap_engine = "ComplexHeatmap",
    color_palette = "RdBu",
    figure_width = 12,
    figure_height = 10,
    verbose = TRUE
  )
}
