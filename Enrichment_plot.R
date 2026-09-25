plot_go_enrichment <- function(
    data_list,
    dataset_names = NULL,
    output_path = NULL,
    pathway_col = NULL,
    fold_enrichment_col = NULL, 
    fdr_col = NULL,
    keywords = NULL,
    n_pathways = NULL,
    selection_method = NULL,
    min_fold_diff = NULL,
    significance_threshold = NULL,
    colors = NULL,
    point_size_range = NULL,
    plot_width = NULL,
    plot_height = NULL,
    plot_dpi = NULL,
    title = NULL,
    subtitle = NULL,
    x_label = NULL,
    y_label = NULL,
    show_significance_background = NULL,
    dodge_width = NULL,
    theme_base = NULL,
    font_size_title = NULL,
    font_size_axis = NULL,
    font_size_axis_title = NULL,
    legend_position = NULL,
    remove_underscores = NULL,
    verbose = NULL
) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(stringr)

  load_data <- function(data_input, name) {
    if (is.character(data_input) && length(data_input) == 1) {
      if (!file.exists(data_input)) stop(paste("File not found:", data_input))
      data <- read.csv(data_input, stringsAsFactors = FALSE)
    } else if (is.data.frame(data_input)) {
      data <- data_input
    } else {
      stop("Data input must be a file path (character) or data frame")
    }
    
    required_cols <- c(pathway_col, fold_enrichment_col, fdr_col)
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns in", name, ":", paste(missing_cols, collapse = ", ")))
    }
    return(data)
  }
  
  if (!is.list(data_list) || is.data.frame(data_list)) {
    data_list <- list(data_list)
    single_dataset <- TRUE
  } else {
    single_dataset <- FALSE
  }
  
  if (is.null(dataset_names)) {
    if (!is.null(names(data_list))) {
      dataset_names <- names(data_list)
    } else if (single_dataset) {
      dataset_names <- "Dataset"
    } else {
      dataset_names <- paste("Dataset", seq_along(data_list))
    }
  } else if (length(dataset_names) != length(data_list)) {
    stop("Length of dataset_names must match length of data_list")
  }
  
  combined_data_list <- list()
  for (i in seq_along(data_list)) {
    data <- load_data(data_list[[i]], dataset_names[i])
    processed_data <- data %>%
      mutate(
        pathway = if(remove_underscores) gsub("_", " ", .data[[pathway_col]]) else .data[[pathway_col]],
        dataset = dataset_names[i],
        log_p_adj = -log10(.data[[fdr_col]]),
        fold_enrichment = .data[[fold_enrichment_col]],
        fdr = .data[[fdr_col]]
      )
    combined_data_list[[i]] <- processed_data
  }
  
  combined_data <- bind_rows(combined_data_list)
  
  if (verbose) {
    cat("Datasets:", length(data_list), "| Pathways:", nrow(combined_data), "| Unique:", length(unique(combined_data$pathway)), "\n")
  }
  
  if (!is.null(keywords) && length(keywords) > 0) {
    pattern <- paste(keywords, collapse = "|")
    combined_data <- combined_data %>%
      filter(str_detect(tolower(pathway), tolower(pattern)))
    if (nrow(combined_data) == 0) {
      warning("No pathways match keywords. Using all pathways.")
      combined_data <- bind_rows(combined_data_list)
    }
  }
  
  pathway_stats <- combined_data %>%
    group_by(pathway) %>%
    summarise(
      n_datasets = n(),
      enrichment_values = list(fold_enrichment),
      fdr_values = list(fdr),
      avg_enrichment = mean(fold_enrichment, na.rm = TRUE),
      min_fdr = min(fdr, na.rm = TRUE),
      max_enrichment = max(fold_enrichment, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      max_fold_diff = if(length(data_list) > 1) {
        sapply(enrichment_values, function(x) if(length(x) > 1) max(x) - min(x) else 0)
      } else {
        rep(0, n())
      },
      n_significant = sapply(fdr_values, function(x) sum(x < significance_threshold, na.rm = TRUE)),
      is_significant = n_significant > 0,
      priority_score = case_when(
        selection_method == "fold_diff" ~ max_fold_diff + ifelse(is_significant, 10, 0),
        selection_method == "top_fdr" ~ -min_fdr,
        selection_method == "combined" ~ max_fold_diff + (-log10(min_fdr)) + ifelse(is_significant, 5, 0),
        TRUE ~ avg_enrichment
      )
    ) %>%
    arrange(desc(priority_score))
  
  if (selection_method == "top_fdr" && single_dataset) {
    selected_pathways <- combined_data %>% arrange(fdr) %>% head(n_pathways) %>% pull(pathway)
  } else {
    selected_pathways <- pathway_stats %>% head(n_pathways) %>% pull(pathway)
  }
  
  if (verbose) {
    cat("Selected:", length(selected_pathways), "pathways using", selection_method, "\n")
  }
  
  plot_data <- combined_data %>%
    filter(pathway %in% selected_pathways) %>%
    left_join(pathway_stats %>% dplyr::select(pathway, max_fold_diff, is_significant, priority_score), by = "pathway") %>%
    arrange(desc(priority_score)) %>%
    mutate(pathway = factor(pathway, levels = unique(pathway)))
  
  if (is.null(colors)) {
    if (length(dataset_names) <= 3) {
      default_colors <- c("#4878a0", "#945a87", "#5c4d8b")
      colors <- setNames(default_colors[1:length(dataset_names)], dataset_names)
    } else {
      colors <- rainbow(length(dataset_names))
      names(colors) <- dataset_names
    }
  } else if (is.null(names(colors))) {
    names(colors) <- dataset_names[1:length(colors)]
  }
  
  if (is.null(title)) {
    title <- if (single_dataset) paste("GO Enrichment:", dataset_names[1]) else "Comparative GO Enrichment"
  }
  
  p <- ggplot(plot_data, aes(x = fold_enrichment, y = pathway))
  
  if (show_significance_background && !single_dataset) {
    sig_data <- plot_data %>%
      group_by(pathway, is_significant, max_fold_diff) %>%
      summarise(.groups = "drop") %>%
      filter(is_significant | max_fold_diff >= min_fold_diff)
    
    if (nrow(sig_data) > 0) {
      p <- p + geom_rect(data = sig_data, aes(xmin = -Inf, xmax = Inf, ymin = as.numeric(pathway) - 0.4, ymax = as.numeric(pathway) + 0.4), fill = "gold", alpha = 0.2, inherit.aes = FALSE)
    }
  }
  
  p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", alpha = 0.5)
  
  # Determine if x and size should be independent
  same_col <- fold_enrichment_col == fdr_col
  
  if (single_dataset) {
    if (same_col) {
      p <- p + geom_point(aes(x = fold_enrichment), size = mean(point_size_range), color = colors[1], alpha = 0.8)
    } else {
      p <- p + geom_point(aes(size = log_p_adj, x = fold_enrichment), color = colors[1], alpha = 0.8)
    }
  } else {
    if (same_col) {
      p <- p + geom_point(aes(x = fold_enrichment, color = dataset), size = mean(point_size_range),
                          alpha = 0.8, position = position_dodge(width = dodge_width)) +
        scale_color_manual(values = colors, name = "Dataset")
    } else {
      p <- p + geom_point(aes(size = log_p_adj, x = fold_enrichment, color = dataset), alpha = 0.8,
                          position = position_dodge(width = dodge_width)) +
        scale_color_manual(values = colors, name = "Dataset")
    }
  }
  
  
  p <- p + 
    scale_size_continuous(range = point_size_range, name = "-log10(adj.P)") +
    labs(title = title, subtitle = subtitle, x = x_label, y = y_label) +
    theme_base() +
    theme(
      plot.title = element_text(size = font_size_title, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = font_size_title - 2, hjust = 0.5, color = "gray60"),
      axis.text.y = element_text(size = font_size_axis),
      axis.text.x = element_text(size = font_size_axis),
      axis.title = element_text(size = font_size_axis_title),
      legend.position = legend_position,
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  if (!is.null(output_path)) {
    ggsave(output_path, p, width = plot_width, height = plot_height, dpi = plot_dpi)
    if (verbose) cat("Saved to:", output_path, "\n")
  }
  
  return(p)
}

################################################################################
# Example usage (edit paths; does not run on source)
################################################################################
if (FALSE) {
  plot_go_enrichment(
    data_list = "path/to/enrichment.csv",
    pathway_col = "Description",
    fold_enrichment_col = "FoldEnrichment",
    fdr_col = "p.adjust",
    n_pathways = 15,
    selection_method = "top_fdr",
    colors = c("Dataset" = "#945a87"),
    theme_base = theme_minimal,
    title = "GO enrichment",
    remove_underscores = TRUE,
    verbose = TRUE
  )

  # Multi-dataset example:
  # plot_go_enrichment(
  #   data_list = list("path/to/set_a.csv", "path/to/set_b.csv"),
  #   dataset_names = c("Set A", "Set B"),
  #   pathway_col = "Description",
  #   fold_enrichment_col = "pvalue",
  #   fdr_col = "pvalue",
  #   n_pathways = 25,
  #   selection_method = "top_fdr",
  #   colors = c("Set A" = "#5c4d8b", "Set B" = "#4878a0"),
  #   theme_base = theme_minimal
  # )
}
