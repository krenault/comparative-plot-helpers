library(ape); library(ggtree); library(ggtreeExtra); library(ggplot2); library(dplyr); library(ggnewscale)

plot_phylogenetic_tree <- function(tree_file, metadata = NULL, color_file = NULL,
                                   bar_column = NULL, bar_title = NULL, bar_legend = NULL, color_column = NULL,
                                   species_column = NULL, assembly_column = NULL,
                                   output_prefix = NULL, output_dir = NULL,
                                   tree_size = NULL, tip_shape = NULL, tip_legend = NULL, tip_size = NULL, bar_width = NULL,
                                   figure_width = NULL, figure_height = NULL,
                                   layout_preference = NULL, show_color_legend = FALSE, show_shape_legend = TRUE, show_bar_legend = TRUE,
                                   make_ultrametric = FALSE) {
  
  cat("Loading tree from:", tree_file, "\n")
  tree <- read.tree(tree_file)
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  
  # ENHANCED: Make tree ultrametric to align all tips
  if (make_ultrametric) {
    cat("Converting tree to ultrametric (all tips aligned)...\n")
    # Calculate depth of each tip from root
    tip_depths <- node.depth.edgelength(tree)[1:length(tree$tip.label)]
    max_depth <- max(tip_depths)
    
    # Extend terminal branches to reach max depth
    # Find which edges lead to tips
    tip_edges <- which(tree$edge[,2] <= length(tree$tip.label))
    
    for (i in tip_edges) {
      tip_node <- tree$edge[i, 2]
      current_depth <- tip_depths[tip_node]
      extension_needed <- max_depth - current_depth
      tree$edge.length[i] <- tree$edge.length[i] + extension_needed
    }
    
    cat("Extended terminal branches to align all tips\n")
  }
  
  get_jewel_palette <- function(n) {
    jewel_colors <- c(
      "#2c6e5f", "#9b383a", "#4878a0", "#64a590", "#945a87",
      "#b68ab9", "#9c6e5a", "#9fbd8b", "#a85472", "#86acb9",
      "#20854c", "#9c534a", "#2c6e5f", "#31a6ad", "#f0a0a3",
      "#c57251", "#856890", "#7e9960", "#7d9a9e", "#B65B32",
      "#c07a7a", "#8a7c6d", "#AA4839", "#4878a0", "#3D7B68"
    )
    if (n <= length(jewel_colors)) return(jewel_colors[1:n]) else return(rainbow(n))
  }
  
  color_palette <- NULL
  if (!is.null(color_file) && file.exists(color_file)) {
    cat("Loading color mapping from:", color_file, "\n")
    color_data <- read.table(color_file, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(color_data) >= 2) {
      color_palette <- setNames(color_data[,2], color_data[,1])
      cat("Loaded", length(color_palette), "category-color mappings\n")
    } else if (ncol(color_data) == 1) {
      color_palette <- color_data[,1]
      cat("Loaded", length(color_palette), "hex codes for automatic assignment\n")
    } else {
      cat("Warning: Color file should have 1 column (hex codes) or 2 columns (category, hex_code)\n")
    }
  }
  
  plot_data <- NULL; tip_colors <- NULL; assigned_colors <- NULL
  
  if (!is.null(metadata)) {
    cat("Processing metadata...\n")
    if (is.character(metadata) && file.exists(metadata)) {
      metadata <- read.csv(metadata, stringsAsFactors = FALSE)
    }
    
    if (!is.null(species_column) && species_column %in% colnames(metadata)) {
      # Convert underscores to spaces in species names to match tree tip labels
      metadata[[species_column]] <- gsub("_", " ", metadata[[species_column]])
      metadata <- metadata[!duplicated(metadata[[species_column]]), ]
      cat("Removed", sum(duplicated(metadata[[species_column]])), "duplicate species entries\n")
      exact_matches <- intersect(tree$tip.label, metadata[[species_column]])
    } else if (!is.null(assembly_column) && assembly_column %in% colnames(metadata)) {
      # Convert underscores to spaces in assembly names to match tree tip labels
      metadata[[assembly_column]] <- gsub("_", " ", metadata[[assembly_column]])
      metadata <- metadata[!duplicated(metadata[[assembly_column]]), ]
      cat("Removed", sum(duplicated(metadata[[assembly_column]])), "duplicate assembly entries\n")
      exact_matches <- intersect(tree$tip.label, metadata[[assembly_column]])
      species_column <- assembly_column
    } else {
      cat("Warning: Specified species_column or assembly_column not found in metadata\n")
      exact_matches <- character(0)
    }
    
    cat("Found", length(exact_matches), "matching species\n")
    if (length(exact_matches) > 0) {
      plot_data <- metadata %>% filter(.data[[species_column]] %in% exact_matches)
      tree <- drop.tip(tree, setdiff(tree$tip.label, exact_matches))
      
      if (!is.null(color_column) && color_column %in% colnames(plot_data)) {
        cat("Setting up colors based on:", color_column, "\n")
        plot_data[[color_column]] <- trimws(as.character(plot_data[[color_column]]))
        unique_values <- unique(plot_data[[color_column]])
        unique_values <- unique_values[!is.na(unique_values)]
        
        # Check if color_column contains hex codes directly
        is_hex <- all(grepl("^#[0-9A-Fa-f]{6}$", unique_values[!is.na(unique_values)]))
        
        if (is_hex) {
          # Direct hex codes: use them as-is
          cat("Detected hex color codes in column, using them directly\n")
          assigned_colors <- setNames(unique_values, unique_values)
        } else if (!is.null(color_palette)) {
          # Category names with color palette
          if (is.named(color_palette)) {
            assigned_colors <- color_palette
            missing_values <- setdiff(unique_values, names(color_palette))
            if (length(missing_values) > 0) {
              additional_colors <- get_jewel_palette(length(missing_values))
              names(additional_colors) <- missing_values
              assigned_colors <- c(assigned_colors, additional_colors)
            }
          } else {
            if (length(color_palette) >= length(unique_values)) {
              assigned_colors <- setNames(color_palette[1:length(unique_values)], unique_values)
            } else {
              assigned_colors <- setNames(color_palette, unique_values[1:length(color_palette)])
              remaining_values <- unique_values[(length(color_palette) + 1):length(unique_values)]
              if (length(remaining_values) > 0) {
                additional_colors <- get_jewel_palette(length(remaining_values))
                names(additional_colors) <- remaining_values
                assigned_colors <- c(assigned_colors, additional_colors)
              }
            }
            cat("Assigned", length(color_palette), "provided colors to categories\n")
          }
        } else {
          # Category names without color palette: assign jewel colors
          assigned_colors <- get_jewel_palette(length(unique_values))
          names(assigned_colors) <- unique_values
        }
        if (any(is.na(plot_data[[color_column]]))) {
          assigned_colors <- c(assigned_colors, "NA" = "#808080")
          plot_data[[color_column]][is.na(plot_data[[color_column]])] <- "NA"
        }
        tip_colors <- setNames(plot_data[[color_column]], plot_data[[species_column]])
      }
    }
  } else {
    cat("No metadata provided, creating basic tree\n")
  }
  
  cat("Creating tree plot...\n")
  if (!is.null(layout_preference)) {
    p <- ggtree(tree, layout = layout_preference)
  } else {
    p <- ggtree(tree)
  }
  
  if (!is.null(tip_colors)) {
    tip_data <- data.frame(
      label = tree$tip.label,
      color_value = tip_colors[tree$tip.label],
      stringsAsFactors = FALSE
    )
    tip_data$color_value[is.na(tip_data$color_value)] <- "NA"
    
    # DEBUG: Print what we're actually using
    cat("DEBUG: First 10 color values:\n")
    print(head(tip_data$color_value, 10))
    cat("DEBUG: assigned_colors is NULL?", is.null(assigned_colors), "\n")
    cat("DEBUG: Are values hex codes?", all(grepl("^#[0-9A-Fa-f]{6}$", tip_data$color_value[!is.na(tip_data$color_value) & tip_data$color_value != "NA"])), "\n")
    
    if (!is.null(tip_shape) && !is.null(plot_data) && tip_shape %in% colnames(plot_data)) {
      shape_mapping_data <- setNames(plot_data[[tip_shape]], plot_data[[species_column]])
      tip_data$shape_value <- shape_mapping_data[tip_data$label]
    }
    
    if (!is.null(tree_size)) {
      p <- p %<+% tip_data +
        geom_tree(aes(color = color_value), size = tree_size)
    } else {
      p <- p %<+% tip_data +
        geom_tree(aes(color = color_value))
    }
    
    # Use identity scale if assigned_colors is NULL (meaning hex codes) or if color_value contains hex codes
    if (is.null(assigned_colors) || all(grepl("^#[0-9A-Fa-f]{6}$", tip_data$color_value[!is.na(tip_data$color_value) & tip_data$color_value != "NA"]))) {
      cat("Using scale_color_identity\n")
      p <- p + scale_color_identity(guide = if(show_color_legend) "legend" else "none")
    } else {
      cat("Using scale_color_manual\n")
      p <- p + scale_color_manual(name = color_column, values = assigned_colors, na.value = "#504f4f", guide = if(show_color_legend) "legend" else "none")
    }
    
    if (!is.null(tip_size)) {
      p <- p + geom_tiplab(aes(color = color_value), size = tip_size, offset = 1, hjust = 0)
    } else {
      p <- p + geom_tiplab(aes(color = color_value), offset = 0.05, hjust = 0)
    }
  } else {
    if (!is.null(tree_size)) {
      p <- p + geom_tree(size = tree_size)
    } else {
      p <- p + geom_tree()
    }
    if (!is.null(tip_size)) {
      p <- p + geom_tiplab(size = tip_size, offset = 0.05, hjust = 0)
    } else {
      p <- p + geom_tiplab(offset = 0.05, hjust = 0)
    }
  }
  
  tree_depth <- max(node.depth.edgelength(tree))
  if (!is.null(bar_column) && !is.null(plot_data) && bar_column %in% colnames(plot_data)) {
    p <- p + xlim(NA, tree_depth * 2.7)
  } else {
    p <- p + xlim(NA, tree_depth * 1.8)
  }
  
  if (!is.null(bar_column) && !is.null(plot_data) && bar_column %in% colnames(plot_data)) {
    cat("Adding bar plot for:", bar_column, "\n")
    bar_data <- plot_data %>% 
      filter(!is.na(.data[[bar_column]]) & is.finite(.data[[bar_column]]))
    
    if (nrow(bar_data) > 0) {
      if (!is.null(color_column) && !is.null(assigned_colors)) {
        if (!is.null(bar_width)) {
          p <- p + geom_fruit(data = bar_data,
                              geom = geom_col,
                              mapping = aes(x = !!sym(bar_column),
                                            y = !!sym(species_column),
                                            fill = !!sym(color_column)),
                              orientation = "y",
                              width = bar_width,
                              color = "#504f4f",
                              size = 0.3,
                              offset = 1,
                              axis.params = list(
                                axis = "x",
                                text.size = 2.5,
                                hjust = 0.5,
                                vjust = 1,
                                text.angle = -45,
                                nbreak = 3,
                                line.size = 0.5,
                                line.color = "black"
                              ))
        } else {
          p <- p + geom_fruit(data = bar_data,
                              geom = geom_col,
                              mapping = aes(x = !!sym(bar_column),
                                            y = !!sym(species_column),
                                            fill = !!sym(color_column)),
                              orientation = "y",
                              color = "#504f4f",
                              size = 0.3,
                              offset = 0.05,
                              axis.params = list(
                                axis = "x",
                                text.size = 1.5,
                                hjust = 0.5,
                                vjust = 1,
                                text.angle = -45,
                                nbreak = 3,
                                line.size = 0.5,
                                line.color = "black"
                              ))
        }
        p <- p + scale_fill_manual(name = ifelse(is.null(bar_legend), bar_title, bar_legend), values = assigned_colors, guide = if(show_bar_legend) "legend" else "none")
        
        p <- p + annotate("text",
                          x = tree_depth * 2.1,
                          y = length(tree$tip.label) + 1,
                          label = paste(bar_title, "(years)"),
                          size = 4,
                          fontface = "bold",
                          hjust = 0.5)
      } else {
        if (!is.null(bar_width)) {
          p <- p + geom_fruit(data = bar_data,
                              geom = geom_col,
                              mapping = aes(x = !!sym(bar_column),
                                            y = !!sym(species_column)),
                              fill = "grey",
                              orientation = "y",
                              width = bar_width,
                              color = "#504f4f",
                              size = 0.01,
                              offset = 0.05,
                              axis.params = list(
                                axis = "x",
                                text.size = 2.5,
                                hjust = 0.5,
                                vjust = 1,
                                text.angle = -45,
                                nbreak = 3,
                                line.size = 0.5,
                                line.color = "black"
                              ))
        } else {
          p <- p + geom_fruit(data = bar_data,
                              geom = geom_col,
                              mapping = aes(x = !!sym(bar_column),
                                            y = !!sym(species_column)),
                              fill = "grey",
                              orientation = "y",
                              color = "#504f4f",
                              size = 0.3,
                              offset = 0.05,
                              axis.params = list(
                                axis = "x",
                                text.size = 2.5,
                                hjust = 0.5,
                                vjust = 1,
                                text.angle = -45,
                                nbreak = 6,
                                line.size = 0.5,
                                line.color = "black"
                              ))
        }
        p <- p + annotate("text",
                          x = tree_depth * 2.1,
                          y = length(tree$tip.label) + 1,
                          label = paste(bar_column, "(years)"),
                          size = 4,
                          fontface = "bold",
                          hjust = 0.5)
      }
    }
  }
  
  if (!is.null(tip_shape) && !is.null(plot_data) && tip_shape %in% colnames(plot_data)) {
    cat("Adding shapes based on:", tip_shape, "\n")
    unique_shapes <- unique(plot_data[[tip_shape]])
    unique_shapes <- unique_shapes[!is.na(unique_shapes)]
    available_shapes <- c(16, 17, 18, 15, 19, 8, 7, 9, 10, 11, 12, 13, 14)
    if (length(unique_shapes) <= length(available_shapes)) {
      shape_mapping <- setNames(available_shapes[1:length(unique_shapes)], unique_shapes)
    } else {
      shape_mapping <- setNames(rep(available_shapes, length.out = length(unique_shapes)), unique_shapes)
    }
    
    if (!is.null(tip_colors) && !is.null(assigned_colors)) {
      p <- p + geom_tippoint(aes(shape = shape_value, color = color_value),
                             size = 2.5) +
        scale_shape_manual(values = shape_mapping,
                           name = ifelse(is.null(tip_legend), tip_shape, tip_legend),
                           guide = if(show_shape_legend) "legend" else "none")
    } else {
      p <- p + geom_tippoint(aes(shape = shape_value),
                             size = 2.5) +
        scale_shape_manual(values = shape_mapping,
                           name = ifelse(is.null(tip_legend), tip_shape, tip_legend),
                           guide = if(show_shape_legend) "legend" else "none")
    }
  }
  
  p <- p +
    theme(legend.position = "right",
          legend.box = "vertical",
          plot.margin = margin(10, 10, 10, 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))
  
  print(p)
  
  if (!is.null(output_dir) && !is.null(output_prefix)) {
    cat("Saving outputs...\n")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    tree_output <- file.path(output_dir, paste0(output_prefix, "_pruned.nwk"))
    write.tree(tree, file = tree_output)
    pdf_output <- file.path(output_dir, paste0(output_prefix, "_plot.pdf"))
    png_output <- file.path(output_dir, paste0(output_prefix, "_plot.png"))
    if (!is.null(figure_width) && !is.null(figure_height)) {
      ggsave(pdf_output, p, width = figure_width, height = figure_height, dpi = 300)
      ggsave(png_output, p, width = figure_width, height = figure_height, dpi = 300, bg = "white")
    } else {
      ggsave(pdf_output, p, dpi = 300)
      ggsave(png_output, p, dpi = 300, bg = "white")
    }
    cat("Files saved:\n - Tree:", tree_output, "\n - PDF:", pdf_output, "\n - PNG:", png_output, "\n")
  }
  
  cat("\n=== PHYLOGENETIC TREE SUMMARY ===\n")
  cat("Tree tips:", length(tree$tip.label), "\n")
  if (!is.null(plot_data)) {
    cat("Species with metadata:", nrow(plot_data), "\n")
    if (!is.null(color_column)) {
      cat("Color groups:", paste(unique(plot_data[[color_column]]), collapse = ", "), "\n")
    }
    if (!is.null(bar_column)) {
      bar_stats <- plot_data[[bar_column]][!is.na(plot_data[[bar_column]])]
      if (length(bar_stats) > 0) {
        cat("Bar plot range:", round(min(bar_stats), 2), "-", round(max(bar_stats), 2), "\n")
        cat("Bar plot mean:", round(mean(bar_stats), 2), "\n")
      }
    }
  }
  
  return(list(plot = p,
              tree = tree,
              data = plot_data,
              colors = assigned_colors))
}

################################################################################
# Example usage (edit paths; does not run on source)
################################################################################
if (FALSE) {
  result <- plot_phylogenetic_tree(
    tree_file = "path/to/tree.nwk",
    metadata = "path/to/metadata.csv",
    color_column = "Lifestyle",
    species_column = "Species",
    bar_column = "maximum_longevity_y",
    bar_title = "Maximum longevity",
    bar_legend = "Lifestyle",
    tip_legend = "Available data",
    tip_shape = "Data_Type",
    tip_size = 5,
    bar_width = 0.8,
    figure_width = 20,
    figure_height = 20,
    show_color_legend = FALSE,
    show_shape_legend = TRUE,
    show_bar_legend = TRUE
  )
  ggsave("path/to/output.png", result$plot, width = 12, height = 12, dpi = 300)
}
