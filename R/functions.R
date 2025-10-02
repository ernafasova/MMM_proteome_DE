profile_boxplot <- function(
        protein_id,
        quant_data,
        metadata_df,
        new_id_col = "new_sample_id",
        group_col = "diet",
        title_size = 16,
        axis_title_size = 14,
        axis_text_size = 12,
        custom_colors = NULL
) {

    # --- 1) Prepare the data ---------------------------------------------------
    plot_data <- quant_data %>%
        tibble::rownames_to_column("protein") %>%
        dplyr::filter(protein == protein_id) %>%
        tidyr::pivot_longer(
            cols = -protein,
            names_to = new_id_col,
            values_to = "log2_intensity"
        ) %>%
        dplyr::left_join(metadata_df, by = new_id_col)

    if (nrow(plot_data) == 0) {
        stop(paste("Protein", protein_id, "not found in the dataset."))
    }

    # --- 2) Set colors if not provided -----------------------------------------
    if (is.null(custom_colors)) {
        groups <- unique(plot_data[[group_col]])
        custom_colors <- scales::hue_pal()(length(groups))
        names(custom_colors) <- groups
    }

    # --- 3) Plot ---------------------------------------------------------------
    p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes_string(x = group_col, y = "log2_intensity", fill = group_col)
    ) +
        ggplot2::geom_boxplot(
            width = 0.6, alpha = 0.6,
            outlier.shape = NA, color = "black"
        ) +
        ggplot2::geom_jitter(
            width = 0.15, size = 2.5, alpha = 0.8,
            color = "black"
        ) +
        ggplot2::scale_fill_manual(values = custom_colors) +
        ggplot2::labs(
            title = protein_id,
            x = group_col,
            y = "Log2 Intensity"
        ) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
            axis.title = ggplot2::element_text(size = axis_title_size),
            axis.text = ggplot2::element_text(size = axis_text_size, color = "black"),
            legend.position = "none",
            panel.border = ggplot2::element_rect(color = "black", size = 0.8),
            panel.grid.major = ggplot2::element_line(color = "grey85", size = 0.3),
            panel.grid.minor = ggplot2::element_blank()
        )

    return(p)
}



################################################################################
plot_pca <- function(
        data,
        color_var,
        shape_var    = NULL,
        color_vals   = NULL,   # manual colors
        shape_vals   = NULL,   # manual shapes
        palette_name = "Set1",
        plot_title   = NULL,
        ellipse      = TRUE
) {
    # Treat aesthetics as factors for discrete scales
    data[[color_var]] <- as.factor(data[[color_var]])
    if (!is.null(shape_var)) {
        data[[shape_var]] <- as.factor(data[[shape_var]])
    }

    # Aesthetic mappings
    mapping <- aes(
        x = PC1,
        y = PC2,
        color = .data[[color_var]],
        text = hover_text
    )
    if (!is.null(shape_var)) {
        mapping <- aes(
            x = PC1,
            y = PC2,
            color = .data[[color_var]],
            shape = .data[[shape_var]],
            text = hover_text
        )
    }

    # Base plot
    p <- ggplot(data, mapping) +
        geom_point(size = 3.5, alpha = 0.7) +
        labs(
            title    = plot_title,
            subtitle = paste("based on", n_nona, "sites out of", n_original),
            x        = paste0("PC1 (", pca_var_perc[1], "%)"),
            y        = paste0("PC2 (", pca_var_perc[2], "%)")
        ) +
        theme_bw(base_size = 14) +
        theme(
            plot.title    = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5)
        )

    # Optional ellipses
    if (ellipse) {
        p <- p +
            stat_ellipse(
                aes(
                    group = .data[[color_var]],
                    fill  = .data[[color_var]]
                ),
                geom        = "path",
                show.legend = FALSE
            )
    }

    # Color scale
    if (!is.null(color_vals)) {
        p <- p + scale_color_manual(values = color_vals)
    } else {
        p <- p + scale_color_brewer(palette = palette_name)
    }

    # Shape scale (manual override)
    if (!is.null(shape_vals)) {
        p <- p + scale_shape_manual(values = shape_vals)
    }

    print(p)
}


###############################################################################
get_limma_results <- function(fit_obj, coef_name, contrast_label, alpha = 0.05) {

    # Extract results table
    top_table <- limma::topTable(
        fit_obj,
        coef          = coef_name,
        adjust.method = "BH",
        sort.by       = "logFC",
        number        = Inf
    ) %>%
        rownames_to_column("Protein")

    # Filter significant hits
    sig <- top_table %>%
        filter(adj.P.Val < alpha) %>%
        drop_na()

    up   <- sig %>% filter(logFC > 0)
    down <- sig %>% filter(logFC < 0)

    # Console summary
    message("\n---- ", contrast_label, " ----")
    message("Total significant: ", nrow(sig))
    message("UP: ", nrow(up), " | DOWN: ", nrow(down))

    # Return result lists
    return(list(
        top_table = top_table,
        sig       = sig,
        up        = up,
        down      = down
    ))
}


############################ volcano
create_volcano_plot <- function(
        df,
        label_proteins = NULL,
        top_n = 10,
        title = "Volcano Plot",
        x_label = "Log2 Fold Change",
        y_label = "-Log10(P-value)"
) {
    # Ensure Protein column exists, create it from rownames if it doesn't
    if (!("Protein" %in% colnames(df))) {
        df <- tibble::rownames_to_column(df, "Protein")
    }

    # Assign significance group and create tooltip
    df <- df %>%
        dplyr::mutate(
            SigGroup = dplyr::case_when(
                adj.P.Val < 0.05 & logFC > 0 ~ "Up",
                adj.P.Val < 0.05 & logFC < 0 ~ "Down",
                TRUE ~ "NotSig"
            ),
            tooltip = paste(
                "Protein:", Protein,
                "<br>logFC:", signif(logFC, 3),
                "<br>P.Value:", signif(P.Value, 3),
                "<br>Adj.P.Val:", signif(adj.P.Val, 3)
            )
        )

    # Calculate significant counts
    up_count <- sum(df$SigGroup == "Up", na.rm = TRUE)
    down_count <- sum(df$SigGroup == "Down", na.rm = TRUE)

    # Select proteins to label
    if (!is.null(label_proteins)) {
        label_df <- dplyr::filter(df, Protein %in% label_proteins)
    } else {
        top_up <- df %>%
            dplyr::filter(SigGroup == "Up") %>%
            dplyr::arrange(P.Value) %>%
            head(top_n)
        top_down <- df %>%
            dplyr::filter(SigGroup == "Down") %>%
            dplyr::arrange(P.Value) %>%
            head(top_n)
        label_df <- dplyr::bind_rows(top_up, top_down)
    }

    xmax <- max(abs(df$logFC), na.rm = TRUE) * 1.1
    ymax <- max(-log10(df$P.Value), na.rm = TRUE) * 1.1

    p <- ggplot2::ggplot(df, aes(x = logFC, y = -log10(P.Value), color = SigGroup, text = tooltip)) +
        geom_point(alpha = 0.4, size = 2, shape = 1) +
        scale_color_manual(values = c("Up" = "#DC4131", "Down" = "#4390DE", "NotSig" = "grey70")) +
        ggrepel::geom_label_repel(
            data = label_df,
            ggplot2::aes(label = Protein, x = logFC, y = -log10(P.Value)), # Changed label to Protein
            fill = "white",
            color = "black",
            box.padding = 0.5,
            segment.color = "grey50",
            size = 3.5,
            show.legend = FALSE,
            inherit.aes = FALSE,
            max.overlaps = Inf
        ) +
        labs(title = title, x = x_label, y = y_label) +
        theme_minimal(base_size = 15) +
        theme(
            legend.position = "none",
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 1)
        ) +
        annotate(
            "text",
            x = xmax * 0.95,
            y = ymax * 0.05,
            label = paste("Up:", up_count),
            color = "#DC4131",
            size = 5,
            hjust = 1
        ) +
        annotate(
            "text",
            x = -xmax * 0.95,
            y = ymax * 0.05,
            label = paste("Down:", down_count),
            color = "#4390DE",
            size = 5,
            hjust = 0
        ) +
        ggplot2::xlim(-xmax, xmax)

    return(p)
}

#############3
create_volcano_plot_FC <- function(
        df,
        label_proteins = NULL,
        top_n = 10,
        logfc_threshold = 1, # <-- New argument for log2 fold change threshold
        p_value_threshold = 0.05,
        title = "Volcano Plot",
        x_label = "Log2 Fold Change",
        y_label = "-Log10(P-value)"
) {

    # Ensure Protein column exists, create it from rownames if it doesn't
    if (!("Protein" %in% colnames(df))) {
        df <- tibble::rownames_to_column(df, "Protein")
    }

    # Assign significance group using BOTH p-value and logFC thresholds
    df <- df %>%
        dplyr::mutate(
            SigGroup = dplyr::case_when(
                adj.P.Val < p_value_threshold & logFC > logfc_threshold  ~ "Up",
                adj.P.Val < p_value_threshold & logFC < -logfc_threshold ~ "Down",
                TRUE                                                    ~ "NotSig"
            ),
            tooltip = paste(
                "Protein:", Protein,
                "<br>logFC:", signif(logFC, 3),
                "<br>P.Value:", signif(P.Value, 3),
                "<br>Adj.P.Val:", signif(adj.P.Val, 3)
            )
        )

    # Calculate significant counts (this will now be lower)
    up_count <- sum(df$SigGroup == "Up", na.rm = TRUE)
    down_count <- sum(df$SigGroup == "Down", na.rm = TRUE)

    # Select proteins to label
    if (!is.null(label_proteins)) {
        label_df <- dplyr::filter(df, Protein %in% label_proteins)
    } else {
        top_up <- df %>% dplyr::filter(SigGroup == "Up") %>% dplyr::arrange(P.Value) %>% head(top_n)
        top_down <- df %>% dplyr::filter(SigGroup == "Down") %>% dplyr::arrange(P.Value) %>% head(top_n)
        label_df <- dplyr::bind_rows(top_up, top_down)
    }

    # Dynamically set plot limits
    xmax <- max(abs(df$logFC), na.rm = TRUE) * 1.1
    ymax <- max(-log10(df$P.Value), na.rm = TRUE) * 1.1

    p <- ggplot2::ggplot(df, aes(x = logFC, y = -log10(P.Value), color = SigGroup, text = tooltip)) +
        # Add lines for thresholds
        geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", color = "grey50") +

        geom_point(alpha = 0.4, size = 2, shape = 1) +
        scale_color_manual(values = c("Up" = "#DC4131", "Down" = "#4390DE", "NotSig" = "grey70")) +
        ggrepel::geom_label_repel(
            data = label_df,
            ggplot2::aes(label = Protein, x = logFC, y = -log10(P.Value)),
            fill = "white", color = "black",
            box.padding = 0.5, segment.color = "grey50",
            size = 3.5, show.legend = FALSE,
            inherit.aes = FALSE, max.overlaps = Inf
        ) +
        labs(title = title, x = x_label, y = y_label) +
        theme_minimal(base_size = 15) +
        theme(
            legend.position = "none",
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 1)
        ) +
        annotate(
            "text", x = xmax * 0.95, y = ymax * 0.05, label = paste("Up:", up_count),
            color = "#DC4131", size = 5, hjust = 1
        ) +
        annotate(
            "text", x = -xmax * 0.95, y = ymax * 0.05, label = paste("Down:", down_count),
            color = "#4390DE", size = 5, hjust = 0
        ) +
        ggplot2::xlim(-xmax, xmax)

    return(p)
}
