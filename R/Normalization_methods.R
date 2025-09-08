# mamual median center (additive, since data is log2)
col_meds <- apply(male_log2_filt, 2, median, na.rm = TRUE)
male_log2_filt_median <- sweep(male_log2_filt, 2, col_meds, FUN = "-")


############################### boxplot ########################################
png(file.path("doc", "Boxplot_male.png"), width = 13, height = 6, units = "in", res = 300, bg = "white")
boxplot(male_log2_filt_median,
        las = 2,                # rotate x-axis labels
        outline = FALSE,        # hide extreme outliers for clarity
        col = "lightgreen",
        main = "After median centering normalization",
        ylab = "Log2 intensity")
dev.off()


############## density plot
# Convert data to long format
df_long <- tidyr::pivot_longer(male_log2_filt_median,
                               cols = everything(),
                               names_to = "new_sample_id",
                               values_to = "intensity")

density_plot_male_norm <- ggplot(df_long, aes(x = intensity, color = new_sample_id)) +
    geom_density() +
    labs(title = "After normalization, male", x = "Intensity", y = "Density") +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
    )

density_plot_male_norm
plotly::ggplotly(density_plot_male_norm)

ggsave(file.path("doc", "density_plot_male_norm.png"), plot = density_plot_male_norm, width = 5, height = 3, dpi = 300, bg = "white")


# center the medians without altering variance. This keeps any real biological variability in spread.
male_filt_medianPhosR_raw <- PhosR::medianScaling(male_raw_filt, scale = FALSE) # data is in original raw intensity scale
male_filt_medianPhosR_raw <- as.data.frame(male_filt_medianPhosR_raw)
male_filt_medianPhosR_raw_log2 <- log2(male_filt_medianPhosR_raw)
all.equal(male_filt_medianPhosR_raw_log2, male_log2_filt_median, tolerance = 1e-10)


male_filt_medianPhosR_log2 <- PhosR::medianScaling(male_log2_filt, scale = FALSE)
male_filt_medianPhosR_log2 <- as.data.frame(male_filt_medianPhosR_log2)



all.equal(male_filt_medianPhosR_raw, male_filt_medianPhosR_log2, tolerance = 1e-10)  # should be TRUE (or extremely close)
all.equal(male_filt_medianPhosR_raw, male_log2_filt_median, tolerance = 1e-10)  # should be TRUE (or extremely close)



png(file.path("doc", "Boxplot_male_norm_PhosR.png"), width = 13, height = 6, units = "in", res = 300, bg = "white")
boxplot(male_filt_medianPhosR_raw,
        las = 2,                # rotate x-axis labels
        outline = FALSE,        # hide extreme outliers for clarity
        col = "orange",
        main = "After median centering PhosR, using raw values",
        ylab = "Log2 intensity")
dev.off()


png(file.path("doc", "Boxplot_male_norm_PhosR_log2.png"), width = 13, height = 6, units = "in", res = 300, bg = "white")
boxplot(male_filt_medianPhosR_log2,
        las = 2,                # rotate x-axis labels
        outline = FALSE,        # hide extreme outliers for clarity
        col = "orange",
        main = "After median centering PhosR, using log2 values",
        ylab = "Log2 intensity")
dev.off()
