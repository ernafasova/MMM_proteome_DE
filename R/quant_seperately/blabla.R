
# Load data and metadata\ -------------------------------------------------

prot_M_renamed_filt

metadata_chow <- metadata_M %>%
    dplyr::filter(
        diet == "LFD"
    )

head(metadata_chow)

data_chow <- prot_M_renamed_filt %>%
    dplyr::select(metadata_chow$sample_id)

head(data_chow)

number_proteins <- data.frame(
    n_proteins = colSums(!is.na(data_chow))) %>%
    tibble::rownames_to_column("sample_id")

data_cvs <- data_chow %>%
    tibble::rownames_to_column("genes") %>%
    tidyr::pivot_longer(
        cols = !genes,
        values_to = "intensities",
        names_to = "sample_id"
    ) %>%
    group_by(genes) %>%
    dplyr::summarise(
        mean = mean(intensities, na.rm = TRUE),
        std = sd(intensities, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(cv = (std/mean)*100)

data_pca_pre_batch <- data_chow %>%
    limma::normalizeBetweenArrays(method = "scale")

tmp <- data_pca_pre_batch %>%
    as.data.frame() %>%
    mutate(
        across(, as.numeric)
    )
    sva::ComBat(batch = metadata_chow$beatbox_batch) %>%
    PhosR::tImpute(m = 1.8, s = 0.3) %>%
    t() %>%
    prcomp()

data_pca_pre_batch <- data_pca_pre_batch$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::inner_join(
        metadata_chow
    ) %>%
    dplyr::inner_join(number_proteins)

data_pca_pre_batch %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = PC1,
            y = PC2,
            color = diet
        )
    ) +
    geom_point()

norm_factor <- data_chow %>%
    tibble::rownames_to_column("genes") %>%
    dplyr::filter(
        genes == "H2ac4;H2ac6;H2ac7;H2ac8;H2ac11;H2ac13;Hist1h2an;Hist1h2ao;Hist1h2ap;H2ac25;Hist1h2af;H2ac12;H2ac15;H2aj"
    ) %>%
    tibble::column_to_rownames("genes")

data_pca_post_batch <- data_chow/as.vector(norm_factor)

data_pca_post_batch <- data_pca_post_batch %>%
    PhosR::tImpute(m =1.8, s = 0.3) %>%
    t() %>%
    prcomp()

data_pca_post_batch <- data_pca_post_batch$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::inner_join(
        metadata_chow
    )

data_pca_post_batch %>%
    ggplot2::ggplot(
        ggplot2::aes(
            x = PC1,
            y = PC2,
            color = weeks
        )
    ) +
    geom_point()
