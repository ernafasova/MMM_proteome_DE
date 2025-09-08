
# load metadata
metadata <- read_excel(here::here("data-raw/MMM_metadata.xlsx"))
metadata <- metadata %>% mutate(Group = paste(sex, diet, strain, weeks, sep = "_"))
metadata <- metadata[order(as.numeric(sub("s", "", metadata$sample_id))), ]
metadata$phos_MS_id <- paste0("S", seq(7, 6 + nrow(metadata)))


usethis::use_data(metadata, overwrite = TRUE)
