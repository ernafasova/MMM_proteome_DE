
# load metadata
metadata <- read_excel(here::here("data-raw/MMM_metadata.xlsx"))
metadata <- metadata %>% mutate(Group = paste(sex, diet, strain, weeks, sep = "_"))
metadata <- metadata[order(as.numeric(sub("s", "", metadata$sample_id))), ]
metadata$phos_MS_id <- paste0("S", seq(7, 6 + nrow(metadata)))


usethis::use_data(metadata, overwrite = TRUE)

getwd()
colnames(metadata)

library(dplyr)
library(writexl)

# columns to remove
cols_to_drop <- c(
    "hist_id","cage_id","Microvesicular","Macrovesicular","Hypertrophy",
    "SteatosisScore","AvgInfl","InflammationScore","MAS","Ishak","...24",
    "phos_MS_id", "Group"
)

# drop them (silently skips any that don’t exist)
metadata_MMM <- metadata %>% dplyr::select(-any_of(cols_to_drop))


# save
out_dir  <- "C:/Users/pqc291/Documents/GitHub/MMM_proteome_DE/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_xlsx(metadata_MMM, file.path(out_dir, "metadata_MMM.xlsx"))
