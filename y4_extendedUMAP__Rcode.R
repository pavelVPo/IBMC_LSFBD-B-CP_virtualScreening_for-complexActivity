library(tidyverse)
library(umap)
library(mirt)
library(praznik)

### Import
# Selected cmpnds
data_main   <- read_tsv(".../13-02-26__selected_cmpnds_described.tab") |> arrange(compound_id)
# Descriptors
desc_cntrl <- read_tsv(".../Smo_controls_SD.tab") |> select(id, mna)
desc_slctd <- read_tsv(".../Smo_selected_SD.tab")  |> select(id, mna)
desc_rndm  <- read_tsv(".../Smo_random_SD.tab")     |> select(id, mna)
descriptors <- bind_rows(desc_cntrl, desc_slctd, desc_rndm)

### Process
## Transform descriptors to binary matrix
desc_bin <- descriptors |> mutate(mna = str_trim(mna), exist = 1) |>
						   separate_longer_delim(mna, delim = "\r\n") |>
						   pivot_wider(names_from = "mna", values_from = "exist") |>
						   mutate_all(~replace_na(., 0)) |>
						   arrange(id) |>
						   column_to_rownames('id') |>
						   mutate_all(~as.integer(.))
names_desc <- paste0("mna__", 1:ncol(desc_bin))
names(desc_bin) <- names_desc
desc_bin <- desc_bin |> rownames_to_column(var = "compound_id")
## Select 15 descriptors allowing to distinguish between the given cmpnds well
n_descriptors <- 15
mna_select <- MRMR( desc_bin |> select(starts_with("mna")), desc_bin |> pull(compound_id) |> as.factor(), n_descriptors)
mna_select_names  <- names(mna_select$selection)
desc_selected <- desc_bin |> select(compound_id, any_of(mna_select_names)) |>
							 unite("mna_id", starts_with("mna__"), sep = "-", remove = FALSE)
## Join with the main data
data <- data_main |> inner_join(desc_selected)
## Check info loss
data_desc <- data |> select(mna_id, compound_id, category, starts_with("mna__")) |>
						group_by(mna_id) |>
						mutate(compound_id = str_flatten(compound_id, collapse = "_"),
								category = str_flatten(category |> sort() |> unique(), collapse = "_")) |>
						slice_head(n = 1) |>
						ungroup()
## Generate all the points from feature space without the already obtained
featureSpace <- thetaComb(theta = c(0,1), nfact = n_descriptors, intercept = FALSE) |> as_tibble(.name_repair = "universal")
names(featureSpace) <- mna_select_names
featureSpace_sample <- featureSpace |> unite("mna_id", starts_with("mna__"), sep = "-", remove = FALSE) |>
					  					mutate(compound_id = row_number()) |>
					  					rowwise() |>
					  					mutate(compound_id = str_c(c("sample", compound_id), collapse = "_")) |>
					  					ungroup() |>
					  					mutate(category = "non_existent") |>
					  					select(compound_id, mna_id, category, starts_with("mna__")) |>
					  					anti_join(data_desc, by = "mna_id") |>
					  					slice_sample(prop = .1)

### UMAP the data
set.seed(78)
data_2umap <- bind_rows(data_desc, featureSpace_sample) |> mutate_at(vars(starts_with("mna__")), as.integer)
data_umap  <- umap(data_2umap |> select(starts_with("mna__")) |> as.data.frame(), preserve.seed = TRUE, metric = "manhattan", n_neighbors = n_descriptors)
umap_coordinates <- data_umap$layout |> bind_cols(data_2umap |>
											select(compound_id), data_2umap |>
											select(mna_id), data_2umap |>
											select(category)) |>
											rename(umap_1 = `...1`, umap_2 = `...2`) |>
											mutate(order = case_when(
													category == "non_existent" ~ as.integer(0),
													category == "random" ~ as.integer(1),
													category == "random_selected" ~ as.integer(2),
													category == "selected" ~ as.integer(3),
													category == "control" ~ as.integer(4)
												)
											)

### Prepare and plot the data
data_2vis <- umap_coordinates |> select(mna_id, compound_id, category, umap_1, umap_2, order) |>
									mutate(category = fct_reorder(category, order)) |>
									mutate(name = NA) |>
									mutate(name = if_else(compound_id == "CHEMBL2043437", "GLASDEGIB", name)) |>
									mutate(name = if_else(compound_id == "CHEMBL2105737", "SONIDEGIB", name)) |>
									mutate(name = if_else(compound_id == "CHEMBL473417", "VISMODEGIB", name))
umap_plot  <- ggplot(data_2vis |> filter(category != "non_existent"), aes(x = umap_1, y = umap_2, label = name)) +
															 geom_point(aes(fill = category, alpha = category), shape=21, size=2.5) +
															 scale_fill_manual(values = c(.4, .6, .8, 1)) +
															 scale_fill_manual(values = c("#8a6e59", "#fec369", "#cb487a", "#7bc576")) +
															 geom_text(vjust = 2, hjust = .6, size = 2.5) +
															 theme_classic() +
															 theme(legend.position = "right") +
															 coord_fixed()
umap_plot

### Export the results
write_tsv(data_2vis, ".../16-02-26__selected_cmpnds_umap.tab")
ggsave(".../16-02-26__extended_umap.png",
			plot = umap_plot, width = 7, height = 5, units = "in", dpi = 300, bg = "white")