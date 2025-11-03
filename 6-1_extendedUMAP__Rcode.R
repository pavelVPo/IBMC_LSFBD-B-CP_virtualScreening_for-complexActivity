library(tidyverse)
library(umap)
library(praznik)
library(mirt)

## Input
studied_raw <- read_file(".../data/sdf_2_umap/smo_studied_cmpnds__2mna_SD.SDF") |> as_tibble() |>
							separate_longer_delim(value, delim = "$$$$") |>
							rowwise() |>
							mutate(compound_id = str_match(value, "END\r\n>  <id>\r\n(.*)\r\n")[2] |> str_trim(),
									probable = "totally",
									probable_lab = "yes",
									source = "chembl" ) |>
							ungroup() |>
					 		separate_wider_delim(value, delim = "\r\n>  <MNA_DESCRIPTORS>\r\n", names = c("value", "mna"), too_many = "merge", too_few = "align_start") |>
					 		mutate(mna = str_trim(mna)) |>
					 		select(compound_id, probable, probable_lab, source, mna) |>
					 		separate_longer_delim(mna, delim="\r\n") |>
					 		pivot_wider(names_from=mna, values_from=mna)
predicted_raw <- read_file(".../data/sdf_2_umap/smo_predicted_cmpnds__2mna_SD.SDF") |> as_tibble() |>
							separate_longer_delim(value, delim = "$$$$") |>
							rowwise() |>
							mutate(compound_id = str_match(value, "END\r\n>  <id>\r\n(.*)\r\n")[2] |> str_trim(),
									probable = str_match(value, "\r\n>  <probable>\r\n(.*)\r\n")[2] |> str_trim(),
									probable_lab = str_match(value, "\r\n>  <probable_lab>\r\n(.*)\r\n")[2] |> str_trim(),
									source = "chembl" ) |>
							ungroup() |>
					 		separate_wider_delim(value, delim = "\r\n>  <MNA_DESCRIPTORS>\r\n", names = c("value", "mna"), too_many = "merge", too_few = "align_start") |>
					 		mutate(mna = str_trim(mna)) |>
					 		select(compound_id, probable, probable_lab, source, mna) |>
					 		separate_longer_delim(mna, delim="\r\n") |>
					 		pivot_wider(names_from=mna, values_from=mna)
cr_raw <- read_file(".../data/sdf_2_umap/smo_cr_cmpnds__2mna_SD.SDF") |> as_tibble() |>
							separate_longer_delim(value, delim = "$$$$") |>
							rowwise() |>
							mutate(compound_id = str_match(value, "END\r\n>  <id>\r\n(.*)\r\n")[2] |> str_trim(),
									probable = str_match(value, "\r\n>  <probable>\r\n(.*)\r\n")[2] |> str_trim(),
									probable_lab = str_match(value, "\r\n>  <probable_lab>\r\n(.*)\r\n")[2] |> str_trim(),
									source = "chemrar" ) |>
							ungroup() |>
					 		separate_wider_delim(value, delim = "\r\n>  <MNA_DESCRIPTORS>\r\n", names = c("value", "mna"), too_many = "merge", too_few = "align_start") |>
					 		mutate(mna = str_trim(mna)) |>
					 		select(compound_id, probable, probable_lab, source, mna) |>
					 		separate_longer_delim(mna, delim="\r\n") |>
					 		pivot_wider(names_from=mna, values_from=mna)
mna_raw <- bind_rows(predicted_raw, studied_raw, cr_raw) |> mutate(probable_lab = fct(probable_lab)) |> filter(!is.na(compound_id))

## Process
names_mna <- paste0("mna_", -3:ncol(mna_raw))
names_mna[1] <- "compound_id"
names_mna[2] <- "probable"
names_mna[3] <- "probable_lab"
names_mna[4] <- "source"
names(mna_raw) <- names_mna
# Convert to int
mna_raw <- mna_raw |> mutate_at(vars(starts_with("mna_")), ~if_else(!is.na(.), 1, 0))
# Prepare for the extended UMAP
n_descriptors <- 15
mna_select <- MRMR( mna_raw |> select(starts_with("mna")), mna_raw |> pull(compound_id) |> as.factor(), n_descriptors)
mna_select_names  <- names(mna_select$selection)
mna_selection <- mna_raw |> select(compound_id, probable, probable_lab, source, any_of(mna_select_names))
# Prepare further
data_2vis_raw <- mna_selection |> unite(mna_id, starts_with("mna"), sep = "-") |>
								group_by(mna_id) |>
								mutate(category = unique(probable) |> sort() |> str_c(collapse = "|")) |>
								mutate(source = unique(source) |> sort() |> str_c(collapse = "|")) |>
								mutate(compound_id = unique(compound_id) |> sort() |> str_c(collapse = "|")) |>
								ungroup() |>
								select(compound_id, category, source, mna_id) |>
								distinct()
# Check
data_2vis_check <- data_2vis_raw |> group_by(category) |> summarise(n = n()) |> arrange(n)
data_2vis_check
# Prepare even further
data_2vis <- data_2vis_raw |> mutate(category = as.factor(category)) |>
							  mutate(mna_vec = mna_id) |>
										separate_wider_delim(mna_vec, delim="-", names=c("mna__1", "mna__2", "mna__3", "mna__4", "mna__5", "mna__6",
										"mna__7","mna__8","mna__9","mna__10","mna__11","mna__12","mna__13","mna__14","mna__15"))
# feature space
featureSpace <- thetaComb(theta = c(0,1), nfact = n_descriptors, intercept = FALSE) |> as_tibble(.name_repair = "universal") |> unite(mna_id, starts_with("..."), sep = "-") |>
					mutate(compound_id = row_number()) |>
					rowwise() |>
					mutate(compound_id = str_c(c("sample", compound_id), collapse = "_")) |>
					ungroup() |>
					mutate(category = as.factor("yet_unknown"), source = "yet unknown") |>
					select(compound_id, category, source, mna_id)
featureSpace_sample <- featureSpace |> anti_join(data_2vis, by = "mna_id") |> slice_sample(prop = .05) |>
										mutate(mna_vec = mna_id, order = 7) |>
										separate_wider_delim(mna_vec, delim="-", names=c("mna__1", "mna__2", "mna__3", "mna__4", "mna__5", "mna__6",
										"mna__7","mna__8","mna__9","mna__10","mna__11","mna__12","mna__13","mna__14","mna__15"))
data_2vis_ready <- bind_rows(data_2vis, featureSpace_sample) |> mutate_at(vars(starts_with("mna__")), as.integer) |> select(-order)
# UMAP the data and plot the results
set.seed(88)
data_umap <- umap( data_2vis_ready |> select(starts_with("mna__")) |> as.data.frame(), preserve.seed = TRUE, metric = "manhattan", n_neighbors = n_descriptors  )
umap_coordinates <- data_umap$layout |> bind_cols(data_2vis_ready |>
										select(category), data_2vis_ready |>
										select(compound_id), data_2vis_ready |>
										select(source)) |>
										rename(UMAP_1 = `...1`, UMAP_2 = `...2`)
umap_plot  <- ggplot(umap_coordinates |> filter(category != "yet_unknown"), aes(x = UMAP_1, y = UMAP_2)) +
															 geom_point(aes(fill = source, alpha = category), shape=21) +
															 scale_fill_manual(values = c("#8FB9A8", "#FAC881FF", "#383788FF")) +
															 scale_alpha_manual(values = c(.5, .75, .75, 1, 1, 1)) +
															 theme_classic() +
															 theme(legend.position = "right") +
															 coord_fixed()
umap_plot

# Export plot
ggsave(".../pics/extUMAP_propensity.png", plot = umap_plot, width = 7, height = 5, units = "in", dpi = 300)