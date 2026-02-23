library(tidyverse)


### Import
# Selected cmpnds
data        <- read_tsv("C:/.../13-02-26__selected_cmpnds_described.tab")
data_vis    <- read_tsv("C:/.../16-02-26__selected_cmpnds_umap.tab", guess_max = 5000) |>
					mutate(compound_id_ = compound_id) |>
					separate_longer_delim(compound_id, delim = "_")

### Select
# Main_d
data_maind <- data |> filter(category == "selected" &
							 0 < logP & logP < 5 &
							 0 < logP & logP < 5 &
							 -5.7 < logSw & logSw < -1) |>
					  arrange(desc(propensity)) |>
					  arrange(desc(main_d))
# Main_i
data_maini <- data |> filter(category == "selected" &
							 0 < logP & logP < 5 &
							 0 < logP & logP < 5 &
							 -5.7 < logSw & logSw < -1) |>
					  arrange(desc(propensity)) |>
					  arrange(desc(main_i))
# local_i
data_locali <- data |> filter(category == "selected" &
							 0 < logP & logP < 5 &
							 0 < logP & logP < 5 &
							 -5.7 < logSw & logSw < -1) |>
					  arrange(desc(propensity)) |>
					  arrange(desc(local_i))

### Prepare to plot
slct_maind <- data_maind |> filter(main_d > .7) |> select(compound_id) |> mutate(type = "main-D")
slct_maini <- data_maind |> filter(main_i > .7) |> select(compound_id) |> mutate(type = "main-I")
slct_locali <- data_maind |> filter(local_i > .7) |> select(compound_id) |> mutate(type = "local-I")
data_2vis <- data_vis |> left_join(slct_maind, by = "compound_id") |> left_join(slct_maini, by = "compound_id") |>
						 left_join(slct_locali, by = "compound_id") |>
						 rowwise() |>
						 mutate(type = str_c(c(type, type.x, type.y) |> sort() |> unique(), collapse = "_")) |>
						 ungroup() |>
						 mutate(type = if_else(type == "", "not_selected", type)) |>
						 select(-compound_id) |>
						 distinct() |>
						 mutate(compound_id = compound_id_) |>
						 select(-compound_id_) |>
						 group_by(compound_id) |>
						 mutate(type = str_c(type |> sort() |> unique(), collapse = "_")) |>
						 slice_head(n = 1) |>
						 ungroup() |>
						 mutate(order_type = case_when(
													type == "not_selected" ~ as.integer(0),
													type == "local-I_not_selected" ~ as.integer(1),
													type == "main-I_not_selected" ~ as.integer(2),
													type == "main-D" ~ as.integer(3),
													type == "main-I" ~ as.integer(4),
													type == "local-I" ~ as.integer(5)
												)
											) |>
						 mutate(type = fct_reorder(type, order_type))

### Plot
umap_plot  <- ggplot(data_2vis |> filter(category != "non_existent"), aes(x = umap_1, y = umap_2, label = name)) +
															 geom_point(aes(alpha = type), fill = "#685D79", shape=21, size=2.5) +
															 scale_alpha_manual(values = c(.2, .5, .5, 1, 1, 1)) +
															 geom_text(vjust = 2, hjust = .6, size = 2.5) +
															 theme_classic() +
															 theme(legend.position = "right") +
															 coord_fixed() +
															 guides(fill=FALSE)
umap_plot
umap_plot_dop  <- ggplot(data_2vis |> filter(category != "non_existent"), aes(x = umap_1, y = umap_2, label = name)) +
															 geom_point(aes(fill = type, alpha = type), shape=21, size=2.5) +
															 scale_alpha_manual(values = c(.2, .5, .5, 1, 1, 1)) +
															 scale_fill_manual(values = c("grey80", "mediumblue", "coral", "mediumorchid", "coral", "mediumblue")) +
															 geom_text(vjust = 2, hjust = .6, size = 2.5) +
															 theme_classic() +
															 theme(legend.position = "right") +
															 coord_fixed()
umap_plot_dop


### Export the results
ggsave("C:/.../16-02-26__extended_umap_dop.png",
			plot = umap_plot_dop, width = 7, height = 5, units = "in", dpi = 300, bg = "white")
write_tsv(data_maind  |> select(compound_id, availability, main_d, propensity, state, availability, color, stereo, n_similars, logP, logD, logSw), "C:/.../17-02-26__ordered_mainD_physchem.tab")
write_tsv(data_maini  |> select(compound_id, availability, main_i, propensity, state, availability, color, stereo, n_similars, logP, logD, logSw), "C:/.../17-02-26__ordered_mainI_physchem.tab")
write_tsv(data_locali |> select(compound_id, availability, local_i, propensity, state, availability, color, stereo, n_similars, logP, logD, logSw), "C:/.../17-02-26__ordered_localI_physchem.tab")