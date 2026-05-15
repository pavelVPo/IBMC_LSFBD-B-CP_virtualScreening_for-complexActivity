library(tidyverse)
library(umap)
library(praznik)
library(mirt)

## Input
# Selected
selected_ids <- read_tsv("C:/.../data/compounds_ids.tab")
sdf_mna_1 <- read_file("C:/.../data/all_CR_raw_SD.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, regex(">  <MNA_DESCRIPTORS>[^\r\n]*\r\n(.*)\r\n\r\n", multiline = TRUE, dotall = TRUE))[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna_2 <- read_file("C:/.../data/all_CR_raw_SD_I.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, regex(">  <MNA_DESCRIPTORS>[^\r\n]*\r\n(.*)\r\n\r\n", multiline = TRUE, dotall = TRUE))[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna_3 <- read_file("C:/.../data/all_CR_raw_SD_II.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, regex(">  <MNA_DESCRIPTORS>[^\r\n]*\r\n(.*)\r\n\r\n", multiline = TRUE, dotall = TRUE))[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna_4 <- read_file("C:/.../data/all_CR_raw_SD_III.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, regex(">  <MNA_DESCRIPTORS>[^\r\n]*\r\n(.*)\r\n\r\n", multiline = TRUE, dotall = TRUE))[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna      <- bind_rows(sdf_mna_1, sdf_mna_2, sdf_mna_3, sdf_mna_4)
selected_mna <- inner_join(selected_ids, sdf_mna) |> mutate(category = "selected") |> select(id, molfile, category, mna)
# Random compounds from CR
set.seed(98)
random_mna <- sdf_mna |> slice_sample(n = 200) |> mutate(category = "random") |> select(id, molfile, category, mna)
# 3 RM_D from DAF
rmd_mna <- read_file("C:/.../data/3_RM_SD.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, "> <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, regex(">  <MNA_DESCRIPTORS>[^\r\n]*\r\n(.*)\r\n\r\n", multiline = TRUE, dotall = TRUE))[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct() |>
						 mutate(category = "rmd") |>
						 select(id, molfile, category, mna)
# 5_Smo_RA_Pareto_SD from DAF
pareto_mna <- read_file("C:/.../data/5_Smo_RA_Pareto_SD.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <MOLID>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, regex(">  <MNA_DESCRIPTORS>[^\r\n]*\r\n(.*)\r\n\r\n", multiline = TRUE, dotall = TRUE))[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct() |>
						 mutate(category = "pareto") |>
						 select(id, molfile, category, mna)

## Process
# gather
mnas_mol <- bind_rows(selected_mna, random_mna, rmd_mna, pareto_mna) |>
					group_by(mna) |>
					mutate(id = str_c(id, collapse = " abc "), 
						   category = str_c(category, collapse = " abc "),
						   mna_str = mna) |>
					slice_head(n = 1) |>
					ungroup()
mnas <- mnas_mol |> select(id, category, mna, mna_str) |>
					separate_longer_delim(mna, delim="\r\n") |>
					pivot_wider(names_from=mna, values_from=mna)
sdf <-  mnas_mol |> mutate(id_field = "\r\n>  <compound_id>\r\n",
						  category_field = "\r\n\r\n>  <category>\r\n",
						  end_rec = "\r\n\r\n$$$$") |>
					select(molfile, id_field, id, category_field, category, end_rec) |>
					unite("record", molfile:end_rec, sep = "")
write_lines(str_c("", sdf[[1]] |> str_replace("\r\n", "")), "C:/.../data/results_for_extUMAP.SDF", sep = "\r\n")

# names
names_mna <- paste0("mna_", -3:ncol(mnas))
names_mna[1] <- "id"
names_mna[2] <- "category"
names_mna[3] <- "mnastr"
names(mnas) <- names_mna
# Convert to int
mnas <- mnas |> mutate_at(vars(starts_with("mna_")), ~if_else(!is.na(.), 1, 0))
# Prepare for the extended UMAP
n_descriptors <- 11
mna_select <- MRMR( mnas |> select(starts_with("mna_")), mnas |> pull(id) |> as.factor(), n_descriptors)
mna_select_names  <- names(mna_select$selection)
mna_selection <- mnas |> select(id, category, mnastr, any_of(mna_select_names)) |>
						 rowwise() |>
						 mutate(id_mna = str_c(c_across(starts_with("mna_")), collapse = "-")) |>
						 ungroup()
mna_c_names <- names(mna_selection)[5:length(names(mna_selection))-1]
# Prepare for the visualization
mna_2vis <- mna_selection |> group_by(id_mna) |>
							  mutate(category = unique(category) |> sort() |> str_c(collapse = " & "),
							  		 id = unique(id) |> sort() |> str_c(collapse = " & "),
							  		 mnastr = unique(mnastr) |> sort() |> str_c(collapse = " & ")) |>
							  slice_head(n = 1) |>
							  ungroup()
# Check
mna_2vis_check <- mna_2vis |> group_by(category) |> summarise(n = n()) |> arrange(n)
mna_2vis_check
# Extend the feature space coverage
# feature space
featureSpace <- thetaComb(theta = c(0,1), nfact = n_descriptors, intercept = FALSE) |> as_tibble(.name_repair = "universal") |> unite(id_mna, starts_with("..."), sep = "-") |>
					mutate(id = row_number()) |>
					rowwise() |>
					mutate(id = str_c(c("sample", id), collapse = "_"),
							mnastr = NA) |>
					ungroup() |>
					mutate(category = as.factor("yet_unknown")) |>
					select(id, category, mnastr, id_mna)
featureSpace_sample <- featureSpace |> anti_join(mna_2vis, by = "id_mna") |> slice_sample(prop = .2) |>
							separate_wider_delim(id_mna, delim="-", names = mna_c_names, cols_remove = FALSE) |>
							select(id, category, mnastr, starts_with("mna_"), id_mna) |>
							mutate(across(starts_with("mna_"), ~ as.integer(.x)))
# Gather
data <- bind_rows(mna_2vis, featureSpace_sample) |> arrange(category)

# Extended UMAP
set.seed(88)
data_umap <- umap( data |> select(starts_with("mna_")) |> as.data.frame(),
					preserve.seed = TRUE,
					metric = "manhattan",
					n_neighbors = n_descriptors - 1 )
umap_coordinates <- data_umap$layout |> bind_cols(data |> select(category),
												  data |> select(id)) |>
										rename(UMAP_1 = `...1`, UMAP_2 = `...2`)
# Plot the results
umap_plot  <- ggplot(umap_coordinates |> filter(category != "yet_unknown"), aes(x = UMAP_1, y = UMAP_2)) +
															 geom_point(aes(fill = category, size = category), color = "black", shape=21) +
															 scale_fill_manual(values = c("#568b87", "#b5d1ae", "grey80",
															 							  "#90caf9", "#FAC881FF", "#42a5f5", "#973B2B")) +
															 scale_size_manual(values = c(3, 3, 1, 3, 3, 3, 3)) +
															 theme_classic() +
															 theme(legend.position = "right") +
															 coord_fixed()
umap_plot

# Export plot
ggsave("C:/.../pics/15-05-26_SMO-search_extUMAP.png", width = 6, height = 5, units = "in", dpi = 300)
# Export coordinates
write_tsv(umap_coordinates, "C:/.../data/15-05-26_SMO-search_extUMAP.tsv")