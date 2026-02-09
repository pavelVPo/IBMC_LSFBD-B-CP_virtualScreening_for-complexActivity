library(tidyverse)


### Import structures
## Import the original data after MNA
sdf_mna_1 <- read_file("C:/.../all_CR_raw_SD.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, ">  <MNA_DESCRIPTORS>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna_2 <- read_file("C:/.../all_CR_raw_SD_I.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, ">  <MNA_DESCRIPTORS>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna_3 <- read_file("C:/.../all_CR_raw_SD_II.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, ">  <MNA_DESCRIPTORS>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna_4 <- read_file("C:/.../all_CR_raw_SD_III.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, ">  <MNA_DESCRIPTORS>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna <- bind_rows(sdf_mna_1, sdf_mna_2, sdf_mna_3, sdf_mna_4)

### Import selected IDs
cs_selected <- read_tsv("C:/.../smo_predictions/smo_data_selected.tsv")

### Get structures
cs_random <- sdf_mna |> slice_sample(n = 207) |> select(molfile, id_fld, id, end_rec)
cs_predicted  <- cs_selected |> inner_join(sdf_mna, by = "id") |>
						mutate(local_fld = "\r\n>  <activity_local>\r\n", i_fld = "\r\n>  <activity_i>\r\n",
								d_fld = "\r\n>  <activity_d>\r\n", propensity_fld = "\r\n>  <propensity>\r\n") |>
						select(molfile, id_fld, id, local_fld, activity_local_i, i_fld, activity_i, d_fld, activity_d, propensity_fld, propensity, end_rec)

### Export the sets
# Write Predicted structures (207)
sdf_predicted <- cs_predicted |> unite("record", molfile:end_rec, sep = "")
write_lines(sdf_predicted[[1]], "C:/.../207selected_CR_09-02-26.SDF", sep = "")
# Write Random structures (207)
sdf_random <- cs_random |> unite("record", molfile:end_rec, sep = "")
write_lines(sdf_random[[1]], "C:/.../207random_CR_09-02-26.SDF", sep = "")