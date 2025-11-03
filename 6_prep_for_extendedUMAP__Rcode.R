library(tidyverse)
library(arrow)

# In this the SDFs should be created to be used to calculate MNA
# Each SDF should include: structure, source, two-level propensity's classification: yes, no; three-level propensity's classification: totally, yes, no;
# For the ChEMBL-predicted and CR-predicted

## Input
# get the list of actually evaluated compounds to filter the predicted compounds
tested_cmpnds <- read_file(".../data/sdf_2_umap/smo_studied_cmpnds__2mna_SD.SDF") |>
							as_tibble() |>
					 		separate_longer_delim(value, delim = "$$$$") |>
					 		rowwise() |>
					 		mutate(compound_id = str_match(value, "END\r\n>  <id>\r\n(.*)\r\n")[2] |> str_trim()) |>
					 		ungroup() |>
					 		select(compound_id) |>
					 		filter(!is.na(compound_id))
# get the list of all predicted compounds and corresponding Smo values
smo_cmpnd_predicted <- read_tsv(".../data/2F-CV_activities.tsv") |>
									group_by(compound_id) |>
									mutate(smo = mean(smo)) |>
									ungroup() |>
									filter(smo > .5 | smo < -.5) |>
									select(compound_id, smo) |>
									distinct() |>
									mutate(compound_id = as.character(compound_id), probable = if_else(smo < 0, "no", "yes"), probable_lab = if_else(smo < 0, "no", "yes"), source = "ChEMBL") |>
									anti_join(tested_cmpnds, by = "compound_id")
predicted_sample_size <- smo_cmpnd_predicted |> filter(probable != "no") |> nrow()
smo_cmpnd_predicted <- smo_cmpnd_predicted |> group_by(probable) |>
									slice_sample(n = predicted_sample_size) |>
									ungroup()
# get the list of all predicted CR compounds and corresponding Smo values
rslt_cr_raw <- open_delim_dataset(
					sources = ".../data/300k__predicted.csv",
					delim = ";"
				)
rslt_cr_schema <- schema(rslt_cr_raw)
rslt_cr_schema[[1]] <- Field$create("<idnumber>", string())
rslt_cr_schema[[2]] <- Field$create("Substructure Descriptors", string())
rslt_cr_schema[[3]] <- Field$create("New Descriptors", string())
rslt_cr_schema[[4]] <- Field$create("Possible Activities at Pa>Pi", string())
rslt_cr_raw <- open_delim_dataset(
					sources = ".../data/300k__predicted.csv",
					skip = 1,
					delim = ";",
					schema = rslt_cr_schema
			)
smo_cmpnd_predicted_cr <- rslt_cr_raw |> mutate(compound_id = as.character(`<idnumber>`),
										n_descriptors = `Substructure Descriptors`,
										new_descriptors = `New Descriptors`,
										n_activities = `Possible Activities at Pa>Pi`,
										smo = gsub(",", ".", `101400`) |> as.numeric()  ) |>
							 select(compound_id, smo) |>
							 collect() |>
							 filter(smo > .5 | smo < -.5) |>
							 select(compound_id, smo) |>
							 mutate(compound_id = as.character(compound_id), probable = if_else(smo < 0, "no", "yes"), probable_lab = if_else(smo < 0, "no", "yes"), source = "CR") |>
							 distinct()
cr_sample_size <- smo_cmpnd_predicted_cr |> filter(probable != "no") |> nrow()
smo_cmpnd_cr <- smo_cmpnd_predicted_cr |> group_by(probable) |>
									slice_sample(n = cr_sample_size) |>
									ungroup() |>
									mutate(source = "cr")
# SDFs
main_sdf_1 <- read_file(".../data/first_set.SDF") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$")
main_sdf_2 <- read_file(".../data/second_set.SDF") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$")
main_sdf <- bind_rows(main_sdf_1, main_sdf_2)
main_sdf <- main_sdf |> rowwise() |>
					mutate(compound_id = str_match(value, "END\r\n>  <id>\r\n(.*)\r\n")[2] |> str_trim()) |>
					ungroup()
cr_sdf <- read_file(".../data/BMS_300K_base.sdf") |>
					 as_tibble() |>
					 separate_longer_delim(value, delim = "$$$$") |>
					 rowwise() |>
					 mutate(compound_id = str_match(value, "END\n>  <idnumber>.*\n(.*)\n")[2] |> str_trim()) |>
					 ungroup()

## Process SDFs
# ChEMBL
smo_cmpnd_predicted_sdf <- smo_cmpnd_predicted |>
						inner_join(main_sdf) |>
						group_by(compound_id) |>
						slice_head(n = 1) |>
						ungroup()
# CR
smo_cmpnd_cr_sdf <- smo_cmpnd_cr |> 
						inner_join(cr_sdf) |>
						mutate(value = str_replace_all(value, "\n", "\r\n"))
# Prepare for export
export_chembl <- smo_cmpnd_predicted_sdf |>
					mutate(id_lab = "\r\n>  <id>\r\n",
							prob_lab = "\r\n\r\n>  <probable>\r\n",
							lab_lab = "\r\n\r\n>  <probable_lab>\r\n",
							source_lab = "\r\n\r\n>  <source>\r\n",
							end_rec = "\r\n\r\n$$$$",
							value = str_replace(value, "\n", "")) |>
					select(value, id_lab, compound_id, prob_lab, probable_lab, lab_lab, probable, source_lab, source, end_rec) |>
					unite("record", value:end_rec, sep = "")
export_cr <- smo_cmpnd_cr_sdf |>
					mutate(id_lab = "\r\n>  <id>\r\n",
							prob_lab = "\r\n\r\n>  <probable>\r\n",
							lab_lab = "\r\n\r\n>  <probable_lab>\r\n",
							source_lab = "\r\n\r\n>  <source>\r\n",
							end_rec = "\r\n\r\n$$$$",
							value = str_replace(value, regex("END.*", multiline = TRUE, dotall = TRUE), "END") |> str_trim()) |>
					select(value, id_lab, compound_id, prob_lab, probable_lab, lab_lab, probable, source_lab, source, end_rec) |>
					unite("record", value:end_rec, sep = "")

## Export
write_lines(str_c("", str_c(export_chembl[[1]], collapse="\r\n")), ".../data/sdf_2_umap/smo_predicted_cmpnds__2mna.SDF")
write_lines(str_c("\r\n", str_c(export_cr[[1]], collapse="\r\n\r\n")), ".../data/sdf_2_umap/smo_cr_cmpnds__2mna.SDF")