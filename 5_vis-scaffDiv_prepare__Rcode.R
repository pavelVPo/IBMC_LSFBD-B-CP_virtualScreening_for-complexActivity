library(tidyverse)

## Input
# records on testing
smo_cmpnd_tested <- read_tsv(".../data/studied_all__list.tsv") |>
						filter(target_id == 101400) |>
						select(compound_id) |>
						distinct()
smo_cmpnd_predicted <- read_tsv(".../data/2F-CV_activities.tsv") |>
									filter(smo > 0) |>
									select(compound_id) |>
									distinct()
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
							 select(compound_id, n_descriptors, new_descriptors, smo) |>
							 collect() |>
							 #filter(smo > 0) |>
							 select(compound_id, smo) |>
							 distinct()
# SDFs
main_sdf_1 <- read_file(".../data/first_set.SDF") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$")
main_sdf_2 <- read_file(".../data/second_set.SDF") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$")
main_sdf <- bind_rows(main_sdf_1, main_sdf_2)
main_sdf <- main_sdf |> rowwise() |>
					mutate(compound_id = str_match(value, "END\r\n>  <id>\r\n(.*)\r\n")[2]) |>
					ungroup()
cr_sdf <- read_file(".../data/BMS_300K_base.sdf") |>
					 as_tibble() |>
					 separate_longer_delim(value, delim = "$$$$") |>
					 rowwise() |>
					 mutate(compound_id = str_match(value, "END\n>  <idnumber>.*\n(.*)\n")[2]) |>
					 ungroup()

# Select the structures
smo_cmpnd_tested_sdf <- smo_cmpnd_tested |>
						mutate(compound_id = as.character(compound_id)) |>
						inner_join(main_sdf) |>
						group_by(compound_id) |>
						slice_head(n = 1) |>
						ungroup() |>
						select(value)
smo_cmpnd_predicted_sdf <- smo_cmpnd_predicted |>
						mutate(compound_id = as.character(compound_id)) |>
						inner_join(main_sdf) |>
						group_by(compound_id) |>
						slice_head(n = 1) |>
						ungroup() |>
						select(value)
smo_cr_sdf <- smo_cmpnd_predicted_cr |> 
						inner_join(cr_sdf) |>
						filter(smo > 0) |>
						mutate(value = str_replace_all(value, "\n", "\r\n")) |>
						select(value)

## Output
set.seed(548)
smo_cmpnd_tested_sdf_2w <- smo_cmpnd_tested_sdf |>
						   mutate(end_rec = "$$$$",
						   		  value = str_trim(value, side = "left")) |>
						   unite("record", value:end_rec, sep = "")
smo_cmpnd_predicted_sdf_2w <- smo_cmpnd_predicted_sdf |>
						   		slice_sample(n = 9999) |>
						   		mutate(end_rec = "$$$$",
						   		  value = str_trim(value, side = "left")) |>
						   		unite("record", value:end_rec, sep = "")
smo_cr_sdf_2w <- smo_cr_sdf |> slice_sample(n = 9999) |>
						   mutate(end_rec = "$$$$",
						    	value = str_trim(value, side = "left")) |>
						   unite("record", value:end_rec, sep = "")

write_lines(str_c("\r\n", str_c(smo_cmpnd_tested_sdf_2w[[1]], collapse="\r\n\r\n")), ".../data/sdf_2_vdiv/smo_studied_cmpnds.SDF")
write_lines(str_c("\r\n", str_c(smo_cmpnd_predicted_sdf_2w[[1]], collapse="\r\n\r\n")), ".../data/sdf_2_vdiv/smo_predicted_cmpnds__sample.SDF")
write_lines(str_c("\r\n", str_c(smo_cr_sdf_2w[[1]], collapse="\r\n\r\n")), ".../data/sdf_2_vdiv/smo_cr_cmpnds__sample.SDF")

## Changes of diversity on the level of scaffolds are due to the at least 2 factors:
# In silico selection of compounds
# Increase of the set's size