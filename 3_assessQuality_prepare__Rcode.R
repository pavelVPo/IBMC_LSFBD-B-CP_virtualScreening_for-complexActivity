library(tidyverse)
library(arrow)

## Input
# use subset SDFs generated on the previous step, since it allows to check the validity of their structure one more time
main_sdf_1 <- read_file(".../data/first_set.SDF") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$")
main_sdf_2 <- read_file(".../data/second_set.SDF") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$")
main_sdf <- bind_rows(main_sdf_1, main_sdf_2)
#		||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#		|||||	
#		|||||	---------+--------------------+---------------+
#		|||||	|What are the targets of interest:			  |
#		|||||	+--------+--------------------+---------------+
#		|||||	| tid    | pref_name          | chembl_id     |
#		|||||	-----------------------------------------------
#		|||||	| 103218 | Ca-Ski             | CHEMBL1075403 |
#		|||||	-----------------------------------------------
#		|||||	| 106482 | C-33-A             | CHEMBL2366313 |
#		|||||	-----------------------------------------------
#		|||||	| 80472  | SiHa               | CHEMBL612542  |
#		|||||	-----------------------------------------------
#		|||||	| 101400 | Protein smoothened | CHEMBL5971    |
#		|||||	+--------+--------------------+---------------+
#		|||||	
#		||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# bind the results of 2-fold CV, these data are quite big, thus, SEE: https://r4ds.hadley.nz/arrow.html
rslt_1_raw <- open_delim_dataset(
					sources = ".../data/first_set__predicted.csv",
					delim = ";"
				)
rslt_1_schema <- schema(rslt_1_raw)
rslt_1_schema[[1]] <- Field$create("<id>", string())
rslt_1_schema[[2]] <- Field$create("Substructure Descriptors", string())
rslt_1_schema[[3]] <- Field$create("New Descriptors", string())
rslt_1_schema[[4]] <- Field$create("Possible Activities at Pa>Pi", string())
rslt_1_raw <- open_delim_dataset(
					sources = ".../data/first_set__predicted.csv",
					skip = 1,
					delim = ";",
					schema = rslt_1_schema
			)
rslt_2_raw <- open_delim_dataset(
					sources = ".../data/second_set__predicted.csv",
					delim = ";"
				)
rslt_2_schema <- schema(rslt_2_raw)
rslt_2_schema[[1]] <- Field$create("<id>", string())
rslt_2_schema[[2]] <- Field$create("Substructure Descriptors", string())
rslt_2_schema[[3]] <- Field$create("New Descriptors", string())
rslt_2_schema[[4]] <- Field$create("Possible Activities at Pa>Pi", string())
rslt_2_raw <- open_delim_dataset(
					sources = ".../data/second_set__predicted.csv",
					skip = 1,
					delim = ";",
					schema = rslt_2_schema
				)
## Process
# Process the data and get the usable table
data <- main_sdf |> rowwise() |>
					mutate(compound_id = str_match(value, "END\r\n>  <id>\r\n(.*)\r\n")[2],
							  target_id = str_match(value, regex("all_targets>\r\n(.*?)\r\n\r", dotall = TRUE))[2]) |>
					ungroup() |>
					separate_longer_delim(target_id, delim = "\r\n") |>
					mutate(activity = 1) |>
					select(-value)
# Export the testing data from SDFs
write_tsv(data, ".../data/studied_all__list.tsv", quote = "needed")
# Process the results of prediction and convert to the long format for the future assessemnt
# Unfortunately, tidyr is not available to deal with such objects,
# Thus, pivoting should be done after the data's collection or using other approaches (time-consuming)
rslt_1 <- rslt_1_raw |> mutate(compound_id = as.character(`<id>`),
										n_descriptors = `Substructure Descriptors`,
										new_descriptors = `New Descriptors`,
										n_activities = `Possible Activities at Pa>Pi`,
										caski = gsub(",", ".", `103218`) |> as.numeric(),
										c33a = gsub(",", ".", `106482`) |> as.numeric(),
										siha = gsub(",", ".", `80472`) |> as.numeric(),
										smo = gsub(",", ".", `101400`) |> as.numeric()  ) |>
							 select(compound_id, n_descriptors, new_descriptors, n_activities, caski, c33a, siha, smo) |>
							 collect()
rslt_2 <- rslt_2_raw |> mutate(compound_id = as.character(`<id>`),
										n_descriptors = `Substructure Descriptors`,
										new_descriptors = `New Descriptors`,
										n_activities = `Possible Activities at Pa>Pi`,
										caski = gsub(",", ".", `103218`) |> as.numeric(),
										c33a = gsub(",", ".", `106482`) |> as.numeric(),
										siha = gsub(",", ".", `80472`) |> as.numeric(),
										smo = gsub(",", ".", `101400`) |> as.numeric()  ) |>
							 select(compound_id, n_descriptors, new_descriptors, n_activities, caski, c33a, siha, smo) |>
							 collect()
rslt <- bind_rows(rslt_1, rslt_2)
# Get the meta and the actual results to the separate things
rslt_meta <- rslt |> select(compound_id, n_descriptors, new_descriptors, n_activities)
rslt_act <- rslt |> select(compound_id, caski, c33a, siha, smo)
# Export these results
write_tsv(rslt_meta, ".../data/2F-CV_meta.tsv")
write_tsv(rslt_act, ".../data/2F-CV_activities.tsv")