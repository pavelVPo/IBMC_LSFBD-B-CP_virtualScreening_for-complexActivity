library(tidyverse)

## Input
current  <- read_tsv("C:/.../current_results.tab")
cs_s <- read_file(".../207selected_CR_09-02-26_SD.SDF") |>
					str_replace_all("\r", "\n") |>
					str_replace_all("\n\n", "\n") |>
					str_replace_all("\n", "\r\n") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$") |>
					rowwise() |>
					mutate(mol = str_extract(value, regex(".*M  END\r\n", multiline = TRUE, dotall = TRUE)),
						   compound_id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
					select(mol, compound_id) |>
					filter(!is.na(compound_id))
scecrbf_select  <- read_file("C:/.../SC_predictions_results_selected_SLA__17-03-2026.sdf") |>
					str_replace_all("\r", "\n") |>
					str_replace_all("\n\n", "\n") |>
					str_replace_all("\n", "\r\n") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$") |>
					rowwise() |>
					mutate(mol = str_extract(value, regex(".*M  END\r\n", multiline = TRUE, dotall = TRUE)),
							scr = str_match(value, ">  <scr>.*\r\n(.*)\r\n")[2] |> str_trim(),
							scrrbf = str_match(value, ">  <scrrbf>.*\r\n(.*)\r\n")[2] |> str_trim(),
							scec1 = str_match(value, ">  <scec1>.*\r\n(.*)\r\n")[2] |> str_trim(),
							scec2 = str_match(value, ">  <scec2>.*\r\n(.*)\r\n")[2] |> str_trim(),
							scec3 = str_match(value, ">  <scec3>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
					ungroup() |>
					mutate(tempid = row_number(), category = "selected") |>
					filter(!is.na(mol)) |>
					select(mol, tempid, scr, scrrbf, scec1, scec2, scec3)
scecrbf <- scecrbf_select

## Get the IDs for the predicted compounds
scecrbf_ids <- scecrbf |> inner_join(cs_s, by = "mol") |> select(-mol, -tempid)

## Join the results
result <- current |> left_join(scecrbf_ids, by = "compound_id") |>
						arrange(desc(scr), desc(scrrbf), desc(scec1), desc(scec2), desc(scec3), desc(propensity), desc(main_d), desc(main_i), desc(local_i), desc(pass_2024)) |>
						select(compound_id, category, collection, pass_2024, local_i, main_i, main_d, scr, scrrbf, scec1, scec2, scec3, propensity, availability, purity, logP, logD, logSw, pka_ma, color, stereo, state, cas, chembl_similars, n_similars, external_id)

# Output
write_tsv(result, "C:/.../results__03-04-26.tab")