library(tidyverse)

# Input
phpt  <- read_tsv("C:/Users/XPS/Documents/IBMC/2025/VVP/cervicalCancer_SHH_MMP/fromColleagues/data/results_ppv.tab") |>
					mutate(category = if_else(category == "selected", "selected_phpt", category))
p2024 <- read_file("C:/Users/XPS/Documents/IBMC/2025/VVP/cervicalCancer_SHH_MMP/fromColleagues/data/CHR_SMO_2540__DAF.SDF") |>
					str_replace_all("\r", "\n") |>
					str_replace_all("\n\n", "\n") |>
					str_replace_all("\n", "\r\n") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$") |>
					rowwise() |>
					mutate(compound_id = str_match(value, ">  <IDNUMBER>.*\r\n(.*)\r\n")[2] |> str_trim(),
							category = "selected_pass",
							collection = str_match(value, ">  <Collection>.*\r\n(.*)\r\n")[2] |> str_trim(),
							pass_2024_raw = str_match(value, ">  <PASS_ACTIVITY_SPECTRUM>.*\r\n(.*)\r\n")[2] |> str_trim(),
							availability = str_match(value, ">  <available>.*\r\n(.*)\r\n")[2] |> str_trim() |> as.numeric(),
							purity = str_match(value, ">  <Purity>.*\r\n(.*)\r\n")[2] |> str_trim(),
							logP = str_match(value, ">  <logP>.*\r\n(.*)\r\n")[2] |> str_trim() |> as.numeric(),
							logD = str_match(value, ">  <logD>.*\r\n(.*)\r\n")[2] |> str_trim() |> as.numeric(),
							logSw = str_match(value, ">  <logSw>.*\r\n(.*)\r\n")[2] |> str_trim() |> as.numeric(),
							pka_ma = str_match(value, ">  <pKa_ma>.*\r\n(.*)\r\n")[2] |> str_trim() |> as.numeric(),
							color = str_match(value, ">  <Color>.*\r\n(.*)\r\n")[2] |> str_trim(),
							stereo = str_match(value, ">  <STEREO>.*\r\n(.*)\r\n")[2] |> str_trim(),
							state = str_match(value, ">  <State>.*\r\n(.*)\r\n")[2] |> str_trim(),
							cas = str_match(value, ">  <CAS>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
					ungroup() |>
					filter(!is.na(compound_id)) |>
					mutate(pass_2024 = as.numeric(str_sub(pass_2024_raw, 1, 5)) - as.numeric(str_sub(pass_2024_raw, 8, 12))) |>
					select(compound_id, category, collection, pass_2024, availability, purity, logP, logD, logSw, pka_ma, color, stereo, state, cas)

result <- p2024 |> full_join(phpt, keep = FALSE) |>
						mutate(external_id = if_else(is.na(external_id), str_c(row_number(), "pass2024", sep = "_"), external_id)) |>
						arrange(desc(propensity), desc(main_d), desc(main_i), desc(local_i), desc(pass_2024)) |>
						select(compound_id, category, collection, pass_2024, local_i, main_i, main_d, propensity, availability, purity, logP, logD, logSw, pka_ma, color, stereo, state, cas, chembl_similars, n_similars, external_id)

# Output
write_tsv(result, "C:/Users/XPS/Documents/IBMC/2025/VVP/cervicalCancer_SHH_MMP/fromColleagues/data/results_daf_ppv__01-04-26.tab")