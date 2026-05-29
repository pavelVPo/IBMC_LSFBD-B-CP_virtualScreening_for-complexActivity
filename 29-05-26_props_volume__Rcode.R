library(tidyverse)
library(rcdk)

## Input
# Selected
selected <- read_file("C:/.../15_SMO_SELECTED (PASS2024CellLines).SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(compound_id = str_match(value, ">  <IDNUMBER>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mol = str_match(value, regex(".*M  END\r\n", dotall = TRUE))[1] |> unlist()) |>
						 ungroup() |>
						 select(compound_id, mol) |>
						 filter(!is.na(mol))
# Find them in the properties file and calculate the possible molar concentration
all_props <- read_tsv("C:/.../tabular_props.tab")

## Process the results
data <- all_props |> inner_join(selected, by = c("id" = "compound_id")) |>
				rowwise() |>
				mutate(rcdk_exactmass_da = get.exact.mass(parse.smiles(smiles)[[1]]),
						weight_g = .005) |>
				ungroup() |>
				select(id, name, smiles, link, p_sp3, stereo, state, purity, salt,
						H_acceptor, H_donor, B_rotN, N_O, logP, logD,
						logSw, pka_ma, pka_mb, n_chirals, rcdk_exactmass_da, weight_g)
data_concentration <- data |> mutate(conc_M = weight_g / rcdk_exactmass_da) |>
							  mutate(VL_at_100_uM = weight_g / (.0001*rcdk_exactmass_da),
							  		 VL_at_500_uM = weight_g / (.0005*rcdk_exactmass_da))

## Export the results with concentration
write_tsv(data_concentration, "C:/.../SMO_15_selected_VL.tsv")