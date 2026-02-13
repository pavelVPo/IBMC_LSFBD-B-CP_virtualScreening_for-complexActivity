library(tidyverse)
library(reticulate)

# Get the ChEMBL RDKit-based tool to assess similarity: https://github.com/chembl/FPSim2
use_python("C:/.../python.exe")
fpsim2 <- import("FPSim2")

### Import properties
props <- read_tsv("C:/.../data/tabular_props.tab")

### Import SDFs
cs_selected   <- read_file("C:/.../data/207selected_CR_09-02-26_SD.SDF") |>
						str_trim() |>
						str_replace_all("\r", "\n") |>
						str_replace_all("\n\n", "\n") |>
						str_replace_all("\n", "\r\n") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$") |>
						filter(value != "") |>
						rowwise() |>
						mutate(
							id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
							name = str_match(value, ">  <name>.*\r\n(.*)\r\n")[2] |> str_trim(),
							local_i = str_match(value, ">  <activity_local>.*\r\n(.*)\r\n")[2] |> str_trim(),
							main_i = str_match(value, ">  <activity_i>.*\r\n(.*)\r\n")[2] |> str_trim(),
							main_d = str_match(value, ">  <activity_d>.*\r\n(.*)\r\n")[2] |> str_trim(),
							propensity = str_match(value, ">  <propensity>.*\r\n(.*)\r\n")[2] |> str_trim(),
							mna = str_match(value, regex(">  <MNA_DESCRIPTORS>(.*)$", dotall = TRUE))[2] |> str_trim()
						) |>
						ungroup()
cs_random     <- read_file("C:/.../data/207random_CR_09-02-26_SD.SDF") |>
						str_trim() |>
						str_replace_all("\r", "\n") |>
						str_replace_all("\n\n", "\n") |>
						str_replace_all("\n", "\r\n") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$") |>
						filter(value != "") |>
						rowwise() |>
						mutate(
							id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
							name = str_match(value, ">  <name>.*\r\n(.*)\r\n")[2] |> str_trim(),
							local_i = str_match(value, ">  <activity_local>.*\r\n(.*)\r\n")[2] |> str_trim(),
							main_i = str_match(value, ">  <activity_i>.*\r\n(.*)\r\n")[2] |> str_trim(),
							main_d = str_match(value, ">  <activity_d>.*\r\n(.*)\r\n")[2] |> str_trim(),
							propensity = str_match(value, ">  <propensity>.*\r\n(.*)\r\n")[2] |> str_trim(),
							mna = str_match(value, regex(">  <MNA_DESCRIPTORS>(.*)$", dotall = TRUE))[2] |> str_trim()
						) |>
						ungroup()
cs_controls   <- read_file("C:/.../data/Smo_controls_SD.SDF") |>
						str_trim() |>
						str_replace_all("\r", "\n") |>
						str_replace_all("\n\n", "\n") |>
						str_replace_all("\n", "\r\n") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$") |>
						filter(value != "") |>
						rowwise() |>
						mutate(
							id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
							name = str_match(value, ">  <name>.*\r\n(.*)\r\n")[2] |> str_trim(),
							local_i = str_match(value, ">  <activity_local>.*\r\n(.*)\r\n")[2] |> str_trim(),
							main_i = str_match(value, ">  <activity_i>.*\r\n(.*)\r\n")[2] |> str_trim(),
							main_d = str_match(value, ">  <activity_d>.*\r\n(.*)\r\n")[2] |> str_trim(),
							propensity = str_match(value, ">  <propensity>.*\r\n(.*)\r\n")[2] |> str_trim(),
							mna = str_match(value, regex(">  <MNA_DESCRIPTORS>(.*)$", dotall = TRUE))[2] |> str_trim(),
							smiles = str_match(value, ">  <smiles>.*\r\n(.*)\r\n")[2] |> str_trim()
						) |>
						ungroup()

### Mutate the ChEMBL data to allow their further processing
read_tsv("C:/.../data/all_chembl36.tab") |>
mutate(chembl_id = str_replace(chembl_id, "CHEMBL", "") |> as.integer()) |>
select(canonical_smiles, chembl_id) |>
write_tsv("C:/.../data/all_chembl36.smi", col_names = FALSE)

### Create FP-db using Morgan FPs, r: 3, fpSize: 512
fpsim2$io$create_db_file(
    mols_source="C:/.../data/all_chembl36.smi",
    filename='all_chembl36.h5',
    mol_format=py_none,
    fp_type='Morgan',
    mol_id_prop='chembl_id'
)

### Assess similarity between the datasets using RDkit
fpdb <- fpsim2$FPSim2Engine('all_chembl36.h5')

### Search for the similar cmpnds in ChEMBL
cs_selected_sim <- cs_selected |> inner_join(props, by = "id") |>
								  rowwise() |>
								  mutate(chembl_similars = tryCatch(fpdb$similarity(smiles, threshold=r_to_py(0.7),
								  			metric='tanimoto', n_workers=r_to_py(1))$tolist() |> map(function(x) paste0("CHEMBL", x[1])),
								  		 error = function(e) NA_character_) |> list() ) |>
								  mutate(n_similars = chembl_similars |> length()) |>
								  ungroup() |>
								  mutate(type = "selected")
cs_random_sim <- cs_random |> inner_join(props, by = "id") |>
								  rowwise() |>
								  mutate(chembl_similars = tryCatch(fpdb$similarity(smiles, threshold=r_to_py(0.7),
								  			metric='tanimoto', n_workers=r_to_py(1))$tolist() |> map(function(x) paste0("CHEMBL", x[1])),
								  		 error = function(e) NA_character_) |> list() ) |>
								  mutate(n_similars = chembl_similars |> length()) |>
								  ungroup() |>
								  mutate(type = "random")
cs_controls_sim <- cs_controls |> rowwise() |>
								  mutate(chembl_similars = tryCatch(fpdb$similarity(smiles, threshold=r_to_py(0.7),
								  			metric='tanimoto', n_workers=r_to_py(1))$tolist() |> map(function(x) paste0("CHEMBL", x[1])),
								  		 error = function(e) NA_character_) |> list() ) |>
								  mutate(n_similars = chembl_similars |> length()) |>
								  ungroup() |>
								  mutate(type = "control")

cs_sim <- bind_rows(cs_selected_sim, cs_random_sim, cs_controls_sim) |> mutate(mna = str_replace_all(mna, "\r\n", "_|_"))

### Export the results
cs_export <- cs_sim |> select(id, type, collection, local_i,
								main_i, main_d, propensity,
								availability, purity, logP, logD, logSw, pka_ma, color, stereo, state, cas, chembl_similars, n_similars) |>
					   rename(compound_id = id, category = type) |>
					   mutate(chembl_similars = if_else(n_similars == 0, "0" |> as.list(), chembl_similars)) |>
					   rowwise() |>
					   mutate(chembl_similars = chembl_similars |> unlist() |> str_c(collapse = ", ")) |>
					   ungroup()
write_tsv(cs_export, "C:/.../data/13-02-26__selected_cmpnds.tab")