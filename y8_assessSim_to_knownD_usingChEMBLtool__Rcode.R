library(tidyverse)
library(reticulate)

# Get the ChEMBL RDKit-based tool to assess similarity: https://github.com/chembl/FPSim2
use_python("C:/.../python.exe")
fpsim2 <- import("FPSim2")

### Import SDFs
# Selected
cs_selected   <- read_file("C:/.../207selected_CR_09-02-26_SD.SDF") |>
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
# Props
props <- read_tsv("C:/.../tabular_props.tab")

### Create DB for searching
fpsim2$io$create_db_file(
    mols_source="C:/.../Smo_controls.smi",
    filename='controls__db.h5',
    mol_format=py_none,
    fp_type='Morgan',
    mol_id_prop='chembl_id'
)

### Assess similarity between the datasets using RDkit
fpdb <- fpsim2$FPSim2Engine('controls__db.h5')

### Search similars
### Search for the similar cmpnds in ChEMBL
cs_selected_sim <- cs_selected |> inner_join(props, by = "id") |> 
								  rowwise() |>
								  mutate(drug_similars_01 = tryCatch(fpdb$similarity(smiles, threshold=r_to_py(0.1),
								  			metric='tanimoto', n_workers=r_to_py(1))$tolist() |> map(function(x) paste0("CHEMBL", x[1])),
								  		 error = function(e) NA_character_) |> list() ) |>
								  mutate(n_similars_01 = drug_similars_01 |> length()) |>
								  mutate(drug_similars_03 = tryCatch(fpdb$similarity(smiles, threshold=r_to_py(0.3),
								  			metric='tanimoto', n_workers=r_to_py(1))$tolist() |> map(function(x) paste0("CHEMBL", x[1])),
								  		 error = function(e) NA_character_) |> list() ) |>
								  mutate(n_similars_03 = drug_similars_03 |> length()) |>
								  ungroup()