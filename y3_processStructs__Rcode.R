library(tidyverse)
library(reticulate)

# Get the ChEMBL RDKit-based tool to prepare structures: https://github.com/chembl/ChEMBL_Structure_Pipeline
use_python("C:/.../python.exe")
csp <- import("chembl_structure_pipeline")

### Import SDF
cs_selected   <- read_file("C:/.../data/207selected_CR_09-02-26_SD.SDF") |>
						str_trim(side = "right") |>
						str_replace_all("\r", "\n") |>
						str_replace_all("\n\n", "\n") |>
						as_tibble() |> separate_longer_delim(value, delim = "$$$$\n") |>
						filter(value != "") |>
						rowwise() |>
						mutate(
							mol = str_match(value, regex("^(.*)\n>  <compound_id>", dotall = TRUE))[2] |> csp$standardize_molblock(),
							id_rec = "\n>  <compound_id>\n",
							id = str_match(value, ">  <compound_id>.*\n(.*)\n")[2] |> str_trim(),
							end_rec = "\n\n$$$$\n"
						) |>
						ungroup()

### Export SDF to check whether prediction results do not change
cs_export <- cs_selected |> select(mol, id_rec, id, end_rec) |>
							unite("record", mol:end_rec, sep = "")
write_lines(cs_export[[1]], "C:/.../data/207_selected_processed.SDF", sep = "")

### Results of predictions are the same