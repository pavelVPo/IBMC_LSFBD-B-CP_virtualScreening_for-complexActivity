library(tidyverse)
library(arrow)
library(ggridges)

### Import the results of prediction, SMO only, and export back
smo_data_local_raw <- read_csv2_arrow("C:/.../predictions/smoi_all_CR.CSV",
								col_select = all_of(c("<compound_id>", "target_101400")))
write_tsv(smo_data_local_raw, "C:/.../smo_data_local_raw.tsv")
smo_data_d_raw <- read_csv2_arrow("C:/.../predictions/d_all_CR.CSV",
								col_select = all_of(c("<compound_id>", "target_101400")))
write_tsv(smo_data_d_raw, "C:/.../smo_data_d_raw.tsv")
smo_data_i_raw <- read_csv2_arrow("C:/.../predictions/i_all_CR.CSV",
								col_select = all_of(c("<compound_id>", "target_101400")))
write_tsv(smo_data_i_raw, "C:/.../smo_data_i_raw.tsv")
# Import larger than memeory file, much faster than using `read_csv2_arrow` with the large files
smo_data_p_raw <- open_dataset("C:/.../predictions/p_all_CR.CSV", format = "csv", delim = ";") |>
								select("<compound_id>", "target_101400") |>
								collect()
write_tsv(smo_data_p_raw, "C:/.../smo_data_p_raw.tsv")