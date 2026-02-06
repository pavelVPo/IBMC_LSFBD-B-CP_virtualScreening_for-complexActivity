library(tidyverse)

### Import structures downloaded from CR: https://mol.chemrar.ru/download-database
## Read
# List files
sdfs_name <- list.files(path = "C:/.../PASS/", pattern = "*.sdf", full.names = TRUE)
# Read files to tibble
for (i in seq(1:length(sdfs_name))) {
	curr_file <- read_file(sdfs_name[i]) |>
					str_replace_all("\r", "\n") |>
					str_replace_all("\n\n", "\n") |>
					str_replace_all("\n", "\r\n") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$")
	if( i == 1) {
		sdfs_raw <- curr_file
	} else {
		sdfs_raw <- bind_rows(sdfs_raw, curr_file) # generally, it is not a good practice; preallocation is much better; however, here we have only 5 objects to bind together, so, OK.
	}
}
## Process
# Get the vec of fields
fields_sample <- sdfs_raw |> slice_head(n = 1) |> map(function(x) str_extract_all(x, "<.*>") |> unlist()) |> unlist() |> unique() |>
						map_chr(function(x) str_replace(x, "<", "")) |>
						map_chr(function(x) str_replace(x, ">", ""))
fields <- sdfs_raw |> map(function(x) str_extract_all(x, "<.*>") |> unlist()) |> unlist() |> unique() |>
						map_chr(function(x) str_replace(x, "<", "")) |>
						map_chr(function(x) str_replace(x, ">", ""))
# Vector of fileds to the vector of search patterns
# Parse SDF strings to fields
# Slow, but conceptually simple 
sdfs <- sdfs_raw |> distinct() |> rowwise() |>
							mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
								id = str_match(value, ">  <IDNUMBER>.*\r\n(.*)\r\n")[2] |> str_trim(),
								availability = str_match(value, ">  <available>.*\r\n(.*)\r\n")[2] |> str_trim(),
								H_acceptor = str_match(value, ">  <H_acceptor>.*\r\n(.*)\r\n")[2] |> str_trim(),
								H_donor = str_match(value, ">  <H_donor>.*\r\n(.*)\r\n")[2] |> str_trim(),
								B_rotN = str_match(value, ">  <B_rotN>.*\r\n(.*)\r\n")[2] |> str_trim(),
								N_O = str_match(value, ">  <N_O>.*\r\n(.*)\r\n")[2] |> str_trim(),
								logP = str_match(value, ">  <logP>.*\r\n(.*)\r\n")[2] |> str_trim(),
								logD = str_match(value, ">  <logD>.*\r\n(.*)\r\n")[2] |> str_trim(),
								logSw = str_match(value, ">  <logSw>.*\r\n(.*)\r\n")[2] |> str_trim(),
								psa = str_match(value, ">  <psa>.*\r\n(.*)\r\n")[2] |> str_trim(),
								pka_ma = str_match(value, ">  <pKa_ma>.*\r\n(.*)\r\n")[2] |> str_trim(),
								pka_mb = str_match(value, ">  <pKa_mb>.*\r\n(.*)\r\n")[2] |> str_trim(),
								n_chirals = str_match(value, ">  <N_Chirals>.*\r\n(.*)\r\n")[2] |> str_trim(),
								name = str_match(value, ">  <Name>.*\r\n(.*)\r\n")[2] |> str_trim(),
								smiles = str_match(value, ">  <Smile>.*\r\n(.*)\r\n")[2] |> str_trim(),
								color = str_match(value, ">  <Color>.*\r\n(.*)\r\n")[2] |> str_trim(),
								collection = str_match(value, ">  <Collection>.*\r\n(.*)\r\n")[2] |> str_trim(),
								inchi = str_match(value, ">  <InChI>.*\r\n(.*)\r\n")[2] |> str_trim(),
								inchi_key = str_match(value, ">  <InChI Key>.*\r\n(.*)\r\n")[2] |> str_trim(),
								p_sp3 = str_match(value, ">  <PERCENTSP3>.*\r\n(.*)\r\n")[2] |> str_trim(),
								link = str_match(value, ">  <Link>.*\r\n(.*)\r\n")[2] |> str_trim(),
								mfcd_n = str_match(value, ">  <MFCDNUMBER>.*\r\n(.*)\r\n")[2] |> str_trim(),
								cl_n = str_match(value, ">  <CLNUMBER>.*\r\n(.*)\r\n")[2] |> str_trim(),
								stereo = str_match(value, ">  <STEREO>.*\r\n(.*)\r\n")[2] |> str_trim(),
								state = str_match(value, ">  <State>.*\r\n(.*)\r\n")[2] |> str_trim(),
								purity = str_match(value, ">  <Purity>.*\r\n(.*)\r\n")[2] |> str_trim(),
								salt = str_match(value, ">  <Saltdata>.*\r\n(.*)\r\n")[2] |> str_trim(),
								salt_mw = str_match(value, ">  <Saltdata_MW>.*\r\n(.*)\r\n")[2] |> str_trim(),
								cas = str_match(value, ">  <CAS>.*\r\n(.*)\r\n")[2] |> str_trim()
							) |>
							select(-value) 
sdf <- sdfs |> ungroup() |> distinct()
# Save this intermediate results
write_tsv(sdf, "C:/.../data/tabular_structs.tab", quote = "all") 						# 1 614 294
write_tsv(sdf |> select(-molfile), "C:/.../data/tabular_props.tab", quote = "all")

### Export structures
sdf_export <- sdf |> select(molfile, id) |> distinct() |> mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$",
                                      molfile = str_trim(molfile, side = "right")) |>
                              select(molfile, id_fld, id, end_rec) |>
                              unite("record", molfile:end_rec, sep = "")
# Write structures 
write_lines(sdf_export[[1]], "C:/.../data//all_CR_raw.SDF", sep = "")

### Import MNAed structures, read it chunked externally not just to save time, but also to explicitly illustrate the whole idea
sdf_mna_1 <- read_file("C:/.../data//all_CR_raw_SD.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, ">  <MNA_DESCRIPTORS>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna_2 <- read_file("C:/.../data//all_CR_raw_SD_I.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, ">  <MNA_DESCRIPTORS>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna_3 <- read_file("C:/.../data//all_CR_raw_SD_II.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, ">  <MNA_DESCRIPTORS>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna_4 <- read_file("C:/.../data//all_CR_raw_SD_III.SDF") |>
						 str_trim(side = "right") |>
						 as_tibble() |>
						 filter(!is.na(value)) |>
						 separate_longer_delim(value, delim = "$$$$") |>
						 rowwise() |>
						 mutate(molfile = str_match(value, regex(".*END", dotall = TRUE))[1],
						 		id = str_match(value, ">  <compound_id>.*\r\n(.*)\r\n")[2] |> str_trim(),
						 		mna = str_match(value, ">  <MNA_DESCRIPTORS>.*\r\n(.*)\r\n")[2] |> str_trim()) |>
						 ungroup() |>
						 filter(!is.na(mna)) |>
						 mutate(id_fld = "\r\n>  <compound_id>\r\n",  end_rec = "\r\n\r\n$$$$") |>
						 select(-value) |>
						 distinct()
sdf_mna <- bind_rows(sdf_mna_1, sdf_mna_2, sdf_mna_3, sdf_mna_4)

### Export filtered structures and properties
# Write OK structures (1614294)
sdf_export <- sdf_mna |> select(molfile, id_fld, id, end_rec) |>
                         unite("record", molfile:end_rec, sep = "")

write_lines(sdf_export[[1]], "C:/.../data//all_CR.SDF", sep = "")
