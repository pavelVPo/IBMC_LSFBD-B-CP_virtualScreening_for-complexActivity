library(tidyverse)

## Input
# Structures
cs_cr <- 		read_tsv("C:/.../tabular_structs.tab") |>
					select(molfile, id) |>
					rename(compound_id = id) |>
					select(molfile, compound_id) |>
					mutate(molfile = molfile |> str_replace_all("\r\n", "\n") |> str_replace_all("\n", "\r\n"))
cs_control <-   read_tsv("C:/.../tabular_controls.tab") |>
					select(molfile, id) |>
					rename(compound_id = id) |>
					select(molfile, compound_id) |>
					mutate(molfile = molfile |> str_replace_all("\r\n", "\n") |> str_replace_all("\n", "\r\n"))
cs_all <- bind_rows(cs_cr, cs_control)
# Predictions
data_raw <- read_tsv("C:/.../data.tsv") |>
					filter(category != "random")
max_scr <- data_raw |> pull(scr) |> max(na.rm = TRUE)
max_scrrbf <- data_raw |> pull(scrrbf) |> max(na.rm = TRUE)
min_scr <- data_raw |> pull(scr) |> min(na.rm = TRUE)
min_scrrbf <- data_raw |> pull(scrrbf) |> min(na.rm = TRUE)

# Map the quantitative values to [-1, 1], SEE: https://stackoverflow.com/questions/12931115/algorithm-to-map-an-interval-to-a-smaller-interval
data <- data_raw |> mutate(	scr_map = (scr - min_scr) * (1 + 1) / (max_scr - min_scr) - 1,
							scrrbf_map = (scrrbf - min_scrrbf) * (1 + 1) / (max_scrrbf - min_scrrbf) - 1
						 ) |>
					select(compound_id, external_id, propensity, pass_2024, main_i, main_d, local_i, scr_map, scrrbf_map, scec1, scec2, scec3,
						category, collection, availability, purity, logP, logD, logSw, pka_ma, color, stereo, state, cas, chembl_similars, n_similars) |>
					mutate()

# Select the results of prediction and code NAs as `2`
data_filled <- data |> select(compound_id, external_id, propensity, pass_2024, main_i, main_d, local_i, scr_map, scrrbf_map, scec1, scec2, scec3) |>
							mutate_at(c("propensity", "pass_2024", "main_i", "main_d", "local_i", "scr_map", "scrrbf_map", "scec1", "scec2", "scec3"), ~ if_else(is.na(.x), 2, .x)) |>
							mutate(type = "result")

# Round the results and categorize them
data_filled_rounded <- data_filled |> mutate_at(c("propensity", "pass_2024", "main_i", "main_d", "local_i", "scr_map", "scrrbf_map", "scec1", "scec2", "scec3"), ~ round(.x, digits = 1)) |>
								mutate(type = case_when(
									type != "maybe" & between(propensity, .001, 1) & between(pass_2024, .5, 1) & (between(main_i, .3, 1) | between(main_d, .3, 1) | between(local_i, .3, 1)) & (between(scr_map, .001, 1) | between(scrrbf_map, .001, 1) ) & (between(scec1, .001, 1) | between(scec2, .001, 1) | between(scec3, .001, 1)) ~ "all",
									type != "maybe" & between(propensity, .001, 1) & between(pass_2024, .5, 1) & (between(main_i, .3, 1) | between(main_d, .3, 1) | between(local_i, .3, 1)) & (between(scr_map, .001, 1) | between(scrrbf_map, .001, 1) ) ~ "minus_scec",
									type != "maybe" & between(propensity, .001, 1) & between(pass_2024, .5, 1) & (between(main_i, .3, 1) | between(main_d, .3, 1) | between(local_i, .3, 1)) & (between(scec1, .001, 1) | between(scec2, .001, 1) | between(scec3, .001, 1)) ~ "minus_quant",
									type != "maybe" & between(propensity, .001, 1) & between(pass_2024, .5, 1) & (between(scr_map, .001, 1) | between(scrrbf_map, .001, 1) ) & (between(scec1, .001, 1) | between(scec2, .001, 1) | between(scec3, .001, 1)) ~ "minus_phpt",
									type != "maybe" & between(propensity, .001, 1) & (between(main_i, .3, 1) | between(main_d, .3, 1) | between(local_i, .3, 1)) & (between(scr_map, .001, 1) | between(scrrbf_map, .001, 1) ) & (between(scec1, .001, 1) | between(scec2, .001, 1) | between(scec3, .001, 1)) ~ "minus_pass2024",
									type != "maybe" & between(pass_2024, .5, 1) & (between(main_i, .3, 1) | between(main_d, .3, 1) | between(local_i, .3, 1)) & (between(scr_map, .001, 1) | between(scrrbf_map, .001, 1) ) & (between(scec1, .001, 1) | between(scec2, .001, 1) | between(scec3, .001, 1)) ~ "minus_propensity",
									.default = type
									))

# Add this to the main data body
result <- data_raw |> select(-propensity, -pass_2024, -main_i, -main_d, -local_i, -scec1, -scec2, -scec3) |>
							left_join(data_filled_rounded, by = "compound_id", suffix=c("",".r")) |>
							select(-ends_with(".r")) |>
							select(compound_id, category, collection, propensity, pass_2024, main_i, main_d, local_i, scr, scrrbf, scr_map, scrrbf_map, scec1, scec2, scec3, type, external_id)

# Create SDF to work with
result_sdf <- result |> inner_join(cs_all, by = "compound_id") |>
						select(molfile, compound_id, type) |>
						mutate(id_fld = "\r\n>  <compound_id>\r\n", val_fld = "\r\n\r\n>  <type>\r\n", end_rec = "\r\n\r\n$$$$") |>
						select(molfile, id_fld, compound_id, val_fld, type, end_rec) |>
						unite("record", molfile:end_rec, sep = "")
result_sdf[1,1] <- result_sdf[1,1] |> str_replace("\r\n", "")
## Output
write_tsv(result, "C:/.../result_06-04-2026.tsv")
write_lines(result_sdf[[1]], "C:/.../result_06-04-2026.SDF", sep = "")