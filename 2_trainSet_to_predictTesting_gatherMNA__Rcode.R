library(tidyverse)

# Input the records on activities
data <- read_tsv(".../data/studied_cmpnds.tsv") |>
              mutate(compound_id = as.character(compound_id))
# Input, the file with structures is divided into the two halves; since its too large as a single file
structures_raw_1 <- read_file(".../data/studied_cmpnds_SD_1.SDF") |> as_tibble()
structures_raw_2 <- read_file(".../data/studied_cmpnds_SD_2.SDF") |> as_tibble()
structures_raw_all <- bind_rows(structures_raw_1, structures_raw_2) |>
                          separate_longer_delim(value, delim = "\r\n\r\n$$$$") |>
                          separate_wider_delim(value, delim = "\r\n>  <ID>\r\n", names = c("value", "compound_id"), too_many = "merge", too_few = "align_start") |>
                          filter(!is.na(compound_id) & !str_detect(compound_id, "PASS_ERR")) |>
                          separate_wider_delim(compound_id, delim = "\r\n>  <MNA_DESCRIPTORS>\r\n", names = c("compound_id", "mna"), too_many = "merge", too_few = "align_start") |>
                          mutate(compound_id = str_trim(compound_id), mna = str_trim(mna))

# Gather structures described by the same MNA to the single record
structures_all <- structures_raw_all |> group_by(mna) |>
                                mutate(all_compound_ids = str_c(compound_id, collapse = ", ")) |>
                                slice_head(n = 1) |>
                                ungroup()

# What to do next?
# Combine structures with the records on IC50 measurements
# Use this data set to build the classifier
# Like primordial PASS, see: 
# Also it would be nice to validate the whole thing as externally as possible.
# To provide an opportunity to do that:
# - Use only the targets having at least 10 distinct compounds
# - Reserve the 50 pcts of the records for validation, it is basically 2-fold CV
# - But first, prepare the set without MNA aggregation, just for the clarity
structures_long <- structures_all |> separate_longer_delim(all_compound_ids, delim = ", ")
main_set_long <- structures_long |> inner_join(data, by = c("all_compound_ids" = "compound_id")) |>
                                          group_by(compound_id, target_id) |>
                                          mutate(n_in_pair = n()) |>
                                          ungroup() |>
                                          group_by(compound_id) |>
                                          mutate(n_for_cmpnd = n()) |>
                                          ungroup() |>
                                          group_by(target_id) |>
                                          mutate(n_for_trgt = n()) |>
                                          ungroup() |>
                                          filter(n_for_trgt > 9)
# Prepare the first set with compounds aggregated using MNA
main_set_mna <- structures_long |> inner_join(data, by = c("all_compound_ids" = "compound_id")) |>
                                      group_by(compound_id) |>
                                      mutate(all_compound_ids = str_c(all_compound_ids, collapse = ", ")) |>
                                      ungroup() |>
                                      group_by(compound_id, target_id) |> # leave the single record compound-target pair
                                      slice_head(n = 1) |>
                                      ungroup() |>
                                      group_by(target_id) |>
                                      mutate(n_for_trgt = n()) |>
                                      ungroup() |>
                                      filter(n_for_trgt > 9)
# Prepare the first subset containing 50pcts of data points for each target
first_set_mna <- main_set_mna |> group_by(target_id) |>
                                    sample_frac(size = .5) |>
                                    ungroup()
# Prepare the second subset not containing ext_testset_mna
second_set_mna <- main_set_mna |> anti_join(first_set_mna, by = c("compound_id"="compound_id", "target_id"="target_id"))

# Prepare the first subset to be exported as SDF
first_set_sdf <- first_set_mna |>  group_by(compound_id) |>
                                    mutate(all_compound_ids = all_compound_ids |> str_c(collapse = ", ") |> str_split(", ") |> unlist() |> unique() |> str_c(collapse = ", "),
                                            all_assay_ids = assay_id |> str_c(collapse = ", ") |> str_split(", ") |> unlist() |> unique() |> str_c(collapse = ", "),
                                            all_target_ids = str_c(target_id, collapse = "\r\n"),
                                            value = str_replace(value, "\r\n", "")
                                            ) |>
                                    select(value, compound_id, all_compound_ids, all_target_ids, all_assay_ids) |>
                                    slice_head(n = 1) |>
                                    ungroup() |>
                                    mutate(id_lab = "\r\n>  <id>\r\n", 
                                            all_ids_lab = "\r\n\r\n>  <all_ids>\r\n",
                                            all_targets_lab = "\r\n\r\n>  <all_targets>\r\n",
                                            all_assays_lab = "\r\n\r\n>  <all_assay_ids>\r\n",
                                            end_rec = "\r\n\r\n$$$$") |>
                                    select(value, id_lab, compound_id, all_ids_lab, all_compound_ids, all_targets_lab, all_target_ids, all_assays_lab, all_assay_ids, end_rec) |>
                                    unite("record", value:end_rec, sep = "")
# Prepare the second subset to be exported as SDF
second_set_sdf <- second_set_mna |>  group_by(compound_id) |>
                                    mutate(all_compound_ids = all_compound_ids |> str_c(collapse = ", ") |> str_split(", ") |> unlist() |> unique() |> str_c(collapse = ", "),
                                            all_assay_ids = assay_id |> str_c(collapse = ", ") |> str_split(", ") |> unlist() |> unique() |> str_c(collapse = ", "),
                                            all_target_ids = str_c(target_id, collapse = "\r\n"),
                                            value = str_replace(value, "\r\n", "")
                                            ) |>
                                    select(value, compound_id, all_compound_ids, all_target_ids, all_assay_ids) |>
                                    slice_head(n = 1) |>
                                    ungroup() |>
                                    mutate(id_lab = "\r\n>  <id>\r\n", 
                                            all_ids_lab = "\r\n\r\n>  <all_ids>\r\n",
                                            all_targets_lab = "\r\n\r\n>  <all_targets>\r\n",
                                            all_assays_lab = "\r\n\r\n>  <all_assay_ids>\r\n",
                                            end_rec = "\r\n\r\n$$$$") |>
                                    select(value, id_lab, compound_id, all_ids_lab, all_compound_ids, all_targets_lab, all_target_ids, all_assays_lab, all_assay_ids, end_rec) |>
                                    unite("record", value:end_rec, sep = "")
# Prepare the main set to be exported as SDF
main_set_sdf <- main_set_mna |> group_by(compound_id) |>
                                    mutate(all_compound_ids = all_compound_ids |> str_c(collapse = ", ") |> str_split(", ") |> unlist() |> unique() |> str_c(collapse = ", "),
                                            all_assay_ids = assay_id |> str_c(collapse = ", ") |> str_split(", ") |> unlist() |> unique() |> str_c(collapse = ", "),
                                            all_target_ids = str_c(target_id, collapse = "\r\n"),
                                            value = str_replace(value, "\r\n", "")
                                            ) |>
                                    select(value, compound_id, all_compound_ids, all_target_ids, all_assay_ids) |>
                                    slice_head(n = 1) |>
                                    ungroup() |>
                                    mutate(id_lab = "\r\n>  <id>\r\n", 
                                            all_ids_lab = "\r\n\r\n>  <all_ids>\r\n",
                                            all_targets_lab = "\r\n\r\n>  <all_targets>\r\n",
                                            all_assays_lab = "\r\n\r\n>  <all_assay_ids>\r\n",
                                            end_rec = "\r\n\r\n$$$$") |>
                                    select(value, id_lab, compound_id, all_ids_lab, all_compound_ids, all_targets_lab, all_target_ids, all_assays_lab, all_assay_ids, end_rec) |>
                                    unite("record", value:end_rec, sep = "")
# Export the first subset
write_lines(str_c("", first_set_sdf[[1]]), ".../data/first_set.SDF")
# Export the second set
write_lines(str_c("", second_set_sdf[[1]]), ".../data/second_set.SDF")
# Export the main set
write_lines(str_c("", main_set_sdf[[1]]), ".../data/main_set.SDF")