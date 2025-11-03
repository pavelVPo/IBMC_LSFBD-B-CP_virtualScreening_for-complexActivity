library(tidyverse)
library(tidymodels)

## Input
predicted_labels <- read_tsv(".../data/2F-CV_activities.tsv") |>
								pivot_longer(cols = !compound_id)
true_labels <- read_tsv(".../data/studied_all__list.tsv") |>
						filter( target_id == 103218 | target_id == 106482 | target_id == 80472 | target_id == 101400 ) |>
						mutate(target_id = as.character(target_id)) |>
						mutate( target_id = case_when(
								target_id == '103218' ~ "caski",
								target_id == '106482' ~ "c33a",
								target_id == '80472' ~ "siha",
								target_id == '101400' ~ "smo",
								.default = as.character(target_id)
							) )

## Process
data <- predicted_labels |> left_join(true_labels, by = c("compound_id"="compound_id", "name"="target_id")) |>
								mutate(activity = if_else(!is.na(activity), '1', '0')) |>
								rename(predicted_label = value, true_label = activity) |>
								select(compound_id, name, true_label, predicted_label) |>
								mutate(true_label = fct(true_label, levels = c('1','0')))

## Assess
auc <- data |> group_by(name) |>
			   mutate(n_records = sum(2-as.integer(true_label)),
			   auc = roc_auc_vec(true_label, predicted_label, estimator = "binary")) |>
			   slice_head(n = 1) |>
			   ungroup() |>
			   rename(label = name) |>
			   select(label, n_records, auc)
# The quality of prediction is high

## Export the results
write_tsv(auc, ".../data/studied_selected__auc.tsv")

# It is safe to say that giventhe high ROC AUC values, prediction of testing could be used further to improve the results of the virtual screening