library(tidyverse)
library(tidymodels)
library(ggupset)

# Use ggupset to visualize the data, SEE: https://github.com/const-ae/ggupset
## Input
predicted_meta <- read_tsv(".../data/2F-CV_meta.tsv")
# Since some compounds could be placed into the both subsets on condition that assigned activities are different,
# this case should be proccessed here -> average the predicted values
predicted_labels <- read_tsv(".../data/2F-CV_activities.tsv") |>
							pivot_longer(cols = !compound_id) |>
							group_by(compound_id, name) |>
							mutate(value = mean(value)) |>
							slice_head(n=1) |>
							ungroup() |>
							filter(value > 0) |>
							select(-value) |>
							mutate(name = str_trim(name)) |>
							group_by(compound_id) |>
							summarize(target = list(sort(name, decreasing = TRUE)) )
true_labels <- read_tsv(".../data/studied_all__list.tsv") |>
						filter( target_id == 103218 | target_id == 106482 | target_id == 80472 | target_id == 101400 ) |>
						mutate(target_id = as.character(target_id)) |>
						mutate( target_id = case_when(
								target_id == '103218' ~ "caski",
								target_id == '106482' ~ "c33a",
								target_id == '80472' ~ "siha",
								target_id == '101400' ~ "smo",
								.default = as.character(target_id)
							) ) |>
						select(-activity) |>
						group_by(compound_id) |>
						summarize(target = list(sort(target_id, decreasing = TRUE) ))
## Visualize
true_plot <- ggplot(true_labels, aes(x = target)) +
    				geom_bar() +
    				scale_x_upset() +
    				theme_minimal()
true_plot
predicted_plot <- ggplot(predicted_labels, aes(x = target)) +
    				geom_bar() +
    				scale_x_upset() +
    				theme_minimal()
predicted_plot



## Export the results
ggsave(".../pics/UpSet_tested.png", plot = true_plot, width = 7, height = 5, dpi = 300, units = "in", bg = "white")
ggsave(".../pics/UpSet_predictedTested.png", plot = predicted_plot, width = 7, height = 5, dpi = 300, units = "in", bg = "white")