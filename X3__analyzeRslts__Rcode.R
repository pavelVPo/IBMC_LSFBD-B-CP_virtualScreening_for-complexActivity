library(tidyverse)
library(ggridges)

### Import
smo_data_local_raw  <- read_tsv("C:/.../smo_data_local_raw.tsv")
smo_data_d_raw		<- read_tsv("C:/.../smo_data_d_raw.tsv")
smo_data_i_raw		<- read_tsv("C:/.../smo_data_i_raw.tsv")
smo_data_p_raw		<- read_tsv("C:/.../smo_data_p_raw.tsv")


### Process and analyze
## Process
smo_data_local  <- smo_data_local_raw |> rename(id = `<compound_id>`, activity_local_i = target_101400)
smo_data_d 		<- smo_data_d_raw |> rename(id = `<compound_id>`, activity_d = target_101400)
smo_data_i 		<- smo_data_i_raw |> rename(id = `<compound_id>`, activity_i = target_101400)
smo_data_p 		<- smo_data_p_raw |> rename(id = `<compound_id>`, propensity = target_101400) |> mutate(order_ = 1) |>
							mutate(propensity = str_replace(propensity, ",", ".")) |>
							mutate(propensity = as.numeric(propensity))

## Join
smo_data <- smo_data_local |> inner_join(smo_data_d, by = "id") |> inner_join(smo_data_i, by = "id") |> inner_join(smo_data_p, by = "id")

## Plot the distributions using ggridges
smo_data_long <- smo_data |> pivot_longer(cols = 2:5, names_to = "type", values_to = "score") |>
								mutate(order_ = case_when(
									type == 'propensity' ~ as.integer(0),
									type == 'activity_i' ~ as.integer(1),
									type == 'activity_d' ~ as.integer(2),
									type == 'activity_local_i' ~ as.integer(3),
									.default = NA_integer_)
								) |>
								arrange(order_) |>
								mutate(type = fct_reorder(type, order_))
mean_score <- smo_data_long |> pull(score) |> mean()
ridge_plot <- ggplot(smo_data_long, aes(x = score, y = type)) +
					geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)mean(x), fill = "grey85") +
					labs(y="Type\n", x = "\nScore") +
					geom_vline(xintercept = mean_score, color = "red") +
					theme_minimal() +
					theme(text = element_text(family = "Roboto", size = 16))
ridge_plot
ggsave("C:/.../ridges_scores_all_join.png", plot = ridge_plot, scale = 1.5, width = 7, height = 5, units="in", dpi = 300)

## Work with the score >= .2
smo_data_selected <- smo_data |> filter( activity_local_i >= .2 & activity_i >= .2 & activity_d >= .2 & propensity >= .2)
smo_data_long_selected <- smo_data_3 |> pivot_longer(cols = 2:5, names_to = "type", values_to = "score") |>
								mutate(order_ = case_when(
									type == 'propensity' ~ as.integer(0),
									type == 'activity_i' ~ as.integer(1),
									type == 'activity_d' ~ as.integer(2),
									type == 'activity_local_i' ~ as.integer(3),
									.default = NA_integer_)
								) |>
								arrange(order_) |>
								mutate(type = fct_reorder(type, order_))
mean_score_selected <- smo_data_long_selected |> pull(score) |> mean()
ridge_plot_selected <- ggplot(smo_data_long_selected, aes(x = score, y = type)) +
					geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)mean(x), fill = "grey85") +
					labs(y="Type\n", x = "\nScore") +
					geom_vline(xintercept = mean_score_selected, color = "red") +
					theme_minimal() +
					theme(text = element_text(family = "Roboto", size = 16))
ridge_plot_selected
ggsave("C:/.../ridges_scores_selected.png", plot = ridge_plot_selected, scale = 1.5, width = 7, height = 5, units="in", dpi = 300)

### Export selected
write_tsv(smo_data_selected, "C:/.../smo_data_selected.tsv")