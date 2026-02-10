library(tidyverse)
library(ggridges)

### Import
smo_data_local_raw  <- read_tsv(".../smo_data_local_raw.tsv")
smo_data_d_raw		<- read_tsv(".../smo_data_d_raw.tsv")
smo_data_i_raw		<- read_tsv(".../smo_data_i_raw.tsv")
smo_data_p_raw		<- read_tsv(".../smo_data_p_raw.tsv")

### Import Pa-Pi for controls
data_controls <- read_tsv("C:/.../shh_controls.tsv") |>
							select(-compound_id) |>
							pivot_longer(cols = 2:5, names_to = "type", values_to = "drug_score")
vismodegib <- data_controls |> filter(name == "VISMODEGIB") |> select(-name) |>
								mutate(type = as_factor(type)) |>
								mutate(type = fct_relevel(type, "propensity", "activity_i", "activity_d", "activity_local_i"))


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
									type == 'propensity' ~ as.integer(3),
									type == 'activity_i' ~ as.integer(2),
									type == 'activity_d' ~ as.integer(1),
									type == 'activity_local_i' ~ as.integer(0),
									.default = NA_integer_)
								) |>
								mutate(type = as_factor(type)) |>
								mutate(type = fct_relevel(type, "propensity", "activity_i", "activity_d", "activity_local_i")) |>
								inner_join(vismodegib, by = "type")
mean_score <- smo_data_long |> pull(score) |> mean()
vismodegib_long <- smo_data_long |> select(type, order_, drug_score) |> distinct()
#Plot
ridge_plot <- ggplot(smo_data_long, aes(x = score, y = type)) +
					geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)mean(x), fill = "grey85") +
					labs(y="Type\n", x = "\nScore") +
					geom_segment(data = vismodegib_long, aes(x = drug_score, xend = drug_score, y = c(5,4,3,2), yend = c(4,3,2,1)), color = "#80ae9a", linewidth = 1) +
					theme_minimal() +
					theme(text = element_text(family = "Roboto", size = 16)) +
					theme(legend.position = "none")
ridge_plot
ggsave("C:/.../ridges_scores_all_join.png", plot = ridge_plot, scale = 1.5, width = 7, height = 5, units="in", dpi = 300)

## Work with the score >= .2
smo_data_selected <- smo_data |> filter( activity_local_i >= .2 & activity_i >= .2 & activity_d >= .2 & propensity >= .2)
smo_data_long_selected <- smo_data_selected |> pivot_longer(cols = 2:5, names_to = "type", values_to = "score") |>
								mutate(order_ = case_when(
									type == 'propensity' ~ as.integer(3),
									type == 'activity_i' ~ as.integer(2),
									type == 'activity_d' ~ as.integer(1),
									type == 'activity_local_i' ~ as.integer(0),
									.default = NA_integer_)
								) |>
								mutate(type = as_factor(type)) |>
								mutate(type = fct_relevel(type, "propensity", "activity_i", "activity_d", "activity_local_i"))
mean_score_selected <- smo_data_long_selected |> pull(score) |> mean()
ridge_plot_selected <- ggplot(smo_data_long_selected, aes(x = score, y = type)) +
					geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)mean(x), fill = "grey85") +
					geom_segment(data = vismodegib_long, aes(x = drug_score, xend = drug_score, y = c(5,4,3,2), yend = c(4,3,2,1)), color = "#80ae9a", linewidth = 1) +
					labs(y="Type\n", x = "\nScore") +
					theme_minimal() +
					theme(text = element_text(family = "Roboto", size = 16)) +
					theme(legend.position = "none")
ridge_plot_selected
ggsave("C:/.../ridges_scores_selected.png", plot = ridge_plot_selected, scale = 1.5, width = 7, height = 5, units="in", dpi = 300)

### Export selected
write_tsv(smo_data_selected, ".../smo_data_selected.tsv")