library(tidyverse)
library(patchwork)

### Import
# Selected cmpnds
data_main   <- read_tsv(".../13-02-26__selected_cmpnds_described.tab") |> arrange(compound_id)

### Plot some properties
data_slctd <- data_main |> filter(category == "selected")
collection_summ <- data_slctd |> group_by(collection) |> summarize(count = n())
color_summ <- data_slctd |> group_by(color) |> summarize(count = n())
stereo_summ <- data_slctd |> group_by(stereo) |> summarize(count = n())
state_summ <- data_slctd |> group_by(state) |> summarize(count = n())
availability <- data_slctd |> select(availability)
similars <- data_slctd |> select(n_similars) |> filter(n_similars > 0)
logp <- data_slctd |> select(logP)
logd <- data_slctd |> select(logD)
logsw <- data_slctd |> select(logSw)

plot_collection <- ggplot(collection_summ, aes(collection, count)) +
						geom_col(color = "black", fill = "grey80") +
						scale_y_sqrt(limits = c(0, 200)) +
						theme_minimal() +
						theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, size = 4))
plot_collection

plot_stereo <- ggplot(stereo_summ, aes(stereo, count)) +
						geom_col(color = "black", fill = "grey80") +
						scale_y_sqrt(limits = c(0, 200)) +
						theme_minimal() +
						theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, size = 4))
plot_stereo

plot_state <- ggplot(state_summ, aes(state, count)) +
						geom_col(color = "black", fill = "grey80") +
						scale_y_sqrt(limits = c(0, 200)) +
						theme_minimal() +
						theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, size = 4))
plot_state

plot_color <- ggplot(color_summ, aes(color, count, fill = color)) +
						geom_col(color = "black") +
						scale_y_continuous(limits = c(0, 200)) +
						scale_fill_manual(values = c("brown", "grey80", "orange", "violet", "yellow")) +
						theme_minimal() +
						guides(fill = "none") +
						theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, size = 4))
plot_color

plot_avail <- ggplot(availability, aes(availability)) +
						geom_histogram(color = "black", fill = "grey80") +
						scale_x_sqrt(breaks = c(1, 2, 5, 10, 25, 50, 100, 200, 390)) +
						scale_y_continuous(breaks = c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55)) +
						theme_minimal() +
						theme(axis.title = element_blank())
plot_avail

plot_similars <- ggplot(similars, aes(n_similars)) +
						geom_histogram(color = "black", fill = "grey80") +
						scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
						scale_y_sqrt(breaks = c(1, 2, 6, 45)) +
						theme_minimal() +
						theme(axis.title = element_blank())
plot_similars

plot_logp <- ggplot(logp, aes(logP)) +
				geom_histogram(color = "black", fill = "grey80") +
				scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7)) +
				scale_y_continuous(limits = c(0, 30)) +
				annotate("rect", xmin = 1, xmax = 3, ymin = 0, ymax = 25, alpha = 0.3, fill = "#228b22") +
				annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 25, alpha = 0.1, fill = "#228b22") +
				annotate("rect", xmin = 3, xmax = 5, ymin = 0, ymax = 25, alpha = 0.1, fill = "#228b22") +
				theme_minimal() +
				theme(axis.title = element_blank())
plot_logp

plot_logd <- ggplot(logd, aes(logD)) +
				geom_histogram(color = "black", fill = "grey80") +
				scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7)) +
				scale_y_continuous(limits = c(0, 30)) +
				annotate("rect", xmin = 1, xmax = 3, ymin = 0, ymax = 25, alpha = 0.3, fill = "#228b22") +
				annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 25, alpha = 0.1, fill = "#228b22") +
				annotate("rect", xmin = 3, xmax = 5, ymin = 0, ymax = 25, alpha = 0.1, fill = "#228b22") +
				theme_minimal() +
				theme(axis.title = element_blank())
plot_logd

plot_logsw <- ggplot(logsw, aes(logSw)) +
				geom_histogram(color = "black", fill = "grey80") +
				scale_x_continuous(breaks = c(-8, -7, -6, -5, -4, -3, -2, -1, 0)) +
				scale_y_continuous(limits = c(0, 30)) +
				annotate("rect", xmin = -4, xmax = -1, ymin = 0, ymax = 25, alpha = 0.3, fill = "#228b22") +
				annotate("rect", xmin = -5.7, xmax = -4, ymin = 0, ymax = 25, alpha = 0.1, fill = "#228b22") +
				theme_minimal() +
				theme(axis.title = element_blank())
plot_logsw

design <- "1234
		   5666"
plot_general <- plot_collection + plot_stereo + plot_state + plot_color + plot_similars + plot_avail + plot_layout(design = design)
plot_general

### Export plots
ggsave(".../16-02-26__gen_props.png",
			plot = plot_general, width = 9, height = 7, units = "in", dpi = 300, bg = "white")
ggsave(".../16-02-26__logp.png",
			plot = plot_logp, width = 7, height = 5, units = "in", dpi = 300, bg = "white")
ggsave(".../16-02-26__logd.png",
			plot = plot_logd, width = 7, height = 5, units = "in", dpi = 300, bg = "white")
ggsave(".../16-02-26__logsw.png",
			plot = plot_logsw, width = 7, height = 5, units = "in", dpi = 300, bg = "white")