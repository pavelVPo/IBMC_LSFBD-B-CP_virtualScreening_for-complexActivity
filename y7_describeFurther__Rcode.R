library(tidyverse)
library(ggridges)

### Import
rslt_d        <- read_tsv(".../24-02-26__ordered_mainD_physchem.tab")
rslt_i        <- read_tsv(".../24-02-26__ordered_mainI_physchem.tab")
rslt_li       <- read_tsv(".../24-02-26__ordered_localI_physchem.tab")
rslt    <- rslt_d |> inner_join(rslt_i) |> inner_join(rslt_li)
liquids <- rslt |> filter(state == "Liquid") 

### Count ideal physchem
ideal_maind <- rslt_d |> filter(1 < logP & logP < 3 &
							 1 < logD & logD < 3 &
							 -4 < logSw & logSw < -1) |>
					  arrange(desc(propensity)) |>
					  arrange(desc(main_d))
ideal_maini <- rslt_i |> filter(1 < logP & logP < 3 &
							 1 < logD & logD < 3 &
							 -4 < logSw & logSw < -1) |>
					  arrange(desc(propensity)) |>
					  arrange(desc(main_i))
ideal_locali <- rslt_li |> filter(1 < logP & logP < 3 &
							 1 < logD & logD < 3 &
							 -4 < logSw & logSw < -1) |>
					  arrange(desc(propensity)) |>
					  arrange(desc(local_i))


### Plot Pa-Pi by color
color_papi <- rslt |> select(compound_id, color, main_d, main_i, local_i)
plot_d <- ggplot(color_papi, aes(main_d)) +
			geom_histogram(aes(fill = color), color = "black") +
			facet_grid(rows = vars(color)) +
			scale_fill_manual(values = c("brown", "grey80", "orange", "violet", "yellow")) + 
			theme_classic() +
			theme(legend.position = "none")
plot_d
plot_i <- ggplot(color_papi, aes(main_i)) +
			geom_histogram(aes(fill = color), color = "black") +
			facet_grid(rows = vars(color)) +
			scale_fill_manual(values = c("brown", "grey80", "orange", "violet", "yellow")) + 
			theme_classic() +
			theme(legend.position = "none")
plot_i
plot_li <- ggplot(color_papi, aes(local_i)) +
			geom_histogram(aes(fill = color), color = "black") +
			facet_grid(rows = vars(color)) +
			scale_fill_manual(values = c("brown", "grey80", "orange", "violet", "yellow")) + 
			theme_classic() +
			theme(legend.position = "none")
plot_li

plot_ridges_d <-  ggplot(color_papi, aes(x = main_d, y = color)) +
				geom_density_ridges(aes(fill = color), quantile_lines=TRUE, quantile_fun=function(x,...)mean(x)) +
				scale_fill_manual(values = c("brown", "grey80", "orange",  "yellow", "black")) +
				theme_classic() +
				theme(legend.position = "none")
plot_ridges_d
plot_ridges_i <-  ggplot(color_papi, aes(x = main_i, y = color)) +
				geom_density_ridges(aes(fill = color), quantile_lines=TRUE, quantile_fun=function(x,...)mean(x)) +
				scale_fill_manual(values = c("brown", "grey80", "orange",  "yellow", "black")) +
				theme_classic() +
				theme(legend.position = "none")
plot_ridges_i
plot_ridges_li <-  ggplot(color_papi, aes(x = local_i, y = color)) +
				geom_density_ridges(aes(fill = color), quantile_lines=TRUE, quantile_fun=function(x,...)mean(x)) +
				scale_fill_manual(values = c("brown", "grey80", "orange",  "yellow", "black")) +
				theme_classic() +
				theme(legend.position = "none")
plot_ridges_li

ggsave("C:/Users/XPS/Documents/IBMC/2025/VVP/cervicalCancer_SHH_MMP/up_to_24-02-26/slides/pics/24-02-26__papi-d_color.png",
			plot = plot_ridges_d, width = 9, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("C:/Users/XPS/Documents/IBMC/2025/VVP/cervicalCancer_SHH_MMP/up_to_24-02-26/slides/pics/24-02-26__papi-i_color.png",
			plot = plot_ridges_i, width = 9, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("C:/Users/XPS/Documents/IBMC/2025/VVP/cervicalCancer_SHH_MMP/up_to_24-02-26/slides/pics/24-02-26__papi-li_color.png",
			plot = plot_ridges_li, width = 9, height = 7, units = "in", dpi = 300, bg = "white")