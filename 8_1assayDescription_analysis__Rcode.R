library(tidyverse)

# Input
assays_labeled <- read_tsv(".../data/assays_against_selectedCL.tsv")
tid <- read_tsv(".../data/acts_against_selectedCL.tsv") |> select(assay_id, tid) |> distinct()
smo_assays <- tid |> filter(tid == "101400")

# Separate to the descriptions of assays related to cancer cell-lines and protein (SMO)
smo_assays <- assays_labeled |> inner_join(smo_assays)
cell_assays <- assays_labeled |> anti_join(smo_assays)

# Summarize by comments for SMO
smo_assays_cmnt_summ <- smo_assays |> group_by(comment, comment_type) |> summarise(n=n()) |> arrange(desc(n))
# It is clearly seen that there are two types of assays:
# one, in which the direct interaction of the SMO and compound is assessed.
# second, in which interaction of the tested compound is assessed via the resulting effect on the cell, in which the SMO is expressed.
# This is in line with the classification of the assays in ChEMBL into Binding and Functional ones.

# Summarize by the time of exposure and comment for SMO
smo_assays_exp_cmnt_summ <- smo_assays |> group_by(exposure_hrs, comment) |> summarise(n=n()) |> arrange(exposure_hrs) |> arrange(comment)
# From this summary it is seen that the Binding assays have shorter and more well-defined expossure

# Summarize by comments for Cells
cell_assays_cmnt_summ <- cell_assays |> group_by(comment, comment_type) |> summarise(n=n()) |> arrange(desc(n))
# The most of assays are cearly dedicated to the cytotoxicity assessemnt, no comments on them
# Some of the assays are dedicated to the effect of compounds on metabolism, cytotoxic effect in presence of the light or other radiance, in presence of bioactive substances.
# Some assays are dedicated to the assessemnt of the cytotoxic effect of the compounds towards the cells grafted onto the mouse, etc.

# Summarize by the time of exposure for Cells
cell_assays_exp_summ <- cell_assays |> group_by(exposure_hrs) |> summarise(n=n()) |> arrange(desc(n))
# Studies conducted in the course of 24-74 hrs constitute the majority (1497 cases),
# however, the number of studies without the indication on the time of exposure is also significant (446 cases).

# Output
write_tsv(cell_assays_cmnt_summ, ".../data/cellAssays_cmnts.tsv")