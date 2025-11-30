library(tidyverse)
library(RMariaDB)
library(DBI)
library(quanteda)
library(umap)
library(praznik)
library(mirt)
library(Rcpp)

# Connect
mysql_password = 'mysql_2024'
con <- dbConnect(
  drv = RMariaDB::MariaDB(),
  dbname = 'chembl_35',
  username = 'root',
  password = mysql_password,
  host = NULL, 
  port = 3306
)

# Select IDs of the compounds studied against the interesting targets:
# 	TID 		NAME
#	103218		caski
#	106482		c33a
#	80472		  siha
#	101400		smo
#	81020		  HepG2 		# addition by request of our colleagues
tid_selected <- tibble(target_id=c('103218', '106482', '80472', '101400', '81020'))
# Read the data on compounds, which actually made it to the training set
main_sdf <- read_file(".../data/second_set.SDF") |> as_tibble() |> separate_longer_delim(value, delim = "$$$$")
data <- main_sdf |> rowwise() |>
					mutate(compound_id = str_match(value, "END\r\n>  <id>\r\n(.*)\r\n")[2],
							  all_ids = str_match(value, regex("all_ids>\r\n(.*?)\r\n\r", dotall = TRUE))[2],
							  all_target_ids = str_match(value, regex("all_targets>\r\n(.*?)\r\n\r", dotall = TRUE))[2],
                all_assay_ids = str_match(value, regex("all_assay_ids>\r\n(.*?)\r\n\r", dotall = TRUE))[2]) |>
					ungroup()

# Filter to eliminate the compounds, which were not evaluated against the targets of interest
# It would be great to transform some columns to sets
data_selected <- data |> rowwise() |>
                  mutate(all_ids = all_ids |> str_split(", ")) |>
                  mutate(all_assay_ids = all_assay_ids |> str_split(", ")) |>
                  ungroup() |>
                  separate_longer_delim(all_target_ids, delim="\r\n") |>
                  inner_join(tid_selected, by=c("all_target_ids"="target_id"))

# Prepare the vector of assay IDs and target IDs for the records to be extracted from the database
target_all <- data_selected |> pull(all_target_ids) |> unique() |> str_c(collapse = ",")
assay_all <- data_selected |> pull(all_assay_ids) |> unlist() |> unique() |> str_c(collapse = "\",\"")
assay_all <- str_c(c("\"", assay_all, "\""), collapse = "")

# Prepare the query
activity__query <- dbSendQuery(con, str_glue('SELECT a.assay_id, a.description, a.assay_type, a.assay_test_type, a.assay_category,
                        a.assay_strain, a.assay_tissue, a.variant_id, a.tid,
                        ac.molregno, ac.standard_type, ac.standard_relation, ac.standard_value, ac.standard_units, ac.data_validity_comment,
                        ac.potential_duplicate, ac.action_type, cp.full_mwt
                      FROM assays as a JOIN
                         activities as ac JOIN
                         compound_properties as cp
                      WHERE ac.molregno = cp.molregno AND
                          a.assay_id = ac.assay_id AND
                          a.tid IN ({target_all}) AND
                          a.chembl_id IN ({assay_all})'))
activity__result <- dbFetch(activity__query)
dbClearResult(activity__query)
dbDisconnect(con)

# Stats to work with only the major units
unit_stats <- activity__result |> select(standard_units) |> group_by(standard_units) |> summarize(n=n())
# nM and ug.mL-1 should be selected

# Some basic filtering
activities <- activity__result |> filter(standard_relation == "=" &
                                          is.na(data_validity_comment) &
                                          potential_duplicate == 0 &
                                          (standard_units == 'nM' | standard_units == 'ug.mL-1') & 
                                          is.na(variant_id))

# Check the number of distinct assays in the set
assays <- activities |> select(assay_id, description) |> unique() # 2084 assays, it is possible to check them manually

# Sort the assays according to the similarity of their descriptions using extended UMAP
# Create corpus
corpus <- corpus(assays, docid_field = 'assay_id',
                      text_field = 'description')
# Generate tokens
toks <- tokens(corpus, remove_punct = TRUE, remove_number = FALSE) |>
                tokens_remove(pattern = stopwords("en")) |>
                tokens_wordstem()

# Get document-feature matrix
dfm <- dfm(toks)
dfm_d <- as.matrix(dfm) |> as.data.frame()|> mutate(across(everything(),  ~as.integer(.))) |>
                                             mutate(across(everything(),  ~if_else(.>0, 1,0)))

# Select the most distinguishing tokens
n_toks <- 200
toks_select <- MRMR( dfm_d, row.names(dfm_d) |> as.factor(), n_toks)
toks_select_names  <- names(toks_select$selection)
toks_selection <- dfm_d |> select(any_of(toks_select_names))

# Tokens to single column and back
docs_points <- toks_selection |> unite("toks", everything(), sep = "-", remove = TRUE) |>
                                  rownames_to_column(var = "assay_id") |>
                                  group_by(toks) |>
                                  mutate(assay_id = str_c(assay_id, collapse = ",")) |>
                                  slice_head(n=1) |>
                                  ungroup() |>
                                  separate_wider_delim(toks, delim = "-", names_sep = "") |>
                                  mutate(across(starts_with("toks"),  ~as.integer(.)))




# Generate background data points
sourceCpp(".../8-1_generateShortestPath.cpp")
# Generate shortest paths, hopefully
#red 
# Very basic procedure without efficiency in mind, tottally not tested
toks_list <- split(docs_points |> column_to_rownames(var = "assay_id") |> as.matrix(), seq(nrow(docs_points)))
background <- shortest_path(toks_list)
background_points <- data.frame(toks = background) |> separate_wider_delim(toks, delim = "-", names_sep = "") |>
                                mutate_all(.funs = ~as.integer(.)) |>
                                mutate(assay_id = "background") |>
                                relocate(assay_id) |>
                                distinct()

# Add
points <- bind_rows(docs_points, background_points)

# UMAP the docs
set.seed(881)
data_umap <- umap( points |> select(starts_with("toks")) |> as.data.frame(), preserve.seed = TRUE, metric = "manhattan", n_neighbors = 15  )
umap_coordinates <- data_umap$layout |> bind_cols(points |> select(assay_id)) |>
                    rename(UMAP_1 = `...1`, UMAP_2 = `...2`)
umap_plot  <- ggplot(umap_coordinates |> filter(assay_id != "background"), aes(x = UMAP_1, y = UMAP_2)) +
                               geom_point(shape=21) +
                               theme_classic() +
                               coord_fixed()
umap_plot

# UMAP the docs to 1 dimension
one_umap <- umap( points |> select(starts_with("toks")) |> as.data.frame(), preserve.seed = TRUE, metric = "manhattan", n_neighbors = 15, n_components = 1  )
umap_coordinate <- one_umap$layout |> bind_cols(points |> select(assay_id)) |>
                    rename(umap = `...1`) |>
                    separate_longer_delim(assay_id, delim = ",") |>
                    inner_join(assays |> mutate(assay_id = as.character(assay_id))) |>
                    arrange(umap)

# Export the results to study them
write_tsv(activities |> select(-description), ".../data/acts_against_selectedCL.tsv")

write_tsv(umap_coordinate, ".../data/assays_against_selectedCL.tsv")
