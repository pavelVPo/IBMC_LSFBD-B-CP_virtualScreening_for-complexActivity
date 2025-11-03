library(tidyverse)
library(RMariaDB)
library(DBI)

# Connect
mysql_password = '***'
con <- dbConnect(
  drv = RMariaDB::MariaDB(),
  dbname = 'chembl_35',
  username = 'root',
  password = mysql_password,
  host = NULL, 
  port = 3306
)

# Prepare the query to extract the training data from ChEMBL
# NB, Since these data will be used to classify compounds into the
# 1) ones, expected to be tested against the particular biological targets
# 2) ones, unexpected to be tested against the particular biological targets
# In my opinion, extensive filtering is not necessary here, tha data will be used as is
# The only requrement to the measurements: IC50 should be used as an indicator of activity (lower values -> higher activity)
# The following will be extracted from ChEMBL:
# Parent structure of chemical ( MOL )
# Corresponding IDs of chemical structures
# ID of the target
# ChEMBL IDs of assays, in which the interaction in this pair ( compound-target ) was studied ( to allow to check smth fast )
biotest__query <- dbSendQuery(con, 'SELECT cs.molfile, mh.parent_molregno as compound_id, a.molregno as structure_id, ay.tid as target_id, ay.chembl_id as assay_id
											FROM compound_structures as cs JOIN
												 molecule_hierarchy as mh JOIN
												 activities as a JOIN
												 assays as ay
											WHERE a.molregno = mh.molregno AND
												  mh.parent_molregno = cs.molregno AND
												  a.assay_id = ay.assay_id AND
												  a.standard_type = "IC50"')
biotest__result <- dbFetch(biotest__query)
dbClearResult(biotest__query)
dbDisconnect(con)

# Merge the records describing the compound-target interactions for compounds having the same parent structure
data_raw <- biotest__result |> group_by(compound_id, target_id) |>
														mutate(structure_id = str_c(structure_id |> unique(), collapse = ", "),
																		assay_id = str_c(assay_id |> unique(), collapse = ", ")) |>
														slice_head(n = 1) |>
														ungroup() |>
														distinct()

# Corrections
data_raw <- read_tsv(".../data/studied_all__raw.tsv")

# Prepare to export
structures <- data_raw |> select(molfile, compound_id) |>
														distinct() |>
														mutate(id = "\r\n>  <ID>\r\n",
																		end_rec = "\r\n\r\n$$$$",
																		molfile = str_trim(molfile, side = "right")) |>			# since some molfiles somehow have additional \r\n on the right side
														select(molfile, id, compound_id, end_rec) |>
														unite("record", molfile:end_rec, sep = "")
data <- data_raw |> select(-molfile) |> distinct()

# Export the data
write_tsv(data_raw, ".../data/studied_all__raw.tsv", quote = "needed")
write_tsv(data, ".../data/studied_cmpnds.tsv", quote = "needed")
write_lines(str_c("", structures[[1]]), ".../data/studied_cmpnds.SDF")