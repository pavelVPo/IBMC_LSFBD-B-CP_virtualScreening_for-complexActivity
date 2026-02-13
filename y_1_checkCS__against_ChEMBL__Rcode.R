library(RMariaDB)
library(DBI)
library(tidyverse)

### Connect to DB
mysql_password = '...'
con <- dbConnect(
  drv = RMariaDB::MariaDB(),
  dbname = 'chembl_36',
  username = 'root',
  password = mysql_password,
  host = NULL, 
  port = 3306
)

### Select all SMILES from ChEMBL v36 and calculate MNA descriptors for them
cs__query <- dbSendQuery(con, 'SELECT cs.molfile, cs.canonical_smiles, md.chembl_id FROM compound_structures cs JOIN
									molecule_dictionary md 
									WHERE cs.molregno = md.molregno AND
									molfile IS NOT NULL')
# Execute the query
cs <- dbFetch(cs__query)
dbClearResult(cs__query)
# Close the connection
dbDisconnect(con)

# Export the results to SDF
cs_export <-  cs |> mutate(id_rec = "\r\n> <id>\r\n", end_rec = "\r\n\r\n$$$$\r\n") |>
						select(molfile, id_rec, chembl_id, end_rec) |>
						unite("record", molfile:end_rec, sep = "")
cs_smiles <- cs |> select(chembl_id, canonical_smiles)
# Write the results
write_tsv(cs_smiles, "C:/.../all_chembl36.tab")
write_lines(cs_export[[1]], "C:/.../all_chembl36.SDF", sep = "")