library("httr")
library("jsonlite")

# query operon information from the DOOR database
door_url <- httr::GET("http://csbl.bmb.uga.edu/DOOR/search_ajax.php",
                       query = list(keyword = "pao1",
                                    iDisplayStart = 0,
                                    iDisplayLength = 10000))
door_content <- suppressWarnings(httr::content(door_url, as = "text"))
door_json <- jsonlite::fromJSON(door_content)

# get the operons data.frame
operon_df <- door_json$aaData

# remove operons that only have one gene
operon_df <- operon_df[operon_df[,3] >1, ]

# save operons as a list with each element being a vector storing genes in
# that operon
operons <- lapply(operon_df[, 4], function(x) unlist(strsplit(x, "; ")))
names(operons) <- operon_df[, 4]

devtools::use_data(operons, overwrite = TRUE, compress = "bzip2")
