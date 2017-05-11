#!/usr/bin/env Rscript
library(docopt)
"Usage: create_files OUTPUT_DIR

-h --help    show this

This utility script converts the google spreadsheet here:

  https://docs.google.com/spreadsheets/d/1n3hXzzhrHZgClLD8P3cyIrK_6YgdjtdGlLswNAuoKSI/edit?usp=sharing

to a set of files including:

  - OUTPUT/single-cell-software-tidy.csv
  - OUTPUT/single-cell-software.json
" -> doc
opts = docopt(doc)
print(opts)

library(googlesheets)
library(readr)
library(jsonlite)
library(dplyr)
library(tidyr)
library(lubridate)
# this is the google spreadsheet ID to use
gskey = "1n3hXzzhrHZgClLD8P3cyIrK_6YgdjtdGlLswNAuoKSI"
gs_auth()


#' Create tidy sheet from the google sheet
#' @export
get_tidy_sw_list = function(gskey) {
  swsheet = gs_key(gskey) %>% 
    gs_read() %>%
    mutate(Preprint = (`Pub Date`=="PREPRINT")) %>%
    mutate(`Pub Date`=as_date(`Pub Date`)) %>%
    mutate(Preprint = ifelse(Preprint==TRUE,TRUE,NA)) %>%
    mutate(Added = as_date(Added)) %>%
    mutate(Updated = as_date(Updated)) %>% 
    mutate(DOI_url = ifelse(is.na(DOI),NA,paste0('http://dx.doi.org/',DOI)))
  gather(swsheet,key='category',value='val',-Description,-Name,-Platform,-DOI,-`Pub Date`,-Updated,-Added,-Preprint,-Code,-DOI_url,-License) %>%
    mutate(Github = grepl('github',Code)) %>%
    mutate(Bioconductor = grepl('bioconductor',Code,ignore.case = TRUE)) %>%
    mutate(CRAN = grepl('cran\\.r-project',Code)) %>%
    filter(val==TRUE) %>%
    select(-val)
}

tidysw_to_list_df <- function(tidysw) {
  catlist = split(tidysw$category,f=tidysw$Name)
  tidyswl = tidysw %>% select(-category) %>% unique()
  tidyswl[['categories']] = catlist[tidyswl$Name]
  tidyswl
}


#' write out json and csv files
#' 
#' @export
write_files = function(destdir) {
  dir.create(destdir, recursive = TRUE)
  swsheet = get_tidy_sw_list(gskey)
  write_csv(swsheet,path=file.path(destdir,'single-cell-software_tidy.csv'))
  writeLines(toJSON(tidysw_to_list_df(swsheet),pretty=TRUE),file.path(destdir,'single-cell-software.json'))
}

write_files(opts$OUTPUT_DIR)