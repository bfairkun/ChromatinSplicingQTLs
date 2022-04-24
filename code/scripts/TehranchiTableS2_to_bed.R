#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : TehranchiTableS2_to_bed
# @created     : Friday Mar 25, 2022 14:25:31 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "Tehranchi/TableS2.xlsx Tehranchi/TableS2.hg19.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


library(tidyverse)
library(readxl)

xls_f_in <- args[1]
bed_f_out <- args[2]

read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}

dat <- read_excel_allsheets(xls_f_in) %>%
    bind_rows(.id="TF")

dat %>%
    mutate(endpos = position+1) %>%
    select(`#Chr`=Chr, position, endpos, TF, ALTallele, REFallele, higherBindingAllele, pvalue) %>%
    arrange(`#Chr`, position) %>%
    write_tsv(bed_f_out, col_names=F)
