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
                 "Tehranchi/TableS2.xlsx Tehranchi/TableS2.hg19.bed", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


library(tidyverse)
library(readxl)
