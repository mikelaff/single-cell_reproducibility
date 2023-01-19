# curate WTK5 and WTK6 Metadata

library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)

# OUTPUT FILES #########################################################################################################
#

# INPUT FILES ##########################################################################################################
# WTK5 Mapping csv
wtk5.mapping.csv <- "/proj/steinlab/projects/IVIV_scRNA/iddrc_hco_reproducibility/WTK5Mapping.csv"
# WTK6 Mapping csv
wtk6.mapping.csv <- "/proj/steinlab/projects/IVIV_scRNA/iddrc_hco_reproducibility/WTK6Mapping.csv"

#
# GLOBALS ##############################################################################################################
#

tmp <- readRDS("/proj/steinlab/projects/IVIV_scRNA/InfoVerifiedID/RObjects/WTK5Verified.rds")

