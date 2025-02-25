# script meant to correct errors in management types
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(data.table)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(htmlwidgets)
library(rstudioapi)
library(stringr)
library(lme4)

# define computation mode, i.e. "local" or "cluster"
computation_mode <- "cluster"

# if comutations are local in rstudio, detect and set script path
# automatically using rstudioapi
if (identical(computation_mode, "local")) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}

# source functions
source("../functions.R")

# set paths
pheno_dir_path <- "../../data/phenotype_data/"

# get raw phenotype data with management type errors
df_raw <- droplevels(as.data.frame(fread(paste0(
  pheno_dir_path,
  "raw_data_phenotype_uncorrected_management_types.csv"
))))

# define environments for which managements must be modified
env_list_2_prime_ <- c("BEL", "ESP", "FRA", "CHE")
as_from_for_2_prime_ <- 2020

env_list_3_prime_ <- c("ITA", "POL")
as_from_for_3_prime_ <- 2019

# rename type 2 into 2 prime for "BEL", "ESP", "FRA", "CHE"
df_raw[df_raw$Year >= as_from_for_2_prime_ &
  df_raw$Country %in% env_list_2_prime_ &
  df_raw$Management == "2", "Management"] <- "2_prime_"

# rename type 2 into 3 prime for ITA and POL
df_raw[df_raw$Year >= as_from_for_3_prime_ &
  df_raw$Country %in% env_list_3_prime_ &
  df_raw$Management == "2", "Management"] <- "3_prime_"

# rename the rest of type 2 into 1
df_raw$Management[df_raw$Management == "2"] <- "1"

# correct type names
df_raw$Management <- str_replace(
  df_raw$Management,
  pattern = "_prime_",
  replacement = ""
)

# change colnames
colnames(df_raw)[match("Envir", colnames(df_raw))] <- "Country_Year"

# create Country_Management and new Envir columns
df_raw$Country_Management <- paste(df_raw$Country,
  df_raw$Management,
  sep = "_"
)
df_raw$Envir <- paste(df_raw$Country,
  df_raw$Year, df_raw$Management,
  sep = "_"
)

# write raw phenotype with correct management types
fwrite(df_raw, paste0(
  pheno_dir_path,
  "raw_data_phenotype_corrected_management_types.csv"
))
