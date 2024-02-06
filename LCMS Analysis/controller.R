library(ggplot2)
library(dplyr)
library(gtools)
library(fANCOVA)
library(pbapply)
library(readxl)
library(data.table)
library(stringr)
library(gt)
library(rlist)
library(purrr)
library(cowplot)
library(scales)
library(tidyverse)
library(ggtext)

# set the appropriate working directory
setwd(".")

source("common/DataExtraction.R")
source("common/LOESS_Functions_simple.R")
source("common/getWTdf.R")
source("common/Comparison_Functions.R")
source("common/plotFunctions.R")
source("common/Preprocessor.R")

# USER INPUTS

# Input File paths
NEGATIVE_RESULTS_FILE = ""
NEGATIVE_RUN_ORDER_FILE = ""
POSITIVE_RESULTS_FILE = ""
POSITIVE_RUN_ORDER_FILE = ""

# file to hold intermediate data after preprocessing, before rendering charts
OUTPUT_DATA_FILE = "data.rda"
# CSV file with intermediate data from RMD rendering to be used for further analysis
OUTPUT_CSV_FILE = "data.csv"
# file for rendered html output
OUTPUT_HTML_FILE = "latest.html"

# Results file column indexes
NEGATIVE_RESULT_FIRST_COLUMN = 28
NEGATIVE_RESULT_LAST_COLUMN = 135
POSITIVE_RESULT_FIRST_COLUMN = 28
POSITIVE_RESULT_LAST_COLUMN =132

# import csv files
neg_data <- read.csv(NEGATIVE_RESULTS_FILE, sep=",", dec=",")
neg_run_order <- read.csv(NEGATIVE_RUN_ORDER_FILE, sep=",", skip=1)
pos_data <- read.csv(POSITIVE_RESULTS_FILE, sep=",", dec=",")
pos_run_order <- read.csv(POSITIVE_RUN_ORDER_FILE, sep=",", skip=1)

# correct name notation of strains
create_corrected <- function(df, runorder, first, last){
  strains <- detect_strains(runorder)
  data_cropped <- df[1: last]
  data_cleaned <- clean_data(strains, data_cropped, runorder, first)
  run_order_cleaned <- data_cleaned[[2]]
  data_cleaned <- data_cleaned[[1]]
  last <- ncol(data_cleaned)
  data_cleaned[first:last] <- lapply(data_cleaned[first:last], function(x) as.numeric(as.character(x)))
  
  names = data_cleaned$Name
  count = 1
  dup = 1
  for (i in 1:length(names)){
    if(names[[i]]==""){
      names[[i]] = "Var"
    }
    names[[i]] = paste(names[[i]], ".", count, sep="")
    count = count + 1
  }
  data_cleaned$Name = names
  
  return (list(data_cleaned, strains, run_order_cleaned))
}

# correct the negative and positive data column names
data <- create_corrected(neg_data, neg_run_order, NEGATIVE_RESULT_FIRST_COLUMN, NEGATIVE_RESULT_LAST_COLUMN)
neg_data_cleaned <- data[[1]]
strains <- data[[2]]
neg_run_order <- data[[3]]
NEGATIVE_RESULT_LAST_COLUMN <- ncol(neg_data_cleaned)

pos_corrected <- create_corrected(pos_data, pos_run_order, POSITIVE_RESULT_FIRST_COLUMN, POSITIVE_RESULT_LAST_COLUMN)
pos_data_cleaned <- data[[1]]
pos_run_order <- data[[3]]
POSITIVE_RESULT_LAST_COLUMN <- ncol(pos_data_cleaned)

# prepare data by performing loess correction and wt fitting
source("common/PrepareData.R")

# save preprocessed data for future use
save(neg_data, neg_data.corrected, neg_data_levels, wt_fits_sigs_neg, pos_data, pos_data.corrected, pos_data_levels, wt_fits_sigs_pos, strains,
       file=OUTPUT_DATA_FILE)

# generate html output
library(rmarkdown)
render('Lines/Stage1HTMLGenerator.Rmd',
       output_file = OUTPUT_HTML_FILE,
       params = list(
         data_file=OUTPUT_DATA_FILE,
         csv_file=OUTPUT_CSV_FILE
       )
     )

