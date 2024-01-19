# import libraries

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

source("DataExtraction.R")
source("LOESS_Functions.R")
source("getWTdf.R")
source("Comparison_Functions.R")
source("plotFunctions.R")

valid_names <- function(patterns, target_list){
  ignore_columns <- c()
  for (i in c(1:length(target_list))){
    n <- target_list[[i]]
    n <- str_replace(n, "SeurmOnly_washglucose", "Ser")
    detected <- NA
    for (pattern in patterns){
      detected <- str_extract(n, pattern)
      if(!is.na(detected)){
        break
      }
    }
    target_list[i] <- detected
  }
  return (target_list)
}

remove_na_column <- function(df){
  ignore_columns <- c()
  col_names <- names(df)
  for(i in c(1:length(col_names))){
    col_name <- col_names[[i]]
    if(is.na(col_name)){
      ignore_columns <- append(ignore_columns, i)
    }
  }
  return (df[-ignore_columns])
}

detect_strains <- function(ro){
  detected <- c()
  for(n in ro$File.Name){
    strain_detected <- str_match(n, "\\d+_([a-z]+[A-Z]?)")[[2]]
    if(!is.na(strain_detected)){
      detected <- append(detected, strain_detected)
    }
  }
  return (unique(detected))
}

add_missing_time <- function(df, filenames){
  occ <- c()
  for(i in c(1:length(filenames))){
    strain = str_extract(filenames[[i]], "\\d+_\\w+")
    if(!is.na(strain)){
      for (nn in names(df)){
        if(str_detect(nn, strain)){
          filenames[[i]] <- nn
          break
        }
      }
    }
  }
  return (filenames)
}

normalize_data <- function(strains, df, run_order, first){
  col_names <- names(df)
  patterns <- c("QC\\d+")
  run_order_patterns <- c("QC\\d+")
  target_names = c(c("WT", "Ser"), strains)
  
  for(target_name in target_names){
    pattern <- paste("\\d+_", target_name, "_\\d+", sep="")
    run_order_pattern <- paste("\\d+_", target_name, sep="")
    patterns <- append(patterns, pattern)
    run_order_patterns <- append(run_order_patterns, run_order_pattern)
  }
  
  # clean names for results
  col_names[first:ncol(df)] <- valid_names(patterns, col_names[first:ncol(df)])
  colnames(df) <- col_names
  
  # remove columns with NA in data frame
  df <- remove_na_column(df)
  
  run_order$File.Name <- valid_names(run_order_patterns, run_order$File.Name)
  run_order <- run_order %>% drop_na(File.Name)
  
  # add missing time in run order
  run_order$File.Name <- add_missing_time(df, run_order$File.Name)
  run_order <- run_order[!duplicated(run_order$File.Name),]
  return (list(df, run_order))

}

