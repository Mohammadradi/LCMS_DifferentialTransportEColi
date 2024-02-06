# This pipeline completes full data processing and analysis for a batch of Serum LCMS Uptake experiments. 
# Import Everything: 

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

setwd("/Users/moosabonomali/Freelancer/SmartBio/DataAnalysis")
source("common/DataExtraction.R")
source("common/LOESS_Functions_simple.R")
source("common/getWTdf.R")
source("common/Comparison_Functions.R")
source("common/plotFunctions.R")
source("common/Preprocessor.R")

#Auto Extract Range Info

auto_get_range <- function(df){
  cols <- names(df)
  for (i in c(1:ncol(df))){
    if(str_detect(cols[i], "Plate\\d+_")==TRUE){
      break
    }
  }
  return (c(i, ncol(df)))
}

#Auto Extract Strains

auto_get_strains <- function(df){
  pattern = "Plate\\d+_\\d+_([A-Za-z]+)_\\d+"
  strains <- c()
  cols <- names(df)
  for (i in c(1:ncol(df))){
    strain<-str_extract(cols[i], pattern, 1)
    if(strain %in% c("WT", "Ser", "WashGlucose", NA,"Pooled") == FALSE){
      if(strain %in% strains == FALSE){
        strains[length(strains)+1] <- strain
      }
    }
  }
  return(strains)
}


# Get column first, last and strains
print('Extracting number of columns and strains')
neg_columns<-auto_get_range(neg_data)
pos_columns<-auto_get_range(pos_data)
neg_strains<-auto_get_strains(neg_data)
pos_strains<-auto_get_strains(pos_data)

neg_first <- neg_columns[1]
neg_last <-  neg_columns[2]
pos_first <- pos_columns[1]
pos_last <-  pos_columns[2]
neg_data[neg_first:neg_last] <- apply(neg_data[neg_first:neg_last], 2, FUN=function(x){as.numeric(as.character(x))})
pos_data[pos_first:pos_last] <- apply(pos_data[pos_first:pos_last], 2, FUN=function(x){as.numeric(as.character(x))})
strains <- auto_get_strains(pos_data)

# Compute Levels
print('Computing levels')
neg_data_levels <- neg_data %>% mutate(level = case_when(mzVault.Best.Match > 70 ~ 1,
                                                         mzCloud.Best.Match > 70 ~ 2, 
                                                         Annotation.Source..Predicted.Compositions == "Full match" ~ 3,
                                                         Formula != "" ~ 4,
                                                         Molecular.Weight > 0 ~ 5))

neg_data_levels <- neg_data_levels[,c(1, 2, neg_last+1, c(3:neg_last))]

pos_data_levels <- pos_data %>% mutate(level = case_when(mzVault.Best.Match > 70 ~ 1,
                                                         mzCloud.Best.Match > 70 ~ 2, 
                                                         Annotation.Source..Predicted.Compositions == "Full match" ~ 3,
                                                         Formula != "" ~ 4,
                                                         Molecular.Weight > 0 ~ 5))

pos_data_levels <- pos_data_levels[,c(1, 2, pos_last+1, c(3:pos_last))]


# Perform LOESS correction
print('Perfoming loess correction')
neg_data.corrected <- loess_correct(neg_data_cleaned, r_order = neg_run_order, first=NEGATIVE_RESULT_FIRST_COLUMN, last=NEGATIVE_RESULT_LAST_COLUMN)
pos_data.corrected <- loess_correct(pos_data_cleaned, r_order = pos_run_order, first=POSITIVE_RESULT_FIRST_COLUMN, last=POSITIVE_RESULT_LAST_COLUMN)

# Generate fits for wt strains to determine which are being imported and exported as well as how:
print('Perfoming WT fitting')
neg_wtFits <- get_import_export_df(neg_data.corrected)
pos_wtFits <- get_import_export_df(pos_data.corrected)


# Generate plots of significance testing for 
print('Perfoming significance testing for strains')
wt_fits_sigs_pos <- add_sigs_or_twofold_to_df(strains, fits = pos_wtFits, dat=pos_data.corrected)
wt_fits_sigs_pos <- remove_sigs(wt_fits_sigs_pos)

wt_fits_sigs_neg <- add_sigs_or_twofold_to_df(strains, fits = neg_wtFits, dat=neg_data.corrected)
wt_fits_sigs_neg <- remove_sigs(wt_fits_sigs_neg)

fits_summary_neg <- wt_fits_sigs_neg %>%
  group_by(Direction, Fit) %>%
  summarize(n = n())

print(fits_summary_neg)


fits_summary_pos <- wt_fits_sigs_pos %>%
  group_by(Direction, Fit) %>%
  summarize(n = n())

print(fits_summary_pos)

# Number of hits for each strain - Neg. 
transported_neg <- wt_fits_sigs_neg[wt_fits_sigs_neg$Fit != "No Transport",]
for(s in strains){
  print(paste(s, ": ", sum(na.omit(transported_neg[,s]))))
}

# Number of hits for each strain - Pos
transported_pos <- wt_fits_sigs_pos[wt_fits_sigs_pos$Fit != "No Transport",]

for(s in strains){
  print(paste(s, ": ", sum(na.omit(transported_pos[,s]))))
}

