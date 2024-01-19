# For processing a single batch of either NEG or POS
# Important inputs
# SIGN : Type of data, POS or NEG
# DATA_FOLDER : The folder for saving generated files and intermediate data
# data_file : results csv file
# runorder_file: run order file
# first : 1st column in the results file with Area data
# last : last column in the results file with Area data

# Type of data (POS/NEG/ALL)
DATA_TYPE <- "ALL"

# set the working directory
setwd('D:/R_Analysis/LCMS_SerumAnalysis-main - simple')

# Master input and output data folders
BATCH_DATA_INPUT_FOLDER <- "D:/R_Analysis/LC_MS_Data/ASKA"
BATCH_DATA_OUTPUT_FOLDER <- "D:/R_Analysis/output_files/ASKA"

# batch number
BATCH = 6


# common batch data output folder
DATA_FOLDER <- paste(BATCH_DATA_OUTPUT_FOLDER,"/Batch_",BATCH, "/",sep="")

# Set the filenames
neg_data_file <- paste(BATCH_DATA_INPUT_FOLDER,"/ESI_NEG/Batch_", BATCH,"/Results.csv", sep="" )
neg_runorder_file <- paste(BATCH_DATA_INPUT_FOLDER,"/ESI_NEG/Batch_", BATCH,"/RunOrder.csv", sep="" )
pos_data_file <- paste(BATCH_DATA_INPUT_FOLDER,"/ESI_POS/Batch_", BATCH,"/Results.csv", sep="" )
pos_runorder_file <- paste(BATCH_DATA_INPUT_FOLDER,"/ESI_POS/Batch_", BATCH,"/RunOrder.csv", sep="" )

# Import data
neg_data <- read.csv(neg_data_file, sep=",", dec=",")
neg_run_order <- read.csv(neg_runorder_file, sep=",", dec=",", skip=1)
pos_data <- read.csv(pos_data_file, sep=",", dec=",")
pos_run_order <- read.csv(pos_runorder_file, sep=",", dec=",", skip=1)

# set the columns indexes, important!
neg_first <- 28
neg_last <- 137
pos_first <- 30
pos_last <- 141

# preprocess the files, making necessary corrections, ammendments and  get the strains
source("Preprocessor.KEIO.R")
# get the strains
strains <- detect_strains(neg_run_order)

# slice target columns
neg_data <- neg_data[c(1:neg_last)]
# normalize column names
results <- normalize_data(strains, neg_data, neg_run_order, neg_first)
neg_data = results[[1]]
neg_run_order = results[[2]]
neg_last = ncol(neg_data)

pos_data <- pos_data[c(1:pos_last)]
# normalize column names
results <- normalize_data(strains, pos_data, pos_run_order, pos_first)
pos_data = results[[1]]
pos_run_order = results[[2]]
pos_last = ncol(pos_data)


# prepare the data for the batch above. Runs loess correction and fitting 
source('PrepareData.R')

# save data
# filename = paste(DATA_FOLDER, "Batch", BATCH,"_data.rda", sep="")
# save(data, data.corrected, data_levels, wt_fits_sigs, BATCH, strains, SIGN,
#     file=filename)

# prepare data to be submitted for rendering
# neg_data <- NA
# neg_data.corrected <- NA
# neg_data_levels <- NA
# wt_fits_sigs_neg <- NA
# pos_data <- NA
# pos_data.corrected <- NA
# pos_data_levels <- NA
# wt_fits_sigs_pos <- NA

# if (DATA_TYPE == "NEG"){
#   neg_data <- data
#   neg_data.corrected <- data.corrected
#   neg_data_levels <- data_levels
#   wt_fits_sigs_neg <- wt_fits_sigs
# }
# 
# if (DATA_TYPE == "POS"){
#   pos_data <- data
#   pos_data.corrected <- data.corrected
#   pos_data_levels <- data_levels
#   wt_fits_sigs_pos <- wt_fits_sigs
# }

# save data
SIGN = DATA_TYPE
filename = paste(DATA_FOLDER, "Batch", BATCH,"_data.rda", sep="")
save(neg_data, neg_data.corrected, neg_data_levels, wt_fits_sigs_neg, pos_data, pos_data.corrected, pos_data_levels, wt_fits_sigs_pos, BATCH, strains, SIGN,
       file=filename)



# # generate html output
library(rmarkdown)
output_file = paste(DATA_FOLDER,"Batch", BATCH,"_data.html", sep="")
all_plots_file =  paste(DATA_FOLDER, "allplotdata_", BATCH,".csv",sep="")
render('Stage1HTMLGenerator.Rmd',
       output_file = output_file,
       params = list(
         all_plots_file=all_plots_file,
         neg_data = neg_data,
         neg_data.corrected = neg_data.corrected,
         neg_data_levels = neg_data_levels,
         wt_fits_sigs_neg = wt_fits_sigs_neg,
         pos_data = pos_data,
         pos_data.corrected = pos_data.corrected,
         pos_data_levels = pos_data_levels,
         wt_fits_sigs_pos = wt_fits_sigs_pos,
         strains=strains
       )
     )

