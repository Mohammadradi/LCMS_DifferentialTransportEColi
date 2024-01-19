# For processing a single batch of either NEG or POS
# Important inputs
# DATA_TYPE : Type of data, POS or NEG
# DATA_FOLDER : The folder for saving generated files and intermediate data
# data_file : results csv file
# runorder_file: run order file
# first : 1st column in the results file with Area data
# last : last column in the results file with Area data

# Type of data (POS/NEG)
DATA_TYPE <- "POS"

# set the working directory
setwd('D:/R_Analysis/LCMS_SerumAnalysis-main - simple')

# Master input and output data folders
BATCH_DATA_INPUT_FOLDER <- "D:/R_Analysis/LC_MS_Data/KEIO"
BATCH_DATA_OUTPUT_FOLDER <- "D:/R_Analysis/output_files/KEIO"

# batch number
BATCH = 6


# common batch data output folder
DATA_FOLDER <- paste(BATCH_DATA_OUTPUT_FOLDER,"/Batch_",BATCH, "/",sep="")

# Importing data
data_file <- paste(BATCH_DATA_INPUT_FOLDER,"/ESI_", DATA_TYPE,"/Batch_", BATCH,"/Results.csv", sep="" )
runorder_file <- paste(BATCH_DATA_INPUT_FOLDER,"/ESI_", DATA_TYPE, "/Batch_", BATCH,"/RunOrder.csv", sep="" )

data <- read.csv(data_file, sep=";", dec=",", encoding="latin1")
run_order <- read.csv(runorder_file, sep=",", dec=",", skip=1)

# set the columns indexes, important!
first <- 28
last <- 157

# preprocess the files, making necessary corrections, ammendments and  get the strains
source("Preprocessor.KEIO.R")
# get the strains
strains <- detect_strains(run_order)

# slice target columns
data <- data[c(1:last)]
# normalize column names
results <- normalize_data(strains, data, run_order,first)
data = results[[1]]
run_order = results[[2]]
last = ncol(data)


# prepare the data for the batch above. Runs loess correction and fitting 
source('PrepareData_single_noloess.R')

# save data
# filename = paste(DATA_FOLDER, "Batch", BATCH,"_data.rda", sep="")
# save(data, data.corrected, data_levels, wt_fits_sigs, BATCH, strains, SIGN,
#     file=filename)

# prepare data to be submitted for rendering
neg_data <- NA
neg_data.corrected <- NA
neg_data_levels <- NA
wt_fits_sigs_neg <- NA
pos_data <- NA
pos_data.corrected <- NA
pos_data_levels <- NA
wt_fits_sigs_pos <- NA

if (DATA_TYPE == "NEG"){
  neg_data <- data
  neg_data.corrected <- data.corrected
  neg_data_levels <- data_levels
  wt_fits_sigs_neg <- wt_fits_sigs
}

if (DATA_TYPE == "POS"){
  pos_data <- data
  pos_data.corrected <- data.corrected
  pos_data_levels <- data_levels
  wt_fits_sigs_pos <- wt_fits_sigs
}

# save data
SIGN = DATA_TYPE
filename = paste(DATA_FOLDER, "Batch", BATCH,"_data.rda", sep="")
if (DATA_TYPE == "NEG"){
  save(neg_data, neg_data.corrected, neg_data_levels, wt_fits_sigs_neg, BATCH, strains, SIGN,
       file=filename)
}
if (DATA_TYPE == "POS"){
  save(pos_data, pos_data.corrected, pos_data_levels, wt_fits_sigs_pos, BATCH, strains, SIGN,
       file=filename)
}


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

