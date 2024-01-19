# set the working directory
setwd('D:/R_Analysis/LCMS_SerumAnalysis-main - simple')

# Master input and output data folders
BATCH_DATA_INPUT_FOLDER <- "D:/R_Analysis/LC_MS_Data/KEIO"
BATCH_DATA_OUTPUT_FOLDER <- "D:/R_Analysis/Data/KEIO"

# batch number
BATCH = 4

# Type of data (POS/NEG)
SIGN = "NEG"

# common batch data output folder
DATA_FOLDER <- paste(BATCH_DATA_OUTPUT_FOLDER,"/Batch_",BATCH, "/",sep="")

# Importing data
data_file <- paste(BATCH_DATA_INPUT_FOLDER,"/ESI_NEG/Batch_", BATCH,"/batch_",BATCH,"_NEG_Results.csv", sep="" )
runorder_file <- paste(BATCH_DATA_INPUT_FOLDER,"/ESI_NEG/Batch_", BATCH,"/batch_",BATCH,"_NEG_RunOrder.csv", sep="" )

data <- read.csv(data_file, sep=";", dec=",")
run_order <- read.csv(runorder_file, sep=",", dec=",", skip=1)

# preprocess the files
source("Preprocessor.KEIO.R")

# get the strains
strains <- detect_strains(run_order)

# set the columns indexes, important!
first <- 28
last <- 140

# slice target columns
data <- data[c(1:last)]
# normalize column names
results <- normalize_data(strains, data, run_order)
data = results[[1]]
run_order = results[[2]]
last = ncol(data)

# prepare the data for the batch above
source('PrepareData_single.R')

# save data
filename = paste(DATA_FOLDER, "Batch", BATCH,"_data.rda", sep="")
save(data, data.corrected, data_levels, wt_fits_sigs, BATCH, strains, SIGN,
     file=filename)

# generate html output
library(rmarkdown)
output_file = paste(DATA_FOLDER,"Batch", BATCH,"_data.html", sep="")
all_plots_file =  paste(DATA_FOLDER, "allplotdata_", BATCH,".csv",sep="")
render('Stage1HTMLGenerator_aska.Rmd',
       output_file = output_file,
       params = list(
         filename=filename,
         all_plots_file=all_plots_file
       )
     )

