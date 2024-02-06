# Serum Incubation Experiments Analysis

## Summary
This tool helps perform data analysis by producing line charts detailing trends of WT, Serum, QC and Strains, based on results obtained from Compound Discoverer 3.1, after conducting LCMS Analysis.

## Usage
The entry point for the analysis is the **controller.R** file

The following setings have to be provided, in the **controller.R** file

- Input files
    
    The user has to provide the following file paths required in the analysis
    - Positive and Negative results files
    - Positive and Negative run order files

    These are specified in the controller.R as;

        NEGATIVE_RESULTS_FILE = ""
        NEGATIVE_RUN_ORDER_FILE = ""

        POSITIVE_RESULTS_FILE = ""
        POSITIVE_RUN_ORDER_FILE = ""

    The input files are expected to be comma separated CSV files

- Input file parameters

    Specify the first and last columns with Strain, WT, Serum data. The 1st column is usually the 1st column after MS2. Column numbers start with 1

    Example
        
        NEGATIVE_RESULT_FIRST_COLUMN = 28
        NEGATIVE_RESULT_LAST_COLUMN = 140
        POSITIVE_RESULT_FIRST_COLUMN = 28
        POSITIVE_RESULT_LAST_COLUMN = 142

- Output files

    Specify the following output files;

        # file to hold intermediate data after preprocessing, before rendering charts
        OUTPUT_DATA_FILE = "data.rda"
        # CSV file with intermediate data from RMD rendering to be used for further analysis
        OUTPUT_CSV_FILE = "data.csv"
        # file for rendered html output
        OUTPUT_HTML_FILE = "output.html"

- Working directory

    The entry program will need to reference other packages in the common folder. The appropriate working directory should thus be set;

        # set the appropriate working directory
        setwd(".")





