
data[first:last] <- apply(data[first:last], 2, FUN=function(x){as.numeric(as.character(x))})

# Compute Levels
print('Computing levels')
data_levels <- data %>% mutate(level = case_when(mzVault.Best.Match > 70 ~ 1,
                                                 mzCloud.Best.Match > 70 ~ 2, 
                                                 Annotation.Source..Predicted.Compositions == "Full match" ~ 3,
                                                 Formula != "" ~ 4,
                                                 Molecular.Weight > 0 ~ 5))

data_levels <- data_levels[,c(1, 2, last+1, c(3:last))]

# Perform LOESS correction
print('Perfoming loess correction')
data.corrected <- loess_correct(data, r_order = run_order, first=first, last=last)

# Generate fits for wt strains to determine which are being imported and exported as well as how:
print('Perfoming WT fitting')
wtFits <- get_import_export_df(data.corrected)


# Generate plots of significance testing for 
print('Perfoming significance testing for strains')

wt_fits_sigs<- add_sigs_or_twofold_to_df(strains, fits = wtFits, dat=data.corrected)
wt_fits_sigs <- remove_sigs(wt_fits_sigs)

fits_summary <- wt_fits_sigs %>%
  group_by(Direction, Fit) %>%
  summarize(n = n())

print(fits_summary)
      
# Number of hits for each strain 
transported <- wt_fits_sigs[wt_fits_sigs$Fit != "No Transport",]
for(s in strains){
  print(paste(s, ": ", sum(na.omit(transported[,s]))))
}
      
      