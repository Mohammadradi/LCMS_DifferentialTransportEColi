
neg_data[neg_first:neg_last] <- apply(neg_data[neg_first:neg_last], 2, FUN=function(x){as.numeric(as.character(x))})
pos_data[pos_first:pos_last] <- apply(pos_data[pos_first:pos_last], 2, FUN=function(x){as.numeric(as.character(x))})

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
neg_data.corrected <- loess_correct(neg_data, r_order = neg_run_order, first=neg_first, last=neg_last)
pos_data.corrected <- loess_correct(pos_data, r_order = pos_run_order, first=pos_first, last=pos_last)

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

