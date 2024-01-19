
#Generates Loess Fit using QC Samples
loess_fit <- function(df){
  qcs <- df[grepl("QC", df$Name),]
  list2env(df, envir = .GlobalEnv)
  fit <- loess.as(qcs[,2], qcs[,3], criterion="gcv",degree = 2,plot = F,control=loess.control(surface="direct"))
  return(fit)
}

#Ensures positive area
ensure_positive_area <- function(x) {
  if_else(x<0, 0, x)
}

#Normalize values of a single compound dataframe.
normalize_vals <- function(df){
  gmed <- median(df[,3])
  fit <- loess_fit(df)
  xpred <- predict(fit, df$Run)
  norm <- ensure_positive_area(df[,3] + gmed - xpred)
  df_withNorms <- mutate(df, norms = norm, pred = xpred)
}

#Plot corrected data with a dataframe input
#To Do - show QC points. 

plot_corr <- function(df){
  names(df) <- c("Name", "Run", "Area", "NormArea", "Predicted")
  p <- ggplot(df, aes(x=Run)) + 
    geom_point(aes(y=Area, color="Original"))+
    geom_point(aes(y=NormArea, color="Normalized")) + 
    geom_line(aes(y=Predicted, color="Predicted"))
  print(p)
}


plot_corr_n <- function(n, raw_data, r_order, first, last){
  df <- merged_df(n, raw_data, r_order, first, last)
  dfNorm <- normalize_vals(df)
  plot <- plot_corr(dfNorm)
}

normalize_n <- function(n, raw_data, r_order, first, last){
  df <- merged_df(n, raw_data, r_order, first, last)
  dfNorm <- normalize_vals(df)
  compound <- names(dfNorm[3])
  dfNormVals <- dfNorm[,c(1,4)]
  names(dfNormVals) <- c("Name", compound)
  return(dfNormVals)
}

#loess_correct takes input of raw data, run order dataframe and inidces of where the sample data starts and ends (first, last)
# and returns a Loesscorrected dataframe as well as printing an example graph (default of col. 100)

loess_correct <- function(raw_data, r_order, first, last, n=100){
  # progress bar
  pb <- txtProgressBar(min = 0, max = nrow(raw_data), style = 3, width = 50, char = "=")
  #cleaned_names <- names(raw_data)[c(first:last)] # sapply(names(raw_data)[c(first:last)], clean_names2)
  pos_norm <- data.frame(Name = names(raw_data)[c(first:last)])
  i <- 0
  for(compound in c(1:nrow(raw_data))){
    i <- i + 1
    setTxtProgressBar(pb, i)
    normDF <- suppressWarnings(normalize_n(compound, raw_data, r_order, first, last))
    pos_norm <- merge(pos_norm, normDF, by="Name")
  }
  #plot_corr_n(n, raw_data, r_order, first, last)
  return(pos_norm)
}
