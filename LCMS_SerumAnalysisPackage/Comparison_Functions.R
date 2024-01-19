get_strain_df<- function(x, strain, normData){
  df <- make_compound_df(x, normData)
  df_strain <- df[df$Strain ==strain, c(2, 4)]
  df_strain$Time <- as.numeric(as.character(df_strain$Time))
  df_strain <- na.omit(df_strain)
  df_strain[df_strain$Area == 0,1] <- 0.1
  df_strain$Strain <- rep(strain, 12)
  rownames(df_strain) <- NULL
  return(df_strain)
}

get_strain_df.controls <- function(x, strain, normData){
  df_strain <- get_strain_df()
}

make_compound_df <- function(n, normData){
  strains <- sapply(as.character(normData$Name), get_strain)
  times <- sapply(as.character(normData$Name), get_time)
  df <- data.frame(Names = normData$Name, Area = normData[,n], Strain = strains, Time = times)
  return(df)
}

plotComparisonLin <- function(x, strain, normdata){
  #TODO: Fix coloring.
  df_strain <- get_strain_df(x, strain, normdata)
  df_wt <- get_strain_df(x, "WT", normdata)
  df_ser <- get_strain_df(x, "Ser", normdata)

  tit <- colnames(normdata)[x]
  fit <- lm(Area ~ Time, data=df_strain)
  fit_wt <- lm(Area ~ Time, data=df_wt)
  fit_ser <- lm(Area ~ Time, data=df_ser)
  df_ser$predSer <- predict(fit_ser)
  df_wt$predWT <- predict(fit_wt)
  df_strain$predStrain <- predict(fit)
  
  
  p <- ggplot(data = df_wt, aes(Time, Area)) + 
    geom_point() + 
    geom_point(data=df, color="red") + 
    geom_point(data=df_ser, color="yellow")+
    geom_line(aes(y=predWT)) + 
    geom_line(data=df, aes(y=predStrain, color = "red")) + 
    theme_bw() + 
    ggtitle(paste(tit, " Strain: ", strain))
  return(p)
  
}

ttestT15 <- function(x, s, normdata){
  df <- make_compound_df(x, normdata)
  dfT15 <- na.omit(df[df$Time ==15,])
  tt <- t.test(dfT15[dfT15$Strain=="WT",2], dfT15[dfT15$Strain==s, 2])
  sig <- tt$p.value < 0.05
  return(sig)
}

ttest_time <- function(df, t, s){
  df_t <- na.omit(df[df$Time == t,])
  tt <- t.test(df_t[df_t$Strain=="WT",2], df_t[df_t$Strain==s,2])
  sig <- tt$p.value < 0.01
}




two_fold_time <- function(df, t, s){
  df_t <- na.omit(df[df$Time == t,])
  wtVals <- df_t[df_t$Strain == "WT", 2]
  sVals <- df_t[df_t$Strain == s, 2]
  twofoldhiger <- mean(wtVals) > 2*mean(sVals) & max(sVals) < min(wtVals)
  twofoldlower <- 2*mean(wtVals) < mean(sVals) & max(wtVals) < min(sVals)
  return(twofoldlower|twofoldhiger)
}

two_fold_any <- function(x, s, normdata){
  df <- make_compound_df(x, normdata)
  any_twofold <- two_fold_time(df, 5, s)|two_fold_time(df, 15, s)|two_fold_time(df, 30, s)
}

get_twofold_compounds <- function(colnums, s, normdata){
  twofolds <- pbsapply(colnums, two_fold_any, s=s, normdata=normdata)
}

add_twofold_to_df <- function(strains, fits, dat) {
  for (s in strains){
  fits[,s] <- get_twofold_compounds(c(2:ncol(dat)), s, dat)
}
return(fits)
}

sig_any <- function(x, s, normdata){
  df <- make_compound_df(x, normdata)
  any_sig <- ttest_time(df, 5, s)|ttest_time(df, 15, s)|ttest_time(df, 30, s)
}

getSigCompounds <- function(colnums, s, normdata){
  sigs <- pbsapply(colnums, sig_any, s=s, normdata=normdata)
}

add_sigs_or_twofold_to_df <- function(strains, fits, dat){
  for (s in strains){
    print(paste("Calculating Sigs for ", s))
    sigs <- getSigCompounds(c(2:ncol(dat)), s, dat)
    print(paste("Calculating twofold for ", s))
    twofold <- get_twofold_compounds(c(2:ncol(dat)), s, dat)
    fits[,s] <- sigs|twofold
  }
  return(fits)
}


add_sigs_to_df <- function(strains, fits, dat){
  for (s in strains){
    fits[,s] <- getSigCompounds(c(2:ncol(dat)), s, dat)
  }
  return(fits)
}



get_sig_transported_df <- function(strains, fits, dat){
  f <- add_sigs_to_df(strains, fits, dat)
  f_transported <- f[f$Fit != "No Transport", ]
  return(f_transported)
}

remove_sigs <- function(sigs_df){
#  nsigs <- as.vector(rowSums(sigs_df[,c(6)], na.rm = TRUE))
#for(i in 1:length(nsigs)){
#  if(nsigs[i] > 2){
#    sigs_df[i,c(6)] <- rep(FALSE, 5)
#  }
#}
return(sigs_df)
}

remove_sigs_batch6 <- function(sigs_df){
  nsigs <- as.vector(rowSums(sigs_df[,c(6:9)], na.rm = TRUE))
  for(i in 1:length(nsigs)){
    if(nsigs[i] > 2){
      sigs_df[i,c(6:9)] <- rep(FALSE, 4)
    }
  }
  return(sigs_df)
}



