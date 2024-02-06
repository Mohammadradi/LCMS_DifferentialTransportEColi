#Functions for grabbing WT data

# take care of proper index for strain
get_strain <- function(x){
  #t <- strsplit(x, "DR")[[1]][2]
  #q <- unlist(strsplit(t, "_"))
  #s <- q[STRAIN_INDEX]
  pattern <- "\\d+_(\\w+)\\_(\\d+)"
  s<-str_extract(x, pattern=pattern,1)
  return(s)
}

# take care of proper index for time
get_time <- function(x)  {
  #t <- strsplit(x, "DR")[[1]][2]
  #q <- unlist(strsplit(t, "_"))
  #s <- q[TIME_INDEX]
  pattern <- "\\d+_(\\w+)\\_(\\d+)"
  s<- str_extract(x, pattern=pattern,2)
  return(s)
}



get_wt <- function(x, normData){
  df <- make_compound_df(x, normData)
  df_wt <- df[df$Strain =="WT", c(2, 4)]
  df_wt$Time <- as.numeric(as.character(df_wt$Time))
  df_wt <- na.omit(df_wt)
  df_wt[df_wt$Area == 0,1] <- 0.1
  return(df_wt)
}

gen_wt_df_list <- function(normData){
  print("Generating wt df list")
  wt_dfs <- pbapply::pblapply(c(2:ncol(normData)), get_wt, normData = normData)
}


lin.fit <- function(df_wt){
  fitLin <- lm(Area ~ Time, data=df_wt)
}

exp.fit <- function(df_wt){
  fitExp <- lm(log(Area) ~ Time, data=df_wt)
}

get_adj.r<- function(fit){
  ifelse(is.nan(summary(fit)$adj.r.squared), 0, summary(fit)$adj.r.squared)
}

getPvalue <- function(fit){
  summary(fit)$coefficients[2, 4]
}

getSlope <- function(fit){
  summary(fit)$coefficients[2, 1]
}

get_mean_at_t <- function(df_wt, t){
  mean(df_wt[df_wt$Time == t, 1])
}

uptake.nofit <- function(df_wt){
  decreasing <- get_mean_at_t(df_wt, 0) >= get_mean_at_t(df_wt, 5) & get_mean_at_t(df_wt,5) >= get_mean_at_t(df_wt,15)&get_mean_at_t(df_wt,15) >= get_mean_at_t(df_wt,30)
  five_fold<- mean(df_wt[df_wt$Time == 15,1])/mean(df_wt[df_wt$Time == 0,1])<0.3 & mean(df_wt[df_wt$Time == 30,1])/mean(df_wt[df_wt$Time == 0,1])<0.2
  allPosT0 <- !any((df_wt[df_wt$Time == 0,1]) == 0.1)
  cv <- sd(df_wt[df_wt$Time==0, 1])/mean(df_wt[df_wt$Time==0,1])
  cv50 <- cv< 0.5
  
  return(allPosT0&decreasing&cv50&five_fold)
}

secretion.nofit <- function(df_wt){
  increasing <- get_mean_at_t(df_wt,0)<=get_mean_at_t(df_wt,5)& get_mean_at_t(df_wt,5) <= get_mean_at_t(df_wt,15)& get_mean_at_t(df_wt,15) <= get_mean_at_t(df_wt,30)
  five_fold <- mean(df_wt[df_wt$Time == 30, 1])/mean(df_wt[df_wt$Time == 0, 1]) > 5
  cv <- sd(df_wt[df_wt$Time == 30, 1])/mean(df_wt[df_wt$Time == 30, 1])
  cv50 <- cv < 0.5
  return(increasing&cv50&five_fold)
}

testFits <- function(df_wt){
    exp_fit <- exp.fit(df_wt)
    lin_fit <- lin.fit(df_wt)
    
    s <- getSlope(lin_fit)
    
    
    if(max(get_adj.r(lin_fit), get_adj.r(exp_fit))>0.7){
      fit <- ifelse(get_adj.r(lin_fit) > get_adj.r(exp_fit),"Linear", "Exponential")
      slope <- ifelse(s > 0, "Exported", "Imported")
    }
    else if(secretion.nofit(df_wt)|uptake.nofit(df_wt)){
      fit <- "Unfitted"
      slope <- ifelse(s > 0, "Exported", "Imported")
    }
    else{
      fit <- "No Transport"
      slope <- "ND"
  }
  return(c(fit, slope))
}

#Problem with a weird missing dataframe which causes an error. 

get_import_export_df <- function(normData){
  wtdfs <- gen_wt_df_list(normData)
  wtfits <- pblapply(wtdfs, testFits)
  wtFitDF <- as.data.frame(do.call(cbind, wtfits))
  n <- data.frame("Name" = c("Fit", "Direction"))
  wtFitDF <- cbind(n, wtFitDF)
  names(wtFitDF) <- names(normData)
  
  dfT <- as.data.frame(t(wtFitDF[,-1]))
  dfT$Compound <- names(wtFitDF)[-1]
  dfT <- dfT[,c(3,1,2)]
  names(dfT) <- c("Compound", "Fit", "Direction")
  dfT$index <- as.character(1:nrow(dfT))
  dfT <- mutate(dfT, n_compound = paste(index, Compound, sep="_"))
  return(dfT)
}


