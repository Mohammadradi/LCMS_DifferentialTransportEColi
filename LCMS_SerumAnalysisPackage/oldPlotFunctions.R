

plot_compound_strain <- function(x, strain, normdata){
  df <- na.omit(make_compound_df(x, normdata)[,c(2:4)])
  rownames(df) <- NULL
  wash <- mean(df$Area[df$Strain == "WashGlucose"])
  df <- df[df$Strain %in% c("WT", strain, "Ser"),]
  df$Time <- as.numeric(as.character(df$Time))
  p <- ggplot(data=df, aes(x=Time, y=Area, color=Strain)) + 
    geom_point() + 
    geom_smooth(method="lm", se=FALSE)
  
  return(p)
}





plotComparisonLin <- function(x, strain, normdata){
  #TODO: Fix coloring.
  df_strain <- get_strain_df(x, strain, normdata)
  df_wt <- get_strain_df(x, "WT", normdata)
  df_ser <- get_strain_df(x, "Ser", normdata)
  
  tit <- colnames(normdata)[x]
  fit <- lm(Area ~ Time, data=df)
  fit_wt <- lm(Area ~ Time, data=df_wt)
  fit_ser <- lm(Area ~ Time, data=df_ser)
  df_ser$predSer <- predict(fit_ser)
  df_wt$predWT <- predict(fit_wt)
  df$predStrain <- predict(fit)
  
  
  p <- ggplot(data = df_wt, aes(Time, Area)) + 
    geom_point() + 
    geom_point(data=df, color="red") + 
    geom_point(data=df_ser, color="yellow")+
    geom_line(aes(y=predWT)) + 
    geom_line(data=df, aes(y=predStrain)) + 
    theme_bw() + 
    ggtitle(paste(tit, " Strain: ", strain))
  return(p)
  
}


plotLinFit <- function(x, normdata){
  df_wt <- get_wt(x, normdata)
  tit <- colnames(normdata)[x]
  fitLin <- lm(Area ~ Time, data=df_wt)
  df_wt$pred <- predict(fitLin)
  r2 <- summary(fitLin)$adj.r.squared
  p <- ggplot(data = df_wt, aes(Time, Area)) + geom_point() + geom_line(aes(y=pred)) + theme_bw() + ggtitle(tit) +                            # 
    annotate("text", x = 2, y = 2.2, label = paste("r2 = ", round(r2, 2)))+
    annotate("text", x=10, y = 2.2, label=paste("p = ", signif(getPvalue(fitLin), 3)))
  return(p)
}

plotexpFit <- function(x, normData){
  df_wt <- get_wt(x, normData)
  tit <- colnames(normData)[x]
  fitExp <- lm(log(Area) ~ Time, data=df_wt)
  df <- data.frame(Time = c(0:30))
  df$pred <- exp(predict(fitExp, newdata = df))
  r2 <- summary(fitExp)$adj.r.squared
  p <- ggplot(data = df_wt, aes(Time, Area)) + geom_point() + geom_smooth(data = df, aes(y=pred, x = Time), size=1) + theme_bw() + ggtitle(tit) + 
    annotate("text", x = 2, y = 2.2, label = paste("r2 = ", round(r2, 2)))+
    annotate("text", x=10, y = 2.2, label=paste("p = ", signif(getPvalue(fitExp), 3)))  
  return(p)
}

plotNoFit <- function(x){
  df_wt <- get_wt(x)
  tit <- colnames(pos_norm)[x]
  p <- ggplot(data = df_wt, aes(Time, Area)) + geom_point() + stat_smooth() + theme_bw() + ggtitle(tit)
  return(p)
}

