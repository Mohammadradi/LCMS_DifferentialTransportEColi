
plot_compound(x, strain, dat){
  compound_name <- names(dat)[x]
  
}







plotLinFit <- function(x, dat){
  df_wt <- get_wt(x, dat)
  tit <- colnames(dat)[x]
  fitLin <- lm(Area ~ Time, data=df_wt)
  df_wt$pred <- predict(fitLin)
  r2 <- summary(fitLin)$adj.r.squared
  p <- ggplot(data = df_wt, aes(Time, Area)) + geom_point() + geom_line(aes(y=pred)) + theme_bw() + ggtitle(tit) +                            # 
    annotate("text", x = 2, y = 2.2, label = paste("r2 = ", round(r2, 2)))+
    annotate("text", x=10, y = 2.2, label=paste("p = ", signif(getPvalue(fitLin), 3)))
  return(p)
}

plotexpFit <- function(x, dat){
  df_wt <- get_wt(x, dat)
  tit <- colnames(dat)[x]
  fitExp <- lm(log(Area) ~ Time, data=df_wt)
  df <- data.frame(Time = c(0:30))
  df$pred <- exp(predict(fitExp, newdata = df))
  r2 <- summary(fitExp)$adj.r.squared
  p <- ggplot(data = df_wt, aes(Time, Area)) + geom_point() + geom_smooth(data = df, aes(y=pred, x = Time), size=1) + theme_bw() + ggtitle(tit) + 
    annotate("text", x = 2, y = 2.2, label = paste("r2 = ", round(r2, 2)))+
    annotate("text", x=10, y = 2.2, label=paste("p = ", signif(getPvalue(fitExp), 3)))  
  return(p)
}

plotNoFit <- function(x, dat){
  df_wt <- get_wt(x)
  tit <- colnames(dat)[x]
  p <- ggplot(data = df_wt, aes(Time, Area)) + geom_point() + stat_smooth() + theme_bw() + ggtitle(tit)
  return(p)
}