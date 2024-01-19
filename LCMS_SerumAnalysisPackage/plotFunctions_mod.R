

make_summary_df <- function(x, strain, normdata){
  df <- na.omit(make_compound_df(x, normdata)[,c(2:4)])
  df <- df[df$Strain %in% c("WT", strain, "Ser", "WashGlucose"),]
  s <- df %>%
    group_by(Strain, Time)%>%
    summarise(Avarea = mean(Area), stdevArea = sd(Area, na.rm = TRUE))
  return(s) 
}



make_summary_df_all <- function(x, normdata){
  df <- na.omit(make_compound_df(x, normdata)[,c(2:4)])
  s <- df %>%
    group_by(Strain, Time)%>%
    summarise(Avarea = mean(Area), stdevArea = sd(Area, na.rm = TRUE))
  return(s) 
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


make_graph_all <- function(n,normdata){
  x <- make_summary_df_all(n, normdata)
  target <- c("WashGlucose", "Ser", "WT", levels(x$Strain)[!levels(x$Strain) %in% c("WT", "Ser", "WashGlucose")])
#  x$Strain <- factor(x$Strain, levels=target)
  washA <- x[x$Strain=="WashGlucose", 3][[1]]
  name <- names(normdata)[[n]]
  x$Time <- as.numeric(as.character(x$Time))
  p <- ggplot(data = x[x$Strain != "WashGlucose",], aes(x=Time, y=Avarea, group=Strain, color=Strain)) +
    geom_line() + geom_point(size=5) + 
    geom_errorbar(aes(x = Time, ymin=Avarea - stdevArea, ymax=Avarea+stdevArea), width=0.3)+
    scale_color_brewer(palette = "Set1")+ 
    theme_minimal()+
    geom_hline(yintercept = washA, color="blue", linetype="dashed") + 
    ggtitle(paste(name))
  
  return(p)
}


make_graph_strain <- function(n, normdata, strain){
  x <- make_summary_df(n, strain, normdata)
  #This level rearrangement is to make sure the colors are consistent.
  target <- c("WashGlucose", "Ser", "WT", strain)
  x$Strain <- factor(x$Strain, levels=target)
  washA <- x[x$Strain=="WashGlucose", 3][[1]]
  name <- names(normdata)[[n]]
  x$Time <- as.numeric(as.character(x$Time))
  p <- ggplot(data = x[x$Strain != "WashGlucose",], aes(x=Time, y=Avarea, group=Strain, color=Strain)) +
    geom_line() + geom_point(size=5) + 
    geom_errorbar(aes(x = Time, ymin=Avarea - stdevArea, ymax=Avarea+stdevArea), width=0.3)+
    scale_color_manual(values=c("darkgoldenrod2", "darkblue", "red2"))+ 
    theme_minimal()+
    geom_hline(yintercept = washA, color="blue", linetype="dashed") + 
    ggtitle(paste(n, ".", name))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
  
}

get_sigs <- function(strain, sig_data){
  transported <- which(sig_data$Fit != "No Transport")
  sig <- which(sig_data[,strain] == TRUE)
  return(intersect(sig, transported))
}



plot_strain <- function(strain, normdata, sig_data){
  compound_list <- get_sigs(strain, sig_data)
  #+ 1 to correct for "Name" column in normdata
  plots <- lapply(compound_list+1, make_publication_graph, strain=strain, normdata=normdata)
}

# plot detected compounds
plot_strain_detected <- function(strain, normdata, sig_data){
  compound_list <- get_sigs(strain, sig_data)
  #+ 1 to correct for "Name" column in normdata
  plots <- lapply(compound_list+1, make_publication_graph, strain=strain, normdata=normdata)
}

# plot undetected/unknown compounds
plot_strain_undetected <- function(strain, normdata, sig_data){
  compound_list <- get_sigs(strain, sig_data)
  #+ 1 to correct for "Name" column in normdata
  plots <- lapply(compound_list+1, make_publication_graph, strain=strain, normdata=normdata)
}

plot_strain_bars <- function(strain, normdata, sig_data){
  compound_list <- get_sigs(strain, sig_data)
  #+ 1 to correct for "Name" column in normdata
  plots <- lapply(compound_list+1, make_publication_graph_bar, strain=strain, normdata=normdata)
}


all_transported <- function(sigs, direction){
  dir_fit <- which(sigs$Fit!="No Transport" & sigs$Direction == direction)
}

make_graph_all_list <- function(sigs, direction, normdata){
  res <- all_transported(sigs, direction)
  graphs <- lapply(res+1, make_graph_all, normdata)
}

print_all_sigs <- function(sigs, strains){
  for (s in strains){
    print(paste("Writing", s))
    pdf(paste("M:/Documents/Experimental Data/LCMS_Vers2/ESI+/Graphs/Sigs2/", s, ".pdf", sep=""))
    invisible(lapply(sigs[[s]], print))
    dev.off()
  }
}

make_publication_graph<- function(n, normdata, strain){
  t = paste(n, names(normdata)[n], sep=".")
  x <- make_summary_df(n, strain, normdata)
  #This level rearrangement is to make sure the colors are consistent.
  target <- c("WashGlucose", "Ser", "WT", strain)
  x$Strain <- factor(x$Strain, levels=target)
  washA <- x[x$Strain=="WashGlucose", 3][[1]]
  name <- names(normdata)[[n]]
  x$Time <- as.numeric(as.character(x$Time))
  scaled_Avarea <- x$Avarea/(10^6)
  p <- ggplot(data = x[x$Strain != "WashGlucose",], aes(x=Time, y=scaled_Avarea, group=Strain, color=Strain)) +
    geom_line() + geom_point(size=plot_dot_width) + 
    geom_errorbar(aes(x = Time, ymin=Avarea - stdevArea, ymax=Avarea+stdevArea), width=0.5, size=0.5)+
    scale_color_manual(values=c(serum_line_color, wt_line_color, strain_line_color))+ 
    theme_minimal()+
    #geom_hline(yintercept = washA, color="blue", linetype="dashed") + 
    ggtitle(t)+
    #ylim(0, NA)+
    labs(x = "Time (min)", y= "Peak Area (x10^6)") +
    theme(legend.position="none")+
    scale_x_continuous(breaks = c(0, 5, 15, 30))+ 
    theme_classic()+
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title=element_text(size=30, face="bold"))

  return (list(name, p, n))
}

make_publication_graph_bar<- function(n, normdata, strain){
  t = paste(n, names(normdata)[n], sep=".")
  x <- make_summary_df(n, strain, normdata)
  #This level rearrangement is to make sure the colors are consistent.
  target <- c("WashGlucose", "Ser", "WT", strain)
  x$Strain <- factor(x$Strain, levels=target)
  washA <- x[x$Strain=="WashGlucose", 3][[1]]
  name <- names(normdata)[[n]]
  x$Time <- as.numeric(as.character(x$Time))
  df <- x[x$Strain != "WashGlucose",]
  compound<-rep(t, 12)
  
  # add compound and add to dataframe
  df2 <- data.frame(cbind(compound, df))
  # put column names
  colnames(df2) <- c('Compound', 'Strain', 'Time', 'Avarea', 'stdevArea')
  # convert zero area values to 1, to avoid division by zero
  df2$Avarea[df2$Avarea==0] <- 1
  # perform calculations
  df3 <- df2 %>% filter(Time==0|Time==30) %>% select(Compound, Strain, Time, Avarea) %>% 
    group_by(Strain) %>% mutate(dd=Avarea/Avarea[Time==0])  %>% filter(Time==30) %>% 
    mutate(Ratio=ifelse(dd<1, -log(1/dd), log(dd))) %>% select(Compound, Strain, Ratio)
  return(df3)
}

make_graph_all <- function(n,normdata){
  x <- make_summary_df_all(n, normdata)
  x <- arrange(x)
  x$Strain <- as.factor(x$Strain)
  target <- c("WashGlucose", "Ser", "WT", levels(x$Strain)[!levels(x$Strain) %in% c("WT", "Ser", "WashGlucose")])
  x <- arrange(x)
  print(x)
  x$Strain <- factor(x$Strain, levels=target)
  
  washA <- x[x$Strain=="WashGlucose", 3][[1]]
  name <- names(normdata)[[n]]
  x$Time <- as.numeric(as.character(x$Time))
  p <- ggplot(data = x[x$Strain != "WashGlucose",], aes(x=Time, y=Avarea, group=Strain, color=Strain)) +
    geom_line() + geom_point(size=8) + 
    geom_errorbar(aes(x = Time, ymin=Avarea - stdevArea, ymax=Avarea+stdevArea), width=0.3)+
    scale_color_manual(values= c("red", "green", "blue", "orange", "lightblue", "yellow", "grey" ))+ 
    scale_color_manual(values=c("darkgoldenrod2", "red", "green", "blue", "orange", "lightblue", "yellow", "grey" ))+
    
    
    
    theme_minimal()+
    #geom_hline(yintercept = washA, color="blue", linetype="dashed") + 
    ggtitle(paste(name)) +
    labs(x = "Time (min)", y= "Peak Area")+
    
    scale_x_continuous(breaks = c(0, 5, 15, 30))+ 
    theme_classic()+
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title=element_text(size=18, face="bold")) +
    scale_x_continuous(breaks = c(0, 5, 15, 30))
    #theme(legend.position="none")
  
  return(p)
}



print_all_sigs_Publication <- function(sigs, strains, b){
  for (s in strains){
    print(paste("Writing", s))
    pdf(paste("M:/Documents/Experimental Data/LCMS_Vers2/ESI+/Graphs/", b, "/", s, ".pdf", sep=""))
    invisible(lapply(sigs[[s]], print))
    dev.off()
  }
}

print_all_sigs_Publication_neg <- function(sigs, strains, b){
  for (s in strains){
    print(paste("Writing", s))
    pdf(paste("M:/Documents/Experimental Data/LCMS_Vers2/ESI-/Graphs/", b, "/", s, ".pdf", sep=""))
    invisible(lapply(sigs[[s]], print))
    dev.off()
  }
}



make_publication_graph_WT<- function(n, normdata){
  x <- make_summary_df(n, "WT", normdata)
  #This level rearrangement is to make sure the colors are consistent.
  target <- c("WashGlucose", "Ser", "WT")
  t = names(normdata)[n]
  washA <- x[x$Strain=="WashGlucose", 3][[1]]
  name <- names(normdata)[[n]]
  x$Time <- as.numeric(as.character(x$Time))
  p <- ggplot(data = x[x$Strain != "WashGlucose",], aes(x=Time, y=Avarea, group=Strain, color=Strain)) +
    geom_line() + geom_point(size=8) + 
    geom_errorbar(aes(x = Time, ymin=Avarea - stdevArea, ymax=Avarea+stdevArea), width=0.5, size=0.5)+
    scale_color_manual(values=c("darkgoldenrod2", "darkblue", "red2"))+ 
    theme_minimal()+
    #geom_hline(yintercept = washA, color="blue", linetype="dashed") + 
    ggtitle(t)+
    #ylim(0, NA)+
    labs(x = "Time (min)", y= "Peak Area")+
    theme(legend.position="none")+
    scale_x_continuous(breaks = c(0, 5, 15, 30))+ 
    theme_classic()+
    theme(axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30),
          axis.title=element_text(size=40, face="bold")) +
    scale_x_continuous(breaks = c(0, 5, 15, 30))
  
  return(p)
}


makeWTgraphfromName <- function(compound, normdata){
  compound_n <- which(names(normdata) == compound)
  g <- make_publication_graph_WT(compound_n, normdata = normdata, t = compound)
  return(g)
}




make_wt_summary_df <- function(x, normdata){
  df <- make_summary_df(x, "WT", normdata)
  df$Time <- as.numeric(df$Time)
  df <- df[df$Strain == "WT",]
  df
  
}


make_wt_plot <- function(compound, normdata){
  x <- make_summary_df(which(names(normdata==compound)), strain = "WT", normdata)
  #This level rearrangement is to make sure the colors are consistent.
  target <- c("Ser", "WT")
  x$Strain <- factor(x$Strain, levels=target)
  washA <- x[x$Strain=="WashGlucose", 3][[1]]
  name <- names(normdata)[[n]]
  x$Time <- as.numeric(as.character(x$Time))
  p <- ggplot(data = x[x$Strain != "WashGlucose",], aes(x=Time, y=Avarea, group=Strain, color=Strain)) +
    geom_line() + geom_point(size=5) + 
    geom_errorbar(aes(x = Time, ymin=Avarea - stdevArea, ymax=Avarea+stdevArea), width=0.3)+
    scale_color_manual(values=c("darkgoldenrod2", "darkblue", "red2"))+ 
    theme_minimal()+
    #geom_hline(yintercept = washA, color="blue", linetype="dashed") + 
    ggtitle(t)+
    ylim(0, NA)+
    labs(x = "Time (min)", y= "Peak Area")+
    theme(legend.position="none")
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
  
  
}
  
  
  
  
  
  
  
  