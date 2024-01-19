count_Background <- function(rawdat){
  sum(rawdat[,"Background"][[1]])
}

remove_background <- function(rawdat){
  bg_removed <- rawdat[(rawdat[,"Background"][[1]]==FALSE),]
}


count_MS2_Pref <- function(rawdat){
  sum(rawdat[,"MS2"] == "DDA for preferred ion")
}

count_MS2_other <- function(rawdat){
  sum(rawdat[,"MS2"] == "DDA for other ion")
}

qc_areas_vect <- function(n, raw_data){
  compound <- as.character(raw_data[n,2])
  qcs <- which(grepl("QC", names(raw_data)))
  areas <- raw_data[n, qcs]
  areas_vect <- as.numeric(as.vector(areas[1,]))
  return(areas_vect)
}

get_cv <- function(v){
  cv <- sd(v)/mean(v)
}

qc_cv <- function(n, rawdat){
  cv <- get_cv(qc_areas_vect(n, rawdat))
}

get_full_matches <- function(rawdat){
  sum(rawdat[,"Annotation.Source..Predicted.Compositions"] == "Full match")
}

getChemSpiderFull <- function(rawdat){
  sum(rawdat[,"Annotation.Source..ChemSpider.Search"]== "Full match")
}



makeDescriptiveDF <- function(rawdat){
  dat <- rawdat
  n_compounds <- nrow(dat)
  ms2Pref <- count_MS2_Pref(dat)
  ms2other <- count_MS2_other(dat)
  cvs <- sapply(c(1:nrow(dat)), qc_cv, rawdat = dat)
  cvs_under_15 <- sum(cvs<0.15)
  cvs_under_30 <- sum(cvs<0.3)
  full_match <- get_full_matches(dat)
  mzcloud70 <- sum(na.omit(dat[,"mzCloud.Best.Match"] > 70))
  fullChemspider <- getChemSpiderFull
  
  df <- data.frame("Compounds" = n_compounds, 
                   "MS2" = ms2Pref + ms2other,
                   "No MS2" = n_compounds - ms2Pref-ms2other,
                   "MS2 Preffered Ion" = ms2Pref,
                   "MS2 Other Ion" = ms2other,
                   "QC CV<15" = cvs_under_15, 
                   "QC CV>30" = cvs_under_30,
                   "Full Match Predicted Composition" = full_match,
                   "MZCluod > 70" = mzcloud70)
  }

getCompoundConfidence <- function(rawdat, n){
  data <- rawdat[n, c(c(2:18), 26, 27)]
  
  names(data) <- c("Name","Level", "Formula", "PredictedComposition", "MzCloudSearch", "MzVaultSearch", "ChemspiderSearch", "MassListSearch", "Fish", "MW", "RT", "Areamax", "chemspiderResults", "MZCloudResults", "msVaultResults", "MzCloudBestMatch",
                   "MzVaultBestMatch", "MzCloudSimMatch", "MS2")
  data
}

getConfidenceDF <- function(rawdat, vect){
  df <- as.data.frame(matrix(nrow = 0, ncol = 19))
  names(df) <- c("Name", "Level", "Formula", "PredictedComposition", "MzCloudSearch", "MzVaultSearch", "ChemspiderSearch", "MassListSearch", "Fish", "MW", "RT", "Areamax", "chemspiderResults", "MZCloudResults", "msVaultResults", "MzCloudBestMatch",
                 "MzVaultBestMatch", "MzCloudSimMatch", "MS2")
  for (i in c(1:length(vect))){
    df = rbind(df, getCompoundConfidence(rawdat, vect[i]))
  }
  return(df)
}


