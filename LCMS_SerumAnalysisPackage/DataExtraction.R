#Extracts a clean version of the sample name
clean_names2 <- function(name){
  x1 = strsplit(name, split="[..]")[[1]][3]
  x2 = strsplit(name, split="[..]")[[1]][4]
  return(paste(x1, x2, sep="-"))
}


#Extracts a clean version of the sample name
clean_names <- function(name){
  strsplit(name, split="[..]")[[1]][3]
}

#Returns a dataframe for a single compound
areas_df <- function(n, raw_data, first, last){
  compound <- as.character(raw_data[n,2])
  areas <- as.data.frame(t(raw_data[n, c(first:last)]))
  areas$Name = rownames(areas) #as.vector(sapply(rownames(areas), clean_names2))
  #print(areas$Name)
  rownames(areas) <- c()
  areas <- areas[,c(2,1)]
  names(areas) <- c("Name", compound)
  return(areas)
}

#Returns a dataframe for the run order
run_order_df <- function(n, r_order){
  df <- data.frame(Sample = r_order$File.Name, Run = c(1:nrow(r_order)))
  names(df) <- c("Name", "Run")
  return(df)
}

#Returns a dataframe for a single compound with run order 
merged_df <- function(n, raw_data, r_order, first, last){
  merge(run_order_df(n, r_order), areas_df(n, raw_data, first, last), by.x = "Name")
}
