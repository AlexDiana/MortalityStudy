library(here)
library(readxl)
library(tidyverse)

# Specify the path to your Excel file
file_path <- here("Data","UKdata.csv")

# Read the Excel file
data <- read_csv(file_path) 

years <- unique(data$Year)
ages <- unique(data$Age)
groups <- c("Male","Female")

X <- length(ages)
Y <- length(years)
P <- length(groups)

d <- array(NA, dim = c(X, Y, P), dimnames = list(ages, years, groups))
E <- array(NA, dim = c(X, Y, P), dimnames = list(ages, years, groups))

for (x in 1:X) {
  for (t in 1:Y) {
    
      idx <- which(data$Age == ages[x] &
                     data$Year == years[t])
      
      if(length(idx) > 0){
        d[x,t,1] <- data$`Male Deaths`[idx]
        E[x,t,1] <- data$`Male Exposure`[idx]  
        d[x,t,2] <- data$`Female Deaths`[idx]
        E[x,t,2] <- data$`Female Exposure`[idx]  
      }
    
  }
}

save(d, E, file = "data_sex_UK.rda")
