library(tidyverse)

country <- "France"

setwd(here("Data","Countries"))

deaths <- read.table(file = paste0(country,"-Deaths.txt"), header = T)
exposures <- read.table(file = paste0(country,"-Exposures.txt"), header = T)

data <- deaths %>% 
  inner_join(exposures, by = c("Age","Year")) %>% 
  rename("Deaths" = "Total.x",
         "Exposures" = "Total.y") %>% 
  select("Age","Year","Deaths","Exposures")

years <- unique(data$Year)
ages <- unique(data$Age)

ages[ages == "110+"] <- 111
ages <- as.numeric(ages)

X <- length(ages)
Y <- length(years)

d <- array(NA, dim = c(X, Y), dimnames = list(ages, years))
E <- array(NA, dim = c(X, Y), dimnames = list(ages, years))

for (x in 1:X) {
  for (t in 1:Y) {
    
    idx <- which(data$Age == ages[x] &
                   data$Year == years[t])
    
    if(length(idx) > 0){
      d[x,t] <- data$Deaths[idx]
      E[x,t] <- data$Exposures[idx]  
    }
    
  }
}

save(d, E, file = paste0(country,"-data.rda"))
