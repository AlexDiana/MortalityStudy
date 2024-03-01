setwd("C:/Users/Alex/Google Drive/R Folder/CMI")

library(tidyverse)
library(readxl)
data1 <- read_excel("CMI data.xlsx", sheet = 1)
data2 <- read_excel("CMI data.xlsx", sheet = 2)
data3 <- read_excel("CMI data.xlsx", sheet = 3)
data4 <- read_excel("CMI data.xlsx", sheet = 4)
data5 <- read_excel("CMI data.xlsx", sheet = 5)
data6 <- read_excel("CMI data.xlsx", sheet = 6)
data7 <- read_excel("CMI data.xlsx", sheet = 7)

data1$Product <- "ACI"
data2$Product <- "ACI"
data3$Product <- "DB"
data4$Product <- "DB"
data5$Product <- "DB"
data6$Product <- "SCI"
data7$Product <- "Annuities"

data7$Year[is.na(data7$Year)] <- 2020

data <- rbind(data1, data2, data3, data4, data5, data6, data7)



# data_modified <- data %>% 
#   filter(Gender == "M",
#          Year < 2020,
#          Product != "Annuities" | ProductCategory == "Standard" ) 

data_summarised <- data %>% 
  # filter(Gender == "M",
         # Year < 2020,
         # Product != "Annuities" | ProductCategory == "Standard" ) %>% 
  group_by(Age, Year, Gender) %>% 
  summarise(Exposure = sum(LivesExposure),
            Claim = sum(IncurredClaims),
            ExpClaim = sum(ExpectedClaims)) %>% 
  mutate(Qx = Claim / Exposure,
         ExpQx = ExpClaim / Exposure) %>% 
  mutate(StdQx = sqrt(Qx * (1 - Qx) / Exposure))

save(data_summarised, file = "data_summarised.rda")

library(ggplot2)

data_summarised %>% 
  filter(Age < 70,Age > 60,
         Product != "SCI") %>%
  ggplot(aes(x = Age, 
             y = Qx, 
             ymin = Qx - 1.96 * StdQx, 
             ymax = Qx + 1.96 * StdQx,
             color = Product,
             group = Product)) + geom_line() + 
  geom_errorbar() + 
  facet_wrap(vars(Year))

data_summarised %>%
  filter(Age < 70,Age > 60) %>%
  ggplot(aes(x = Age, 
             y = Qx, 
             ymin = Qx - 1.96 * StdQx, 
             ymax = Qx + 1.96 * StdQx,
             color = Gender,
             group = Gender)) + geom_line() + 
  geom_errorbar() + 
  facet_wrap(vars(Year))
