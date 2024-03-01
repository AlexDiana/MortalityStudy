library(here)
library(readxl)
library(tidyverse)

# Specify the path to your Excel file
file_path <- here("Data","CMI Annuities PAIP all offices datasheets 2020 v01 2022-03-28 (1).xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Individual Internal") %>% 
  select(Age, `Lives exposure`, `Incurred deaths`) %>% 
  group_by(Age) %>% 
  summarise(Exposure = sum(`Lives exposure`),
            Deaths = sum(`Incurred deaths`))

qplot(data$Age, data$Deaths / data$Exposure)

file_path <- here("Data","Exposure_Deaths_UK.xlsx")

data2_exposure <- read_excel(file_path, sheet = "Exposure")
data2_death <- read_excel(file_path, sheet = "Deaths")

data2 <- full_join(data2_exposure, data2_death, by = c("Year","Age")) %>% 
  filter(Year == 2020,
         Age != "110+") %>% 
  mutate(Rate = Total.y / Total.x,
         Age = as.numeric(Age))
  
ggplot() + 
  geom_point(data = data, aes(x = Age, 
                              y = Deaths / Exposure)) +
  geom_point(data = data2, aes(x = Age, 
                               y = Rate), color = "red") + 
  coord_cartesian(xlim = c(60, 100))

