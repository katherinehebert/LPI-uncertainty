# Script to demonstrate the need for a multispecies index of biodiversity change
# using a simple ecological example

library(tidyverse)
library(hrbrthemes)
library(synchrony)

# 1. explore data to find a good example ----

# read LPD dataset
lpd <- read_csv("data_raw/LPR2020data_public.csv")

# subset to example system
ex_sub <- lpd %>% 
  filter(Location == "Kluane Lake Research Station, Yukon, Canada") 

# clean up and convert to long format
ex <- ex_sub %>%
  pivot_longer(cols = "1950":ncol(.), values_to = "size", names_to = "year") %>%
  mutate_at(vars(year), as.integer) %>%
  mutate_at(vars(size), as.numeric) %>%
  # remove rows with NA
  drop_na(size)

# create clean wide version
ex_wide <- ex_sub[,30:ncol(ex_sub)] %>% 
  apply(1:2, as.numeric) %>%
  t()
colnames(ex_wide) = ex_sub$Binomial

# plot time series ----
ggplot(ex) +
  geom_line(aes(x = year, y = size, col = Binomial), lwd = 1.1) +
  scale_y_log10() +  
  theme_ipsum_rc() +
  labs(y = "Population size (log10)",
       title = "Interacting populations",
       subtitle = "Coyote, Snowshoe hare, and Canadian lynx")


# calculate synchrony ----

trio_sync <- community.sync(ex_wide, nrands = 100, use = "complete.obs")
