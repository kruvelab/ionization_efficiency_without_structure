library(ggplot2)
library(cowplot)
library(randomForest)
library(caret)
library(tidyverse)
library(caTools)


setwd("~/Bc excel files/20210621")
dataset <- read_delim("Mix_M_dataset 20210621.csv", delim = ";", col_names = TRUE) 
setwd("~/Bc excel files")
sodium_adducts <- read_delim("Mix M sodium adducts no MW.csv", delim = ";", col_names = TRUE) 

dataset <- dataset%>%
  left_join(sodium_adducts) 


Unique_SMILES <- dataset %>% 
  select(SMILES)%>%
  unique()

set.seed(123)
split <- sample.split(Unique_SMILES$SMILES, SplitRatio = 0.8) 

Unique_SMILES <- Unique_SMILES %>% 
  mutate(split = split)

dataset <- dataset %>% 
  left_join(Unique_SMILES)

dataset <- dataset %>%
  filter(Compound_name != "1-Nitropyrene", Compound_name != "1-hydroxypyrene")







knn_st61_1tR_pH_2.7 <- dataset %>% 
  filter(pH == 2.7 & mode == "positive")



standards <- knn_st61_1tR_pH_2.7 %>%
  group_by(Compound_name, pH, slope,) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup()%>%
  na.omit()

assigning_closest_cal_comp_slope <- function(xRT, cal_compounds) {
  RF_cal <- cal_compounds %>% slice(which.min(abs(xRT - ret_time))) %>%select(slope)
  unlist(RF_cal)
}

knn_st61_1tR_slope_pred <- c()
for(i in 1:423){
  knn_st61_1tR_slope_pred <- c(knn_st61_1tR_slope_pred, assigning_closest_cal_comp_slope(knn_st61_1tR_pH_2.7[i,]$ret_time, standards))
}




knn_st61_1tR_pH_2.7 <- data.frame(knn_st61_1tR_pH_2.7,knn_st61_1tR_slope_pred)

knn_st61_1tR_pH_2.7 <- knn_st61_1tR_pH_2.7 %>%
  mutate(knn_st61_1tR_c_pred = Peak_Area_IC / knn_st61_1tR_slope_pred)



setwd("~/Bc excel files/20210621")
dataset_E <- read_delim("New_Mix_E_dataset 20210802.csv", delim = ";", col_names = TRUE) 

dataset_E <- dataset_E %>%
  mutate(log_slope = log10(slope)) 


sodium_adducts_E <- read_delim("Mix E sodium adducts no MW 20210621.csv", delim = ";", col_names = TRUE) 

dataset_E <- dataset_E%>%
  left_join(sodium_adducts_E) 


knn_1tR_k1_st61_pH_2.7_E <- dataset_E %>% 
  filter(pH == 2.7, mode == "positive")

knn_1tR_k1_st61_slope_pred <- c()
for(i in 1:362){
  knn_1tR_k1_st61_slope_pred <- c(knn_1tR_k1_st61_slope_pred, assigning_closest_cal_comp_slope(knn_1tR_k1_st61_pH_2.7_E[i,]$ret_time, standards))
}



knn_1tR_k1_st61_pH_2.7 <- data.frame(knn_1tR_k1_st61_pH_2.7_E,knn_1tR_k1_st61_slope_pred)

knn_1tR_k1_st61_pH_2.7 <- knn_1tR_k1_st61_pH_2.7 %>%
  mutate(knn_1tR_k1_st61_c_pred = Peak_Area_IC / knn_1tR_k1_st61_slope_pred)


dataset_2 <- knn_1tR_k1_st61_pH_2.7 %>%
  mutate(error_factor = case_when(knn_1tR_k1_st61_c_pred < conc_M ~ conc_M/knn_1tR_k1_st61_c_pred,
                                  knn_1tR_k1_st61_c_pred > conc_M ~ knn_1tR_k1_st61_c_pred/conc_M))%>%
  mutate(less_than_ten = case_when(error_factor < 10 ~ TRUE,
                                   error_factor > 10 ~ FALSE))
mean(dataset_2$error_factor)  
median(dataset_2$error_factor)
max(dataset_2$error_factor)
count(dataset_2, less_than_ten)

