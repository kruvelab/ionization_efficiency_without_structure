library(ggplot2)
library(cowplot)
library(randomForest)
library(caret)
library(tidyverse)
library(caTools)


setwd("~/Bc excel files/20210621")
dataset <- read_delim("Mix_M_dataset 20210621.csv", delim = ";", col_names = TRUE) 
dataset_E <- read_delim("New_Mix_E_dataset 20210802.csv", delim = ";", col_names = TRUE) 

dataset_2 <- dataset %>%
  group_by(pH, mode)%>%
  mutate(MS = mean(log10(slope), na.rm = TRUE))%>%
  ungroup()
dataset_3 <- dataset_E %>%
  mutate(MS = case_when(pH == 2.7 & mode == "positive" ~ 15.87717,
                        pH == 2.7 & mode == "negative" ~ 14.98798,
                        pH == 8 & mode == "positive" ~ 15.80952,
                        pH == 8 & mode == "negative" ~ 14.8875,
                        pH == 10 & mode == "positive" ~ 15.56404,
                        pH == 10 & mode == "negative" ~ 14.98628))

dataset_3 <- dataset_3 %>%
  mutate(c_pred_mean = Peak_Area_IC/10^MS)

dataset_3 <- dataset_3 %>%
  mutate(error_factor = case_when(c_pred_mean < conc_M ~ conc_M/c_pred_mean,
                                  c_pred_mean > conc_M ~ c_pred_mean/conc_M))%>%
  mutate(less_than_ten = case_when(error_factor < 10 ~ TRUE,
                                   error_factor > 10 ~ FALSE))


mean(dataset_3$error_factor)  
median(dataset_3$error_factor)
max(dataset_3$error_factor)
count(dataset_3, less_than_ten)
