library(ggplot2)
library(cowplot)
library(randomForest)
library(caret)
library(tidyverse)
library(caTools)

#setwd("~/Bc excel files/20210621")


Random_forest_regressor_pH27_pos <- readRDS(file = "Random_forest_regressor_pH27_pos.rds")
Random_forest_regressor_pH27_neg <- readRDS(file = "Random_forest_regressor_pH27_neg.rds")
Random_forest_regressor_pH8_pos <- readRDS(file = "Random_forest_regressor_pH8_pos.rds")
Random_forest_regressor_pH8_neg <- readRDS(file = "Random_forest_regressor_pH8_neg.rds")
Random_forest_regressor_pH10_pos <- readRDS(file = "Random_forest_regressor_pH10_pos.rds")
Random_forest_regressor_pH10_neg <- readRDS( file = "Random_forest_regressor_pH10_neg.rds")


#setwd("~/Bc excel files")
IC <- read_delim("IC for SB.csv", delim = ";", col_names = TRUE)
dataset_E <- read_delim("SB dataset.csv", delim = ";", col_names = TRUE)%>%
  select(-Compound)
dataset_E <- dataset_E%>%
  left_join(IC)
dataset_E <- dataset_E %>%
  group_by(Compound_name, pH, mode, type, Dilution)%>%
  mutate(ret_time = mean(RT))%>%
  mutate(Peak_Area_IC = sum(Area)*IC)%>%
    ungroup()%>%
  filter(type == "Tap")


sodium_adducts_E <- dataset_E %>%
  select(Compound_name, `Sodium adducts`)%>%
  mutate(Sodium_adduct_2.7 = `Sodium adducts`)%>%
  select(-`Sodium adducts`)

dataset_E <- dataset_E[!(dataset_E[,8]== "Clotrimazole" & dataset_E[,6]== "negative"),]


pH_separated_columns_E <- dataset_E %>%
  filter(Dilution == 1000)%>%
  select(Compound_name, pH, ret_time, Peak_Area_IC, mode, type) %>%
  group_by(Compound_name, pH, Peak_Area_IC, mode) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup() 


pH_separated_columns_E_A <- pH_separated_columns_E %>% 
  select(-Peak_Area_IC)
pH_separated_columns_E_B <- pH_separated_columns_E %>% 
  select(-ret_time)


pH_separated_columns_E_a <- pH_separated_columns_E_A %>%
  group_by(Compound_name,pH, mode)%>%
  mutate(ret_time == mean(ret_time))%>%
  ungroup()
pH_separated_columns_E_a <- pH_separated_columns_E_a %>%
  pivot_wider(names_from = mode, names_glue ="{mode}_{.value}", values_from = ret_time) 

pH_separated_columns_E_RT_pos <- pH_separated_columns_E_a %>% 
  select(-negative_ret_time)
pH_separated_columns_E_RT_pos <- pH_separated_columns_E_RT_pos %>%
  pivot_wider(names_from = pH, names_glue = "{pH}_{.value}", values_from = positive_ret_time) 

pH_separated_columns_E_RT_neg <- pH_separated_columns_E_a %>% 
  select(-positive_ret_time)
pH_separated_columns_E_RT_neg <- pH_separated_columns_E_RT_neg %>%
  pivot_wider(names_from = pH, names_glue = "{pH}_{.value}", values_from = negative_ret_time) 

pH_separated_columns_E_c <- pH_separated_columns_E_B %>% 
  mutate(pH = str_c("Peak_Area_IC_", pH, sep = ""))
pH_separated_columns_E_b <- pH_separated_columns_E_c %>% 
  pivot_wider(names_from = pH, values_from = Peak_Area_IC)

pH_separated_columns_E_slope_27 <- pH_separated_columns_E_b %>% 
  select(-Peak_Area_IC_10, -Peak_Area_IC_8)
pH_separated_columns_E_slope_27 <- pH_separated_columns_E_slope_27 %>% 
  pivot_wider(names_from = mode, names_glue = "{mode}_{.value}", values_from = "Peak_Area_IC_2.7") 

pH_separated_columns_E_slope_8 <- pH_separated_columns_E_b %>%
  select(-Peak_Area_IC_10, -Peak_Area_IC_2.7)
pH_separated_columns_E_slope_8 <- pH_separated_columns_E_slope_8 %>%
  pivot_wider(names_from = mode, names_glue = "{mode}_{.value}", values_from = "Peak_Area_IC_8") 

pH_separated_columns_E_slope_10 <- pH_separated_columns_E_b %>%
  select(-Peak_Area_IC_2.7, -Peak_Area_IC_8)
pH_separated_columns_E_slope_10 <- pH_separated_columns_E_slope_10 %>%
  pivot_wider(names_from = mode, names_glue = "{mode}_{.value}", values_from = "Peak_Area_IC_10") 

all_pH_separated_columns_E <- pH_separated_columns_E_slope_27 %>% 
  left_join(pH_separated_columns_E_slope_8)%>%
  left_join(pH_separated_columns_E_slope_10)%>%
  left_join(pH_separated_columns_E_RT_pos)%>%
  left_join(pH_separated_columns_E_RT_neg)

molecular_weights_E <- dataset_E %>% 
  select(Compound_name, MW)%>%
  unique()


molecular_weights_E <- molecular_weights_E %>% 
  mutate(integer = as.integer(MW))

is.odd <- function(MW) MW %% 2 !=0 
molecular_weights_E <- molecular_weights_E %>%
  mutate(odd_mass = case_when(is.odd(molecular_weights_E$integer)== TRUE ~ TRUE,
                              is.odd(molecular_weights_E$integer)== FALSE ~ FALSE)) 

all_pH_separated_columns_E <- all_pH_separated_columns_E %>% 
  left_join(molecular_weights_E)%>%
  select(-integer)

all_pH_separated_columns_E <- all_pH_separated_columns_E %>% 
  rename(x2.7_positive_ret_time = '2.7_positive_ret_time',
         x2.7_negative_ret_time = '2.7_negative_ret_time',
         x8_positive_ret_time = '8_positive_ret_time',
         x8_negative_ret_time = '8_negative_ret_time',
         x10_positive_ret_time = '10_positive_ret_time',
         x10_negative_ret_time = '10_negative_ret_time')

all_pH_separated_columns_E <- all_pH_separated_columns_E %>%  
  mutate(tR_pH_27 = case_when(is.na(x2.7_positive_ret_time) ~ x2.7_negative_ret_time,
                              is.na(x2.7_negative_ret_time) ~ x2.7_positive_ret_time,
                              TRUE ~(x2.7_positive_ret_time + x2.7_negative_ret_time)/2))%>%
  
  mutate(tR_pH_8 = case_when(is.na(all_pH_separated_columns_E$x8_positive_ret_time) ~ all_pH_separated_columns_E$x8_negative_ret_time,
                             is.na(all_pH_separated_columns_E$x8_negative_ret_time) ~ all_pH_separated_columns_E$x8_positive_ret_time,
                             TRUE ~ (all_pH_separated_columns_E$x8_positive_ret_time + all_pH_separated_columns_E$x8_negative_ret_time)/2))%>%
  
  mutate (tR_pH_10 = case_when(is.na(all_pH_separated_columns_E$x10_positive_ret_time) ~ all_pH_separated_columns_E$x10_negative_ret_time, 
                               is.na(all_pH_separated_columns_E$x10_negative_ret_time) ~ all_pH_separated_columns_E$x10_positive_ret_time,
                               TRUE ~ (all_pH_separated_columns_E$x10_positive_ret_time + all_pH_separated_columns_E$x10_negative_ret_time)/2))%>%
  select(-x2.7_negative_ret_time, -x2.7_positive_ret_time,-x8_positive_ret_time, -x8_negative_ret_time,-x10_negative_ret_time, -x10_positive_ret_time)


all_pH_separated_columns_E <- all_pH_separated_columns_E %>% 
  mutate(pos_neg_ratio_2.7 = case_when(is.na(positive_Peak_Area_IC_2.7)  ~ -999,
                                       is.na(negative_Peak_Area_IC_2.7)  ~ 999,
                                       TRUE ~ log10(positive_Peak_Area_IC_2.7/negative_Peak_Area_IC_2.7)),
         pos_neg_ratio_8 = case_when(is.na(positive_Peak_Area_IC_8) ~ -999,
                                     is.na(negative_Peak_Area_IC_8) ~ 999,
                                     TRUE ~ log10(positive_Peak_Area_IC_8/negative_Peak_Area_IC_8)),
         pos_neg_ratio_10 = case_when(is.na(positive_Peak_Area_IC_10) ~ -999,
                                      is.na(negative_Peak_Area_IC_10) ~ 999,
                                      TRUE ~ log10(positive_Peak_Area_IC_10/negative_Peak_Area_IC_10)))

all_pH_separated_columns_E <- all_pH_separated_columns_E %>% 
  mutate(Is_detected_neg = case_when(!is.na(negative_Peak_Area_IC_2.7) | !is.na(negative_Peak_Area_IC_8) | !is.na(negative_Peak_Area_IC_10) ~ TRUE,
                                     TRUE ~ FALSE),
         diff_tR_8_2.7 = tR_pH_8 - tR_pH_27,
         diff_tR_8_10 = tR_pH_8 - tR_pH_10)

all_pH_separated_columns_E <- all_pH_separated_columns_E %>% 
  mutate(log_positive_Peak_Area_IC_2.7 = log10(positive_Peak_Area_IC_2.7),
         log_negative_Peak_Area_IC_2.7 = log10(negative_Peak_Area_IC_2.7),
         log_positive_Peak_Area_IC_8 = log10(positive_Peak_Area_IC_8),
         log_negative_Peak_Area_IC_8 = log10(negative_Peak_Area_IC_8),
         log_positive_Peak_Area_IC_10 = log10(positive_Peak_Area_IC_10),
         log_negative_Peak_Area_IC_10 = log10(negative_Peak_Area_IC_10))


model_validation_set_E_27_pos <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_10, -negative_Peak_Area_IC_10, -positive_Peak_Area_IC_8, -negative_Peak_Area_IC_8, -negative_Peak_Area_IC_2.7, 
         -positive_Peak_Area_IC_10, -log_negative_Peak_Area_IC_10, -log_positive_Peak_Area_IC_8, -log_negative_Peak_Area_IC_8, -log_negative_Peak_Area_IC_2.7)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()


model_validation_set_E_27_neg <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_10, -negative_Peak_Area_IC_10, -positive_Peak_Area_IC_8, -negative_Peak_Area_IC_8, -positive_Peak_Area_IC_2.7, 
         -log_positive_Peak_Area_IC_10, -log_negative_Peak_Area_IC_10, -log_positive_Peak_Area_IC_8, -log_negative_Peak_Area_IC_8, -log_positive_Peak_Area_IC_2.7,)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()



model_validation_set_E_8_pos <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_10, -negative_Peak_Area_IC_10, -positive_Peak_Area_IC_2.7, -negative_Peak_Area_IC_8, -negative_Peak_Area_IC_2.7, 
         -log_positive_Peak_Area_IC_10, -log_negative_Peak_Area_IC_10, -log_positive_Peak_Area_IC_2.7, -log_negative_Peak_Area_IC_8, -log_negative_Peak_Area_IC_2.7,)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()



model_validation_set_E_8_neg <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_10, -negative_Peak_Area_IC_10, -negative_Peak_Area_IC_2.7, -positive_Peak_Area_IC_8, -positive_Peak_Area_IC_2.7,
         -log_positive_Peak_Area_IC_10, -log_negative_Peak_Area_IC_10, -log_negative_Peak_Area_IC_2.7, -log_positive_Peak_Area_IC_8, -log_positive_Peak_Area_IC_2.7,)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()


model_validation_set_E_10_pos <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_8, -negative_Peak_Area_IC_10, -positive_Peak_Area_IC_2.7, -negative_Peak_Area_IC_8, -negative_Peak_Area_IC_2.7,
         -log_positive_Peak_Area_IC_8, -log_negative_Peak_Area_IC_10, -log_positive_Peak_Area_IC_2.7, -log_negative_Peak_Area_IC_8, -log_negative_Peak_Area_IC_2.7,)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()


model_validation_set_E_10_neg <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_10, -negative_Peak_Area_IC_8, -negative_Peak_Area_IC_2.7, -positive_Peak_Area_IC_8, -positive_Peak_Area_IC_2.7,
         -log_positive_Peak_Area_IC_10, -log_negative_Peak_Area_IC_8, -log_negative_Peak_Area_IC_2.7, -log_positive_Peak_Area_IC_8, -log_positive_Peak_Area_IC_2.7,)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()



model_validation_set_E_27_pos <- model_validation_set_E_27_pos %>% 
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH27_pos, newdata = model_validation_set_E_27_pos))
model_validation_set_E_27_neg <- model_validation_set_E_27_neg %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH27_neg, newdata = model_validation_set_E_27_neg))
model_validation_set_E_8_pos <- model_validation_set_E_8_pos %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH8_pos, newdata = model_validation_set_E_8_pos))
model_validation_set_E_8_neg <- model_validation_set_E_8_neg %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH8_neg, newdata = model_validation_set_E_8_neg))
model_validation_set_E_10_pos <- model_validation_set_E_10_pos %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH10_pos, newdata = model_validation_set_E_10_pos))
model_validation_set_E_10_neg <- model_validation_set_E_10_neg %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH10_neg, newdata = model_validation_set_E_10_neg))


predicted_slopes_27_pos_val <- model_validation_set_E_27_pos %>% 
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 2.7,
         mode = "positive")
predicted_slopes_27_neg_val <- model_validation_set_E_27_neg %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 2.7,
         mode = "negative")
predicted_slopes_8_pos_val <- model_validation_set_E_8_pos %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 8,
         mode = "positive")
predicted_slopes_8_neg_val <- model_validation_set_E_8_neg %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 8,
         mode = "negative")
predicted_slopes_10_pos_val <- model_validation_set_E_10_pos %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 10,
         mode = "positive")
predicted_slopes_10_neg_val <- model_validation_set_E_10_neg %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 10,
         mode = "negative")

validation_dataset_27_pos <- dataset_E %>% 
  filter(pH ==2.7, Dilution == 1000)%>%
  filter(mode == "positive")
validation_dataset_27_pos <- validation_dataset_27_pos %>%
  left_join(predicted_slopes_27_pos_val)
validation_dataset_27_neg <- dataset_E %>%
  filter(pH ==2.7, Dilution == 1000)%>%
  filter(mode == "negative")
validation_dataset_27_neg <- validation_dataset_27_neg %>%
  left_join(predicted_slopes_27_neg_val)
validation_dataset_8_pos <- dataset_E %>%
  filter(pH ==8, Dilution == 1000)%>%
  filter(mode == "positive")
validation_dataset_8_pos <- validation_dataset_8_pos %>%
  left_join(predicted_slopes_8_pos_val)
validation_dataset_8_neg <- dataset_E %>%
  filter(pH ==8, Dilution == 1000)%>%
  filter(mode == "negative")
validation_dataset_8_neg <- validation_dataset_8_neg %>%
  left_join(predicted_slopes_8_neg_val)
validation_dataset_10_pos <- dataset_E %>%
  filter(pH ==10, Dilution == 1000)%>%
  filter(mode == "positive")
validation_dataset_10_pos <- validation_dataset_10_pos %>%
  left_join(predicted_slopes_10_pos_val)
validation_dataset_10_neg <- dataset_E %>%
  filter(pH ==10, Dilution == 1000)%>%
  filter(mode == "negative")
validation_dataset_10_neg <- validation_dataset_10_neg %>%
  left_join(predicted_slopes_10_neg_val)

concentrations <- read_delim("conc_validation_emma.csv", delim = ";", col_names = TRUE)%>%
  mutate(conc_M = conc_in_sample)%>%
  mutate(Compound_name = Compound,
         type = case_when(type == "tap" ~ "Tap",
                          TRUE ~ type))%>%
  select(Compound_name, Dilution, type, Mix, conc_M)%>%
  unique()


validation_dataset_27_pos <- validation_dataset_27_pos %>%
  mutate(predicted_slope = 10^(predicted_log_slope))%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_27_pos <- validation_dataset_27_pos %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  TRUE ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)

validation_dataset_27_neg <- validation_dataset_27_neg %>%
  mutate(predicted_slope = 10^predicted_log_slope)%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_27_neg <- validation_dataset_27_neg %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  c_pred > conc_M ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)

validation_dataset_8_pos <- validation_dataset_8_pos %>%
  mutate(predicted_slope = 10^predicted_log_slope)%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_8_pos <- validation_dataset_8_pos %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  c_pred > conc_M ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW,  predicted_log_slope, c_pred, conc_M)

validation_dataset_8_neg <- validation_dataset_8_neg %>%
  mutate(predicted_slope = 10^predicted_log_slope)%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_8_neg <- validation_dataset_8_neg %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  c_pred > conc_M ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)

validation_dataset_10_pos <- validation_dataset_10_pos %>%
  mutate(predicted_slope = 10^predicted_log_slope)%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_10_pos <- validation_dataset_10_pos %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  c_pred > conc_M ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)

validation_dataset_10_neg <- validation_dataset_10_neg %>%
  mutate(predicted_slope = 10^predicted_log_slope)%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_10_neg <- validation_dataset_10_neg %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  c_pred > conc_M ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)


all_errors_1000 <- validation_dataset_27_pos %>% 
  rbind(validation_dataset_27_neg) %>%
  rbind(validation_dataset_8_pos)%>%
  rbind(validation_dataset_8_neg)%>%
  rbind(validation_dataset_10_pos)%>%
  rbind(validation_dataset_10_neg)%>%
  mutate(integer = as.integer(MW))

all_errors_1000 <- all_errors_1000 %>% 
  mutate(less_than_ten = case_when(error_factor < 10 ~ TRUE,
                                   error_factor > 10 ~ FALSE))



pH_separated_columns_E <- dataset_E %>%
  filter(Dilution == 500)%>%
  select(Compound_name, pH, ret_time, Peak_Area_IC, mode) %>%
  group_by(Compound_name, pH, Peak_Area_IC, mode) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup()


pH_separated_columns_E_A <- pH_separated_columns_E %>%
  select(-Peak_Area_IC)
pH_separated_columns_E_B <- pH_separated_columns_E %>%
  select(-ret_time)


pH_separated_columns_E_a <- pH_separated_columns_E_A %>%
  pivot_wider(names_from = mode, names_glue ="{mode}_{.value}", values_from = ret_time)

pH_separated_columns_E_RT_pos <- pH_separated_columns_E_a %>%
  select(-negative_ret_time)
pH_separated_columns_E_RT_pos <- pH_separated_columns_E_RT_pos %>%
  pivot_wider(names_from = pH, names_glue = "{pH}_{.value}", values_from = positive_ret_time)

pH_separated_columns_E_RT_neg <- pH_separated_columns_E_a %>%
  select(-positive_ret_time)
pH_separated_columns_E_RT_neg <- pH_separated_columns_E_RT_neg %>%
  pivot_wider(names_from = pH, names_glue = "{pH}_{.value}", values_from = negative_ret_time)

pH_separated_columns_E_c <- pH_separated_columns_E_B %>%
  mutate(pH = str_c("Peak_Area_IC_", pH, sep = ""))
pH_separated_columns_E_b <- pH_separated_columns_E_c %>%
  pivot_wider(names_from = pH, values_from = Peak_Area_IC)

pH_separated_columns_E_slope_27 <- pH_separated_columns_E_b %>%
  select(-Peak_Area_IC_10, -Peak_Area_IC_8)
pH_separated_columns_E_slope_27 <- pH_separated_columns_E_slope_27 %>%
  pivot_wider(names_from = mode, names_glue = "{mode}_{.value}", values_from = "Peak_Area_IC_2.7")

pH_separated_columns_E_slope_8 <- pH_separated_columns_E_b %>%
  select(-Peak_Area_IC_10, -Peak_Area_IC_2.7)
pH_separated_columns_E_slope_8 <- pH_separated_columns_E_slope_8 %>%
  pivot_wider(names_from = mode, names_glue = "{mode}_{.value}", values_from = "Peak_Area_IC_8")

pH_separated_columns_E_slope_10 <- pH_separated_columns_E_b %>%
  select(-Peak_Area_IC_2.7, -Peak_Area_IC_8)
pH_separated_columns_E_slope_10 <- pH_separated_columns_E_slope_10 %>%
  pivot_wider(names_from = mode, names_glue = "{mode}_{.value}", values_from = "Peak_Area_IC_10")

all_pH_separated_columns_E <- pH_separated_columns_E_slope_27 %>%
  left_join(pH_separated_columns_E_slope_8)%>%
  left_join(pH_separated_columns_E_slope_10)%>%
  left_join(pH_separated_columns_E_RT_pos)%>%
  left_join(pH_separated_columns_E_RT_neg)

molecular_weights_E <- dataset_E %>%
  select(Compound_name, MW)%>%
  unique()


molecular_weights_E <- molecular_weights_E %>%
  mutate(integer = as.integer(MW))

is.odd <- function(MW) MW %% 2 !=0
molecular_weights_E <- molecular_weights_E %>%
  mutate(odd_mass = case_when(is.odd(molecular_weights_E$integer)== TRUE ~ TRUE,
                              is.odd(molecular_weights_E$integer)== FALSE ~ FALSE))

all_pH_separated_columns_E <- all_pH_separated_columns_E %>%
  left_join(molecular_weights_E)%>%
  select(-integer)

all_pH_separated_columns_E <- all_pH_separated_columns_E %>%
  rename(x2.7_positive_ret_time = '2.7_positive_ret_time',
         x2.7_negative_ret_time = '2.7_negative_ret_time',
         x8_positive_ret_time = '8_positive_ret_time',
         x8_negative_ret_time = '8_negative_ret_time',
         x10_positive_ret_time = '10_positive_ret_time',
         x10_negative_ret_time = '10_negative_ret_time')

all_pH_separated_columns_E <- all_pH_separated_columns_E %>%
  mutate(tR_pH_27 = case_when(is.na(x2.7_positive_ret_time) ~ x2.7_negative_ret_time,
                              is.na(x2.7_negative_ret_time) ~ x2.7_positive_ret_time,
                              TRUE ~(x2.7_positive_ret_time + x2.7_negative_ret_time)/2))%>%
  
  mutate(tR_pH_8 = case_when(is.na(all_pH_separated_columns_E$x8_positive_ret_time) ~ all_pH_separated_columns_E$x8_negative_ret_time,
                             is.na(all_pH_separated_columns_E$x8_negative_ret_time) ~ all_pH_separated_columns_E$x8_positive_ret_time,
                             TRUE ~ (all_pH_separated_columns_E$x8_positive_ret_time + all_pH_separated_columns_E$x8_negative_ret_time)/2))%>%
  
  mutate (tR_pH_10 = case_when(is.na(all_pH_separated_columns_E$x10_positive_ret_time) ~ all_pH_separated_columns_E$x10_negative_ret_time, 
                               is.na(all_pH_separated_columns_E$x10_negative_ret_time) ~ all_pH_separated_columns_E$x10_positive_ret_time,
                               TRUE ~ (all_pH_separated_columns_E$x10_positive_ret_time + all_pH_separated_columns_E$x10_negative_ret_time)/2))%>%
  select(-x2.7_negative_ret_time, -x2.7_positive_ret_time,-x8_positive_ret_time, -x8_negative_ret_time,-x10_negative_ret_time, -x10_positive_ret_time)


all_pH_separated_columns_E <- all_pH_separated_columns_E %>%
  mutate(pos_neg_ratio_2.7 = case_when(is.na(positive_Peak_Area_IC_2.7)  ~ -999,
                                       is.na(negative_Peak_Area_IC_2.7)  ~ 999,
                                       TRUE ~ log10(positive_Peak_Area_IC_2.7/negative_Peak_Area_IC_2.7)),
         pos_neg_ratio_8 = case_when(is.na(positive_Peak_Area_IC_8) ~ -999,
                                     is.na(negative_Peak_Area_IC_8) ~ 999,
                                     TRUE ~ log10(positive_Peak_Area_IC_8/negative_Peak_Area_IC_8)),
         pos_neg_ratio_10 = case_when(is.na(positive_Peak_Area_IC_10) ~ -999,
                                      is.na(negative_Peak_Area_IC_10) ~ 999,
                                      TRUE ~ log10(positive_Peak_Area_IC_10/negative_Peak_Area_IC_10)))

all_pH_separated_columns_E <- all_pH_separated_columns_E %>%
  mutate(Is_detected_neg = case_when(!is.na(negative_Peak_Area_IC_2.7) | !is.na(negative_Peak_Area_IC_8) | !is.na(negative_Peak_Area_IC_10) ~ TRUE,
                                     TRUE ~ FALSE),
         diff_tR_8_2.7 = tR_pH_8 - tR_pH_27,
         diff_tR_8_10 = tR_pH_8 - tR_pH_10)

all_pH_separated_columns_E <- all_pH_separated_columns_E %>%
  mutate(log_positive_Peak_Area_IC_2.7 = log10(positive_Peak_Area_IC_2.7),
         log_negative_Peak_Area_IC_2.7 = log10(negative_Peak_Area_IC_2.7),
         log_positive_Peak_Area_IC_8 = log10(positive_Peak_Area_IC_8),
         log_negative_Peak_Area_IC_8 = log10(negative_Peak_Area_IC_8),
         log_positive_Peak_Area_IC_10 = log10(positive_Peak_Area_IC_10),
         log_negative_Peak_Area_IC_10 = log10(negative_Peak_Area_IC_10))


model_validation_set_E_27_pos <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_10, -negative_Peak_Area_IC_10, -positive_Peak_Area_IC_8, -negative_Peak_Area_IC_8, -negative_Peak_Area_IC_2.7, 
         -positive_Peak_Area_IC_10, -log_negative_Peak_Area_IC_10, -log_positive_Peak_Area_IC_8, -log_negative_Peak_Area_IC_8, -log_negative_Peak_Area_IC_2.7)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()


model_validation_set_E_27_neg <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_10, -negative_Peak_Area_IC_10, -positive_Peak_Area_IC_8, -negative_Peak_Area_IC_8, -positive_Peak_Area_IC_2.7, 
         -log_positive_Peak_Area_IC_10, -log_negative_Peak_Area_IC_10, -log_positive_Peak_Area_IC_8, -log_negative_Peak_Area_IC_8, -log_positive_Peak_Area_IC_2.7,)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()



model_validation_set_E_8_pos <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_10, -negative_Peak_Area_IC_10, -positive_Peak_Area_IC_2.7, -negative_Peak_Area_IC_8, -negative_Peak_Area_IC_2.7, 
         -log_positive_Peak_Area_IC_10, -log_negative_Peak_Area_IC_10, -log_positive_Peak_Area_IC_2.7, -log_negative_Peak_Area_IC_8, -log_negative_Peak_Area_IC_2.7,)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()



model_validation_set_E_8_neg <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_10, -negative_Peak_Area_IC_10, -negative_Peak_Area_IC_2.7, -positive_Peak_Area_IC_8, -positive_Peak_Area_IC_2.7,
         -log_positive_Peak_Area_IC_10, -log_negative_Peak_Area_IC_10, -log_negative_Peak_Area_IC_2.7, -log_positive_Peak_Area_IC_8, -log_positive_Peak_Area_IC_2.7,)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()


model_validation_set_E_10_pos <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_8, -negative_Peak_Area_IC_10, -positive_Peak_Area_IC_2.7, -negative_Peak_Area_IC_8, -negative_Peak_Area_IC_2.7,
         -log_positive_Peak_Area_IC_8, -log_negative_Peak_Area_IC_10, -log_positive_Peak_Area_IC_2.7, -log_negative_Peak_Area_IC_8, -log_negative_Peak_Area_IC_2.7,)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()


model_validation_set_E_10_neg <- all_pH_separated_columns_E %>%
  select(-positive_Peak_Area_IC_10, -negative_Peak_Area_IC_8, -negative_Peak_Area_IC_2.7, -positive_Peak_Area_IC_8, -positive_Peak_Area_IC_2.7,
         -log_positive_Peak_Area_IC_10, -log_negative_Peak_Area_IC_8, -log_negative_Peak_Area_IC_2.7, -log_positive_Peak_Area_IC_8, -log_positive_Peak_Area_IC_2.7,)%>%
  na.omit()%>%
  left_join(sodium_adducts_E)%>%
  unique()



model_validation_set_E_27_pos <- model_validation_set_E_27_pos %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH27_pos, newdata = model_validation_set_E_27_pos))
model_validation_set_E_27_neg <- model_validation_set_E_27_neg %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH27_neg, newdata = model_validation_set_E_27_neg))
model_validation_set_E_8_pos <- model_validation_set_E_8_pos %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH8_pos, newdata = model_validation_set_E_8_pos))
model_validation_set_E_8_neg <- model_validation_set_E_8_neg %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH8_neg, newdata = model_validation_set_E_8_neg))
model_validation_set_E_10_pos <- model_validation_set_E_10_pos %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH10_pos, newdata = model_validation_set_E_10_pos))
model_validation_set_E_10_neg <- model_validation_set_E_10_neg %>%
  mutate(predicted_log_slope = predict(Random_forest_regressor_pH10_neg, newdata = model_validation_set_E_10_neg))


predicted_slopes_27_pos_val <- model_validation_set_E_27_pos %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 2.7,
         mode = "positive")
predicted_slopes_27_neg_val <- model_validation_set_E_27_neg %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 2.7,
         mode = "negative")
predicted_slopes_8_pos_val <- model_validation_set_E_8_pos %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 8,
         mode = "positive")
predicted_slopes_8_neg_val <- model_validation_set_E_8_neg %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 8,
         mode = "negative")
predicted_slopes_10_pos_val <- model_validation_set_E_10_pos %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 10,
         mode = "positive")
predicted_slopes_10_neg_val <- model_validation_set_E_10_neg %>%
  select(predicted_log_slope, Compound_name)%>%
  mutate(pH = 10,
         mode = "negative")

validation_dataset_27_pos <- dataset_E %>%
  filter(pH ==2.7, Dilution == 500)%>%
  filter(mode == "positive")
validation_dataset_27_pos <- validation_dataset_27_pos %>%
  left_join(predicted_slopes_27_pos_val)
validation_dataset_27_neg <- dataset_E %>%
  filter(pH ==2.7, Dilution == 500)%>%
  filter(mode == "negative")
validation_dataset_27_neg <- validation_dataset_27_neg %>%
  left_join(predicted_slopes_27_neg_val)
validation_dataset_8_pos <- dataset_E %>%
  filter(pH ==8, Dilution == 500)%>%
  filter(mode == "positive")
validation_dataset_8_pos <- validation_dataset_8_pos %>%
  left_join(predicted_slopes_8_pos_val)
validation_dataset_8_neg <- dataset_E %>%
  filter(pH ==8, Dilution == 500)%>%
  filter(mode == "negative")
validation_dataset_8_neg <- validation_dataset_8_neg %>%
  left_join(predicted_slopes_8_neg_val)
validation_dataset_10_pos <- dataset_E %>%
  filter(pH ==10, Dilution == 500)%>%
  filter(mode == "positive")
validation_dataset_10_pos <- validation_dataset_10_pos %>%
  left_join(predicted_slopes_10_pos_val)
validation_dataset_10_neg <- dataset_E %>%
  filter(pH ==10, Dilution == 500)%>%
  filter(mode == "negative")
validation_dataset_10_neg <- validation_dataset_10_neg %>%
  left_join(predicted_slopes_10_neg_val)

concentrations <- concentrations %>%
  select(-type,-Mix)

validation_dataset_27_pos <- validation_dataset_27_pos %>%
  mutate(predicted_slope = 10^(predicted_log_slope))%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_27_pos <- validation_dataset_27_pos %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  TRUE ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)

validation_dataset_27_neg <- validation_dataset_27_neg %>%
  mutate(predicted_slope = 10^predicted_log_slope)%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_27_neg <- validation_dataset_27_neg %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  c_pred > conc_M ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)

validation_dataset_8_pos <- validation_dataset_8_pos %>%
  mutate(predicted_slope = 10^predicted_log_slope)%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_8_pos <- validation_dataset_8_pos %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  c_pred > conc_M ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)

validation_dataset_8_neg <- validation_dataset_8_neg %>%
  mutate(predicted_slope = 10^predicted_log_slope)%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_8_neg <- validation_dataset_8_neg %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  c_pred > conc_M ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)

validation_dataset_10_pos <- validation_dataset_10_pos %>%
  mutate(predicted_slope = 10^predicted_log_slope)%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_10_pos <- validation_dataset_10_pos %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  c_pred > conc_M ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)

validation_dataset_10_neg <- validation_dataset_10_neg %>%
  mutate(predicted_slope = 10^predicted_log_slope)%>%
  mutate(c_pred = Peak_Area_IC/predicted_slope)%>%
  na.omit()%>%
  left_join(concentrations)
validation_dataset_10_neg <- validation_dataset_10_neg %>%
  mutate(error_factor = case_when(c_pred < conc_M ~ conc_M/c_pred,
                                  c_pred > conc_M ~ c_pred/conc_M))%>%
  select(Compound_name, pH, mode, error_factor, ret_time, Dilution, MW, predicted_log_slope, c_pred, conc_M)


all_errors_500 <- validation_dataset_27_pos %>%
  rbind(validation_dataset_27_neg) %>%
  rbind(validation_dataset_8_pos)%>%
  rbind(validation_dataset_8_neg)%>%
  rbind(validation_dataset_10_pos)%>%
  rbind(validation_dataset_10_neg)%>%
  mutate(integer = as.integer(MW))

all_errors_500 <- all_errors_500 %>%
  mutate(less_than_ten = case_when(error_factor < 10 ~ TRUE,
                                   error_factor > 10 ~ FALSE))




all_errors <- all_errors_1000 %>% 
  rbind(all_errors_500)

count(all_errors, less_than_ten) 

mean(all_errors$error_factor)
max(all_errors$error_factor)
median(all_errors$error_factor)

