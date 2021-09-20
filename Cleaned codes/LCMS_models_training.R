library(ggplot2)
library(cowplot)
library(randomForest)
library(caret)
library(tidyverse)
library(caTools)



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
  left_join(Unique_SMILES)%>%
  select(-MW)

setwd("~/Bc excel files/20210621")

MWs <- read_delim("Mix_M_MW.csv", delim = ";", col_names = TRUE)

dataset <- dataset %>% 
  left_join(MWs)

dataset <- dataset %>%
  filter(Compound_name != "1-Nitropyrene", Compound_name != "1-hydroxypyrene")




pH_separated_columns <- dataset %>% 
  select(Compound_name, pH, ret_time, slope, split, mode) %>%
  group_by(Compound_name, pH, slope, split, mode) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup()


pH_separated_columns_A <- pH_separated_columns %>% 
  select(-slope)
pH_separated_columns_B <- pH_separated_columns %>% 
  select(-ret_time)

pH_separated_columns_a <- pH_separated_columns_A %>%
  pivot_wider(names_from = mode, names_glue ="{mode}_{.value}", values_from = ret_time) 

pH_separated_columns_RT_pos <- pH_separated_columns_a %>% 
  select(-negative_ret_time)
pH_separated_columns_RT_pos <- pH_separated_columns_RT_pos %>%
  pivot_wider(names_from = pH, names_glue = "{pH}_{.value}", values_from = positive_ret_time) 

pH_separated_columns_RT_neg <- pH_separated_columns_a %>% 
  select(-positive_ret_time)
pH_separated_columns_RT_neg <- pH_separated_columns_RT_neg %>%
  pivot_wider(names_from = pH, names_glue = "{pH}_{.value}", values_from = negative_ret_time) 

pH_separated_columns_c <- pH_separated_columns_B %>% 
  mutate(pH = str_c("slope_", pH, sep = ""))
pH_separated_columns_b <- pH_separated_columns_c %>% 
  pivot_wider(names_from = pH, values_from = slope)

pH_separated_columns_slope_27 <- pH_separated_columns_b %>% 
  select(-slope_10, -slope_8)
pH_separated_columns_slope_27 <- pH_separated_columns_slope_27 %>%
  pivot_wider(names_from = mode, names_glue = "{mode}_{.value}", values_from = "slope_2.7") 

pH_separated_columns_slope_8 <- pH_separated_columns_b %>% 
  select(-slope_10, -slope_2.7)
pH_separated_columns_slope_8 <- pH_separated_columns_slope_8 %>%
  pivot_wider(names_from = mode, names_glue = "{mode}_{.value}", values_from = "slope_8")

pH_separated_columns_slope_10 <- pH_separated_columns_b %>% 
  select(-slope_2.7, -slope_8)
pH_separated_columns_slope_10 <- pH_separated_columns_slope_10 %>%
  pivot_wider(names_from = mode, names_glue = "{mode}_{.value}", values_from = "slope_10") 

all_pH_separated_columns <- pH_separated_columns_slope_27 %>% 
  left_join(pH_separated_columns_slope_8)%>%
  left_join(pH_separated_columns_slope_10)%>%
  left_join(pH_separated_columns_RT_pos)%>%
  left_join(pH_separated_columns_RT_neg)



molecular_weights <- dataset %>% 
  select(Compound_name, MW)%>%
  unique()


molecular_weights <- molecular_weights %>% 
  mutate(integer = as.integer(MW))

is.odd <- function(MW) MW %% 2 !=0 
molecular_weights <- molecular_weights %>% 
  mutate(odd_mass = case_when(is.odd(molecular_weights$integer)== TRUE ~ TRUE,
                              is.odd(molecular_weights$integer)== FALSE ~ FALSE))

all_pH_separated_columns <- all_pH_separated_columns %>% 
  left_join(molecular_weights)%>%
  select(-integer)

all_pH_separated_columns <- all_pH_separated_columns %>% 
  rename(x2.7_positive_ret_time = '2.7_positive_ret_time',
         x2.7_negative_ret_time = '2.7_negative_ret_time',
         x8_positive_ret_time = '8_positive_ret_time',
         x8_negative_ret_time = '8_negative_ret_time',
         x10_positive_ret_time = '10_positive_ret_time',
         x10_negative_ret_time = '10_negative_ret_time')

all_pH_separated_columns <- all_pH_separated_columns %>% 
  mutate(tR_pH_27 = case_when(is.na(x2.7_positive_ret_time) ~ x2.7_negative_ret_time,
                              is.na(x2.7_negative_ret_time) ~ x2.7_positive_ret_time,
                              TRUE ~(x2.7_positive_ret_time + x2.7_negative_ret_time)/2))%>%  
  
  mutate(tR_pH_8 = case_when(is.na(all_pH_separated_columns$x8_positive_ret_time) ~ all_pH_separated_columns$x8_negative_ret_time,
                             is.na(all_pH_separated_columns$x8_negative_ret_time) ~ all_pH_separated_columns$x8_positive_ret_time,
                             TRUE ~ (all_pH_separated_columns$x8_positive_ret_time + all_pH_separated_columns$x8_negative_ret_time)/2))%>%
  
  mutate (tR_pH_10 = case_when(is.na(all_pH_separated_columns$x10_positive_ret_time) ~ all_pH_separated_columns$x10_negative_ret_time, 
                               is.na(all_pH_separated_columns$x10_negative_ret_time) ~ all_pH_separated_columns$x10_positive_ret_time,
                               TRUE ~ (all_pH_separated_columns$x10_positive_ret_time + all_pH_separated_columns$x10_negative_ret_time)/2))%>%
  select(-x2.7_negative_ret_time, -x2.7_positive_ret_time,-x8_positive_ret_time, -x8_negative_ret_time,-x10_negative_ret_time, -x10_positive_ret_time)

all_pH_separated_columns <- all_pH_separated_columns %>% 
  mutate(pos_neg_ratio_2.7 = case_when(is.na(positive_slope_2.7)  ~ -999,
                                       is.na(negative_slope_2.7)  ~ 999,
                                       TRUE ~ log10(positive_slope_2.7/negative_slope_2.7)),
         pos_neg_ratio_8 = case_when(is.na(positive_slope_8) ~ -999,
                                     is.na(negative_slope_8) ~ 999,
                                     TRUE ~ log10(positive_slope_8/negative_slope_8)),
         pos_neg_ratio_10 = case_when(is.na(positive_slope_10) ~ -999,
                                      is.na(negative_slope_10) ~ 999,
                                      TRUE ~ log10(positive_slope_10/negative_slope_10)))

all_pH_separated_columns <- all_pH_separated_columns %>% 
  mutate(Is_detected_neg = case_when(!is.na(negative_slope_2.7) | !is.na(negative_slope_8) | !is.na(negative_slope_10) ~ TRUE,
                                     TRUE ~ FALSE),
         diff_tR_8_2.7 = tR_pH_8 - tR_pH_27,
         diff_tR_8_10 = tR_pH_8 - tR_pH_10)



model_development_set_27_pos <- all_pH_separated_columns %>%
  select(-positive_slope_10, -negative_slope_10, -positive_slope_8, -negative_slope_8, -negative_slope_2.7, )%>%
  na.omit()%>%
  mutate(positive_log_slope_2.7 = log10(positive_slope_2.7))%>%
  left_join(sodium_adducts)


model_validation_set_27_pos <- all_pH_separated_columns %>%
  select(-positive_slope_10, -negative_slope_10, -positive_slope_8, -negative_slope_8, -negative_slope_2.7, )%>%
  na.omit()%>%
  mutate(positive_log_slope_2.7 = log10(positive_slope_2.7))%>%
  left_join(sodium_adducts)

model_development_set_27_neg <- all_pH_separated_columns %>%
  select(-positive_slope_10, -negative_slope_10, -positive_slope_8, -negative_slope_8, -positive_slope_2.7, )%>%
  na.omit()%>%
  mutate(negative_log_slope_2.7 = log10(negative_slope_2.7))%>%
  left_join(sodium_adducts)

model_validation_set_27_neg <- all_pH_separated_columns %>%
  select(-positive_slope_10, -negative_slope_10, -positive_slope_8, -negative_slope_8, -positive_slope_2.7, )%>%
  na.omit()%>%
  mutate(negative_log_slope_2.7 = log10(negative_slope_2.7))%>%
  left_join(sodium_adducts)

model_development_set_8_pos <- all_pH_separated_columns %>%
  select(-positive_slope_10, -negative_slope_10, -positive_slope_2.7, -negative_slope_8, -negative_slope_2.7, )%>%
  na.omit()%>%
  mutate(positive_log_slope_8 = log10(positive_slope_8))%>%
  left_join(sodium_adducts)

model_validation_set_8_pos <- all_pH_separated_columns %>%
  select(-positive_slope_10, -negative_slope_10, -positive_slope_2.7, -negative_slope_8, -negative_slope_2.7, )%>%
  na.omit()%>%
  mutate(positive_log_slope_8 = log10(positive_slope_8))%>%
  left_join(sodium_adducts)

model_development_set_8_neg <- all_pH_separated_columns %>%
  select(-positive_slope_10, -negative_slope_10, -negative_slope_2.7, -positive_slope_8, -positive_slope_2.7, )%>%
  na.omit()%>%
  mutate(negative_log_slope_8 = log10(negative_slope_8))%>%
  left_join(sodium_adducts)

model_validation_set_8_neg <- all_pH_separated_columns %>%
  select(-positive_slope_10, -negative_slope_10, -negative_slope_2.7, -positive_slope_8, -positive_slope_2.7, )%>%
  na.omit()%>%
  mutate(negative_log_slope_8 = log10(negative_slope_8))%>%
  left_join(sodium_adducts)

model_development_set_10_pos <- all_pH_separated_columns %>%
  select(-positive_slope_8, -negative_slope_10, -positive_slope_2.7, -negative_slope_8, -negative_slope_2.7, )%>%
  na.omit()%>%
  mutate(positive_log_slope_10 = log10(positive_slope_10))%>%
  left_join(sodium_adducts)

model_validation_set_10_pos <- all_pH_separated_columns %>%
  select(-positive_slope_8, -negative_slope_10, -positive_slope_2.7, -negative_slope_8, -negative_slope_2.7, )%>%
  na.omit()%>%
  mutate(positive_log_slope_10 = log10(positive_slope_10))%>%
  left_join(sodium_adducts)

model_development_set_10_neg <- all_pH_separated_columns %>%
  select(-positive_slope_10, -negative_slope_8, -negative_slope_2.7, -positive_slope_8, -positive_slope_2.7, )%>%
  filter(Compound_name != "1-hydroxypyrene")%>%
  na.omit()%>%
  mutate(negative_log_slope_10 = log10(negative_slope_10))%>%
  left_join(sodium_adducts)

model_validation_set_10_neg <- all_pH_separated_columns %>%
  select(-positive_slope_10, -negative_slope_8, -negative_slope_2.7, -positive_slope_8, -positive_slope_2.7, )%>%
  filter(Compound_name != "1-hydroxypyrene")%>%
  na.omit()%>%
  mutate(negative_log_slope_10 = log10(negative_slope_10))%>%
  left_join(sodium_adducts)

fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5)



Random_forest_regressor_pH27_pos <- train(`positive_log_slope_2.7` ~ `MW` + `odd_mass` + `tR_pH_27` + `tR_pH_8` + `tR_pH_10` + 
                                            `Is_detected_neg` + `diff_tR_8_2.7` + `diff_tR_8_10` + `pos_neg_ratio_2.7` + 
                                            `pos_neg_ratio_8` + `pos_neg_ratio_10` + `Sodium_adduct_2.7`,
                                          data = model_development_set_27_pos,
                                          method = "RRF",
                                          trControl = fitControl,
                                          verbose = FALSE)
Random_forest_regressor_pH27_neg <- train(`negative_log_slope_2.7` ~ `MW` + `odd_mass` + `tR_pH_27` + `tR_pH_8` + `tR_pH_10` + 
                                            `Is_detected_neg` + `diff_tR_8_2.7` + `diff_tR_8_10`+ `pos_neg_ratio_2.7` + 
                                            `pos_neg_ratio_8` + `pos_neg_ratio_10`+ `Sodium_adduct_2.7`,
                                          data = model_development_set_27_neg,
                                          method = "RRF",
                                          trControl = fitControl,
                                          verbose = FALSE)
Random_forest_regressor_pH8_pos <- train(`positive_log_slope_8` ~ `MW` + `odd_mass` + `tR_pH_27` + `tR_pH_8` + `tR_pH_10` + 
                                           `Is_detected_neg` + `diff_tR_8_2.7` + `diff_tR_8_10` + `pos_neg_ratio_2.7` + 
                                           `pos_neg_ratio_8` + `pos_neg_ratio_10`+ `Sodium_adduct_2.7`,
                                         data = model_development_set_8_pos,
                                         method = "RRF",
                                         trControl = fitControl,
                                         verbose = FALSE)
Random_forest_regressor_pH8_neg <- train(`negative_log_slope_8` ~ `MW` + `odd_mass` + `tR_pH_27` + `tR_pH_8` + `tR_pH_10` + 
                                           `Is_detected_neg` + `diff_tR_8_2.7` + `diff_tR_8_10` + `pos_neg_ratio_2.7` + 
                                           `pos_neg_ratio_8` + `pos_neg_ratio_10`+ `Sodium_adduct_2.7`,
                                         data = model_development_set_8_neg,
                                         method = "RRF",
                                         trControl = fitControl,
                                         verbose = FALSE)
Random_forest_regressor_pH10_pos <- train(`positive_log_slope_10` ~ `MW` + `odd_mass` + `tR_pH_27` + `tR_pH_8` + `tR_pH_10` + 
                                            `Is_detected_neg` + `diff_tR_8_2.7` + `diff_tR_8_10` + `pos_neg_ratio_2.7` + 
                                            `pos_neg_ratio_8` + `pos_neg_ratio_10`+ `Sodium_adduct_2.7`,
                                          data = model_development_set_10_pos,
                                          method = "RRF",
                                          trControl = fitControl,
                                          verbose = FALSE)
Random_forest_regressor_pH10_neg <- train(`negative_log_slope_10` ~ `MW` + `odd_mass` + `tR_pH_27` + `tR_pH_8` + `tR_pH_10` + 
                                            `Is_detected_neg` + `diff_tR_8_2.7` + `diff_tR_8_10` + `pos_neg_ratio_2.7` + 
                                            `pos_neg_ratio_8` + `pos_neg_ratio_10`+ `Sodium_adduct_2.7`,
                                          data = model_development_set_10_neg,
                                          method = "RRF",
                                          trControl = fitControl,
                                          verbose = FALSE)


saveRDS(Random_forest_regressor_pH27_pos, file = "Random_forest_regressor_pH27_pos.rds")
saveRDS(Random_forest_regressor_pH27_neg, file = "Random_forest_regressor_pH27_neg.rds")
saveRDS(Random_forest_regressor_pH8_pos, file = "Random_forest_regressor_pH8_pos.rds")
saveRDS(Random_forest_regressor_pH8_neg, file = "Random_forest_regressor_pH8_neg.rds")
saveRDS(Random_forest_regressor_pH10_pos, file = "Random_forest_regressor_pH10_pos.rds")
saveRDS(Random_forest_regressor_pH10_neg, file = "Random_forest_regressor_pH10_neg.rds")
