
library(caret) #machine learinign workflow
library(caTools) #sample split is from this package
library(tidyverse) #helps us to write concise code
library(plotly)
library(extrafont)
loadfonts(device = "win")


# Importing the data with read_delim
dataset <-  read_delim('data_ret_time_200514.csv',
                       delim = ",",
                       col_names = TRUE)

#model training parameters
fitControl <- trainControl(
  method = "repeatedcv", 
  number = 5, 
  repeats = 5) 

#Creating a table with the unique SMILES
unique_SMILES <- dataset %>%
  select(SMILES)%>%
  unique()

# Slpittng the unique SMILES in to a method development and validation set
set.seed(123) 
split <- sample.split(unique_SMILES$SMILES, SplitRatio = 0.8)

#Adds the split to the unique_SMILES dataset
unique_SMILES <- unique_SMILES %>%
  mutate(split = split)

#Adds the split to the original dataset
dataset <- unique_SMILES %>%
  left_join(dataset %>% select(Compound_name, SMILES) %>% unique())

compound_splits <- dataset %>%
  select(Compound_name,split)
write_delim(compound_splits,
            "Mix_M_splits.csv",
            delim = ";")

dataset_M <-  read_delim('Dataset_Mix_M_slopes.csv',
                       delim = ";",
                       col_names = TRUE)
dataset_M <- dataset_M %>%
  left_join(compound_splits)
dataset_M <- dataset_M %>%
  filter( Compound_name != "Sudan_III",
          Compound_name != "Chloramphenicol",
          Compound_name != "Cyclamic_acid",
          Compound_name != "1,4-Naphthoquinone",
          slope > 0)

dataset_M <- dataset_M %>%
  mutate(split = case_when(
    is.na(split) ~ FALSE,
    TRUE ~ split
  ))

dataset_M <- dataset_M %>%
  mutate(Buffer = case_when(pH == "2.7" ~ "Formic acid",
                            pH == 5.0 ~ "Formiate",
                            pH == 8.0 ~ "Bicarbonate"))



dataset_M <-  dataset_M %>%
  unite("x", "pH","Buffer", remove = FALSE)
dataset_M <- dataset_M %>%
  unite("MP", "x", "Organic_modifier", remove = FALSE)

dataset_M <- dataset_M %>% 
  filter(MP == "2.7_Formic acid_Acetonitrile"|MP=="5_Formiate_Acetonitrile" | MP=="5_Formiate_Methanol" | MP=="2.3_TFA_Methanol" | MP=="8_Bicarbonate_Acetonitrile" | MP=="8_Bicarbonate_Methanol")



#Hypothesis 2------
# The linear regression model. The models use the slope from one of the  
# mobile phases and then the retention time in the same mobile phase.
# Using 61 standard compounds.


#Selects the desired columns
MPs_separated_columns <- dataset_M %>%
  select(Compound_name, MP, ret_time, slope, split) %>%
  group_by(Compound_name, MP, slope, split) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup()

#Creating the H_two_Take_all_MP dataset
MPs_separated_columns_A <- MPs_separated_columns %>%
  select(-slope)
MPs_separated_columns_B <- MPs_separated_columns %>%
  select(-ret_time)

MPs_separated_columns_a <-  spread(MPs_separated_columns_A, key = "MP", value = "ret_time")
MPs_separated_columns_c <- MPs_separated_columns_B %>%
  mutate(MP = str_c("slope_", MP, sep = ""))
MPs_separated_columns_b <- spread(MPs_separated_columns_c, key = "MP", value = "slope")

MPs_separated_columns <- MPs_separated_columns_a %>%
  left_join(MPs_separated_columns_b)%>%
  drop_na() %>%
  unique() 


# Creating the modeldevelopment and validation datasets as well as the regressor
MPs_separated_columns_model_development <- subset(MPs_separated_columns, split == TRUE) 



Linear_6tR_st_61_regressor_pH3_MeCN <- train(`slope_2.7_Formic acid_Acetonitrile` ~  `2.7_Formic acid_Acetonitrile`, 
                             data = MPs_separated_columns_model_development, 
                             method = "glmStepAIC", 
                             trControl = fitControl, 
                             verbose = FALSE)



Linear_6tR_st_61_regressor_pH5_MeCN <- train(`slope_5_Formiate_Acetonitrile` ~ `5_Formiate_Acetonitrile`, 
                             data = MPs_separated_columns_model_development, 
                             method = "glmStepAIC", 
                             trControl = fitControl, 
                             verbose = FALSE)




Linear_6tR_st_61_regressor_pH8_MeCN <- train(`slope_8_Bicarbonate_Acetonitrile` ~ `8_Bicarbonate_Acetonitrile`, 
                              data = MPs_separated_columns_model_development, 
                              method = "glmStepAIC", 
                              trControl = fitControl, 
                              verbose = FALSE)



# Predicting the slopes and adding then to te dataset
Linear_6tR_st_61_Take_all_MP_pred <- MPs_separated_columns %>% 
  mutate(Linear_6tR_st_61_slope_pred_pH3_MeCN = predict(Linear_6tR_st_61_regressor_pH3_MeCN, newdata = MPs_separated_columns),
         Linear_6tR_st_61_slope_pred_pH5_MeCN = predict(Linear_6tR_st_61_regressor_pH5_MeCN, newdata = MPs_separated_columns),
         Linear_6tR_st_61_slope_pred_pH8_MeCN = predict(Linear_6tR_st_61_regressor_pH8_MeCN, newdata = MPs_separated_columns))

#Selecte the wanted columns
Linear_6tR_st_61_Take_all_MP_pred <- Linear_6tR_st_61_Take_all_MP_pred %>%
  select(Compound_name, Linear_6tR_st_61_slope_pred_pH3_MeCN, Linear_6tR_st_61_slope_pred_pH5_MeCN, Linear_6tR_st_61_slope_pred_pH8_MeCN)

#Gaters the slope predictions
Linear_6tR_st_61_Take_all_MP_pred <- gather(Linear_6tR_st_61_Take_all_MP_pred, key = MP, value = "Linear_6tR_st_61_slope_pred", Linear_6tR_st_61_slope_pred_pH3_MeCN, Linear_6tR_st_61_slope_pred_pH5_MeCN, Linear_6tR_st_61_slope_pred_pH8_MeCN)

#Renameing
Linear_6tR_st_61_Take_all_MP_pred <- Linear_6tR_st_61_Take_all_MP_pred %>%
  mutate(MP = case_when(
    MP == "Linear_6tR_st_61_slope_pred_pH3_MeCN" ~ "2.7_Formic acid_Acetonitrile",
    MP == "Linear_6tR_st_61_slope_pred_pH5_MeCN" ~ "5_Formiate_Acetonitrile",
    MP == "Linear_6tR_st_61_slope_pred_pH8_MeCN" ~ "8_Bicarbonate_Acetonitrile"
  ))


dataset_M <- dataset_M %>%
  left_join(Linear_6tR_st_61_Take_all_MP_pred) %>%
  mutate(Linear_6tR_st_61_c_pred = Peak_Area_IC / Linear_6tR_st_61_slope_pred)

#Hypothesis 1-----
# Using the slope of the median compound or the mean slope of all compounds.
# Retention time from the same mobile phase as the slope


mean_27 <- mean(MPs_separated_columns$`slope_2.7_Formic acid_Acetonitrile`)
mean_5 <- mean(MPs_separated_columns$slope_5_Formiate_Acetonitrile)
mean_8 <- mean(MPs_separated_columns$slope_8_Bicarbonate_Acetonitrile)
median_27 <- median(MPs_separated_columns$`slope_2.7_Formic acid_Acetonitrile`)
median_5 <- median(MPs_separated_columns$slope_5_Formiate_Acetonitrile)
median_8 <- median(MPs_separated_columns$slope_8_Bicarbonate_Acetonitrile)

dataset_M <- dataset_M %>%
  mutate(mean_c = case_when(pH == 2.7 ~ Peak_Area_IC / mean_27,
                            pH == 5.0 ~ Peak_Area_IC / mean_5,
                            pH == 8.0 ~ Peak_Area_IC / mean_8))

dataset_M <- dataset_M %>%
  mutate(median_c = case_when(pH == 2.7 ~ Peak_Area_IC / median_27,
                            pH == 5.0 ~ Peak_Area_IC / median_5,
                            pH == 8.0 ~ Peak_Area_IC / median_8))


#Hypothesis 3 with three retention times----
# Multilinear regression with slopes from one of the mobile phases and retention
# times for the mobile phases with the same organic phase (3 mobile phases).
# Using 61 standard compounds.

MLR_3tR_st61_regressor_pH3_MeCN <- train(`slope_2.7_Formic acid_Acetonitrile` ~  `2.7_Formic acid_Acetonitrile` + `5_Formiate_Acetonitrile` + `8_Bicarbonate_Acetonitrile`,
                                    data = MPs_separated_columns_model_development,
                                    method = "glmStepAIC", 
                                    trControl = fitControl, 
                                    verbose = FALSE)



MLR_3tR_st61_regressor_pH5_MeCN <- train(`slope_5_Formiate_Acetonitrile` ~ `2.7_Formic acid_Acetonitrile` + `5_Formiate_Acetonitrile` + `8_Bicarbonate_Acetonitrile`, 
                                    data = MPs_separated_columns_model_development, 
                                    method = "glmStepAIC", 
                                    trControl = fitControl, 
                                    verbose = FALSE)





MLR_3tR_st61_regressor_pH8_MeCN <- train(`slope_8_Bicarbonate_Acetonitrile` ~  `2.7_Formic acid_Acetonitrile` + `5_Formiate_Acetonitrile` + `8_Bicarbonate_Acetonitrile`, 
                                     data = MPs_separated_columns_model_development, 
                                     method = "glmStepAIC", 
                                     trControl = fitControl, 
                                     verbose = FALSE)



MLR_3tR_st61_Take_all_MP <- MPs_separated_columns %>% 
  mutate(split = split, 
         MLR_3tR_st61_slope_pred_pH3_MeCN = predict(MLR_3tR_st61_regressor_pH3_MeCN, newdata = MPs_separated_columns), 
         MLR_3tR_st61_slope_pred_pH5_MeCN = predict(MLR_3tR_st61_regressor_pH5_MeCN, newdata = MPs_separated_columns), 
         MLR_3tR_st61_slope_pred_pH8_MeCN = predict(MLR_3tR_st61_regressor_pH8_MeCN, newdata = MPs_separated_columns))
        
MLR_3tR_st61_Take_all_MP <- MLR_3tR_st61_Take_all_MP %>%
  select(Compound_name, MLR_3tR_st61_slope_pred_pH3_MeCN, MLR_3tR_st61_slope_pred_pH5_MeCN, MLR_3tR_st61_slope_pred_pH8_MeCN)

MLR_3tR_st61_Take_all_MP <- gather(MLR_3tR_st61_Take_all_MP, key = MP, value = "MLR_3tR_st61_slope_pred", MLR_3tR_st61_slope_pred_pH3_MeCN, MLR_3tR_st61_slope_pred_pH5_MeCN, MLR_3tR_st61_slope_pred_pH8_MeCN)

MLR_3tR_st61_Take_all_MP <- MLR_3tR_st61_Take_all_MP %>%
  mutate(MP = case_when(
    MP == "MLR_3tR_st61_slope_pred_pH3_MeCN" ~ "2.7_Formic acid_Acetonitrile",
    MP == "MLR_3tR_st61_slope_pred_pH5_MeCN" ~ "5_Formiate_Acetonitrile",
    MP == "MLR_3tR_st61_slope_pred_pH8_MeCN" ~ "8_Bicarbonate_Acetonitrile"
  ))


dataset_M <- dataset_M %>%
  left_join(MLR_3tR_st61_Take_all_MP) %>%
  mutate(MLR_3tR_st61_c_pred = Peak_Area_IC / MLR_3tR_st61_slope_pred)




#Hypothesis three 2MPs----
# Multilinear regression with slopes from one mobile phase and retnetion times 
# for the lowest and highest pH value for that organic pahse. (2 mobile pahses).
# Using 61 standard compounds.

MLR_2tR_st61_regressor_pH3_MeCN <- train(`slope_2.7_Formic acid_Acetonitrile` ~  `2.7_Formic acid_Acetonitrile` +  `8_Bicarbonate_Acetonitrile`,
                                    data = MPs_separated_columns_model_development,
                                    method = "glmStepAIC", 
                                    trControl = fitControl, 
                                    verbose = FALSE)



MLR_2tR_st61_regressor_pH5_MeCN <- train(`slope_5_Formiate_Acetonitrile` ~ `2.7_Formic acid_Acetonitrile` +  `8_Bicarbonate_Acetonitrile`, 
                                    data = MPs_separated_columns_model_development, 
                                    method = "glmStepAIC", 
                                    trControl = fitControl, 
                                    verbose = FALSE)





MLR_2tR_st61_regressor_pH8_MeCN <- train(`slope_8_Bicarbonate_Acetonitrile` ~  `2.7_Formic acid_Acetonitrile` +  `8_Bicarbonate_Acetonitrile`, 
                                     data = MPs_separated_columns_model_development, 
                                     method = "glmStepAIC", 
                                     trControl = fitControl, 
                                     verbose = FALSE)




MLR_2tR_st61_Take_all_MP <- MPs_separated_columns %>% 
  mutate(split = split, 
         MLR_2tR_st61_slope_pred_pH3_MeCN = predict(MLR_2tR_st61_regressor_pH3_MeCN, newdata = MPs_separated_columns), 
         MLR_2tR_st61_slope_pred_pH5_MeCN = predict(MLR_2tR_st61_regressor_pH5_MeCN, newdata = MPs_separated_columns), 
         MLR_2tR_st61_slope_pred_pH8_MeCN = predict(MLR_2tR_st61_regressor_pH8_MeCN, newdata = MPs_separated_columns),
        )

MLR_2tR_st61_Take_all_MP <- MLR_2tR_st61_Take_all_MP %>%
  select(Compound_name, MLR_2tR_st61_slope_pred_pH3_MeCN, MLR_2tR_st61_slope_pred_pH5_MeCN, MLR_2tR_st61_slope_pred_pH8_MeCN,)

MLR_2tR_st61_Take_all_MP <- gather(MLR_2tR_st61_Take_all_MP, key = MP, value = "MLR_2tR_st61_slope_pred", MLR_2tR_st61_slope_pred_pH3_MeCN, MLR_2tR_st61_slope_pred_pH5_MeCN, MLR_2tR_st61_slope_pred_pH8_MeCN)

MLR_2tR_st61_Take_all_MP <- MLR_2tR_st61_Take_all_MP %>%
  mutate(MP = case_when(
    MP == "MLR_2tR_st61_slope_pred_pH3_MeCN" ~ "2.7_Formic acid_Acetonitrile",
    MP == "MLR_2tR_st61_slope_pred_pH5_MeCN" ~ "5_Formiate_Acetonitrile",
    MP == "MLR_2tR_st61_slope_pred_pH8_MeCN" ~ "8_Bicarbonate_Acetonitrile",
   
  ))


dataset_M <- dataset_M %>%
  left_join(MLR_2tR_st61_Take_all_MP) %>%
  mutate(MLR_2tR_st61_c_pred = Peak_Area_IC / MLR_2tR_st61_slope_pred)



#Hypothesis 4 61 standards----
# k-nearest neighbour model with euclidean distance. Using the slope in one
# of the mobile phases and the retnetion time from one mobile phase.
# k = 1 and 61 standard compounds.

#MP = 2.7 Ac

knn_st61_1tR_pH_2.7 <- dataset_M %>% 
  filter(MP == "2.7_Formic acid_Acetonitrile")

standards <- subset(knn_st61_1tR_pH_2.7, split == TRUE)


standards <- standards %>%
  group_by(Compound_name, MP, slope,) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup()

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



#MP = ph 5 Ac

knn_st61_1tR_pH_5_Ac <- dataset_M %>% 
  filter(MP == "5_Formiate_Acetonitrile")

standards <- subset(knn_st61_1tR_pH_5_Ac, split == TRUE)


standards <- standards %>%
  group_by(Compound_name, MP, slope,) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup()


knn_st61_1tR_slope_pred <- c()
for(i in 1:403){
  knn_st61_1tR_slope_pred <- c(knn_st61_1tR_slope_pred, assigning_closest_cal_comp_slope(knn_st61_1tR_pH_5_Ac[i,]$ret_time, standards))
}




knn_st61_1tR_pH_5_Ac <- data.frame(knn_st61_1tR_pH_5_Ac,knn_st61_1tR_slope_pred)

knn_st61_1tR_pH_5_Ac <- knn_st61_1tR_pH_5_Ac %>%
  mutate(knn_st61_1tR_c_pred = Peak_Area_IC / knn_st61_1tR_slope_pred)


#MP = pH 8 Ac


knn_st61_1tR_pH_8_Ac <- dataset_M %>% 
  filter(MP == "8_Bicarbonate_Acetonitrile")


standards <- subset(knn_st61_1tR_pH_8_Ac, split == TRUE)


standards <- standards %>%
  group_by(Compound_name, MP, slope,) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup()

knn_st61_1tR_slope_pred <- c()
for(i in 1:391){
  knn_st61_1tR_slope_pred <- c(knn_st61_1tR_slope_pred, assigning_closest_cal_comp_slope(knn_st61_1tR_pH_8_Ac[i,]$ret_time, standards))
}




knn_st61_1tR_pH_8_Ac <- data.frame(knn_st61_1tR_pH_8_Ac,knn_st61_1tR_slope_pred)

knn_st61_1tR_pH_8_Ac <- knn_st61_1tR_pH_8_Ac %>%
  mutate(knn_st61_1tR_c_pred = Peak_Area_IC / knn_st61_1tR_slope_pred)





# New Table

dataset_two = rbind(knn_st61_1tR_pH_2.7, knn_st61_1tR_pH_5_Ac, knn_st61_1tR_pH_8_Ac)%>%
  select(knn_st61_1tR_c_pred,knn_st61_1tR_slope_pred,Compound_name, Peak_Area)


dataset_M <- dataset_M %>%
  left_join(dataset_two)


#Hypothesis 4 manual knn 20 standards----
# k-nearest neighbour model with euclidean distance. Using the slope in one
# of the mobile phases and the retnetion time from one mobile phase.
# k = 1 and 20 standard compounds.

# MP = pH 2.7 Ac

knn_st20_1tR_pH_2.7 <- dataset_M %>% 
  filter(MP == "2.7_Formic acid_Acetonitrile")

standards <- knn_st20_1tR_pH_2.7 %>%
  select(Compound_name, ret_time, slope, MP)

standards <- standards %>%
  filter( Compound_name == "Creatinine" | Compound_name == "Theophylline"| 
            Compound_name == "Gabapentin" |Compound_name == "Cotinine" |
            Compound_name == "4-Chloroaniline" |Compound_name == "Hydrochlorothiazide" |
            Compound_name == "Nicotinamide" |Compound_name == "Quinoxaline" | 
            Compound_name == "Monuron" |Compound_name == "Haloperidol" | 
            Compound_name == "Carbazole" | Compound_name == "Valsartan" | 
            Compound_name == "Alachlor" | Compound_name == "Fipronil" |
            Compound_name == "1-Nitropyrene" | Compound_name == "Abietic acid" |
            Compound_name == "Trifluralin" | Compound_name == "Sudan I" | 
            Compound_name == "Sulfapyridine" | Compound_name == "Metazachlor")

standards <- standards %>%
  group_by(Compound_name, MP, slope,) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup()

assigning_closest_cal_comp_slope <- function(xRT, cal_compounds) {
  RF_cal <- cal_compounds %>% slice(which.min(abs(xRT - ret_time))) %>%select(slope)
  unlist(RF_cal)
}

knn_st20_1tR_slope_pred <- c()
for(i in 1:423){
  knn_st20_1tR_slope_pred <- c(knn_st20_1tR_slope_pred, assigning_closest_cal_comp_slope(knn_st20_1tR_pH_2.7[i,]$ret_time, standards))
}




knn_st20_1tR_pH_2.7 <- data.frame(knn_st20_1tR_pH_2.7,knn_st20_1tR_slope_pred)

knn_st20_1tR_pH_2.7 <- knn_st20_1tR_pH_2.7 %>%
  mutate(knn_st20_1tR_c_pred = Peak_Area_IC / knn_st20_1tR_slope_pred)


#MP = pH5 Ac

knn_st20_1tR_pH_5_Ac <- dataset_M %>% 
  filter(MP == "5_Formiate_Acetonitrile")

standards <- knn_st20_1tR_pH_5_Ac %>%
  select(Compound_name, ret_time, slope, MP)

standards <- standards %>%
  filter( Compound_name == "Creatinine" | Compound_name == "Theophylline"| 
            Compound_name == "Gabapentin" |Compound_name == "Cotinine" |
            Compound_name == "4-Chloroaniline" |Compound_name == "Hydrochlorothiazide" |
            Compound_name == "Nicotinamide" |Compound_name == "Quinoxaline" | 
            Compound_name == "Monuron" |Compound_name == "Haloperidol" | 
            Compound_name == "Carbazole" | Compound_name == "Valsartan" | 
            Compound_name == "Alachlor" | Compound_name == "Fipronil" |
            Compound_name == "1-Nitropyrene" | Compound_name == "Abietic acid" |
            Compound_name == "Trifluralin" | Compound_name == "Sudan I" | 
            Compound_name == "Sulfapyridine" | Compound_name == "Metazachlor")

standards <- standards %>%
  group_by(Compound_name, MP, slope,) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup()


knn_st20_1tR_slope_pred <- c()
for(i in 1:403){
  knn_st20_1tR_slope_pred <- c(knn_st20_1tR_slope_pred, assigning_closest_cal_comp_slope(knn_st20_1tR_pH_5_Ac[i,]$ret_time, standards))
}




knn_st20_1tR_pH_5_Ac <- data.frame(knn_st20_1tR_pH_5_Ac,knn_st20_1tR_slope_pred)

knn_st20_1tR_pH_5_Ac <- knn_st20_1tR_pH_5_Ac %>%
  mutate(knn_st20_1tR_c_pred = Peak_Area_IC / knn_st20_1tR_slope_pred)


#MP = pH 8 Ac


knn_st20_1tR_pH_8_Ac <- dataset_M %>% 
  filter(MP == "8_Bicarbonate_Acetonitrile")

standards <- knn_st20_1tR_pH_8_Ac %>%
  select(Compound_name, ret_time, slope, MP)

standards <- standards %>%
  filter( Compound_name == "Creatinine" | Compound_name == "Theophylline"| 
            Compound_name == "Gabapentin" |Compound_name == "Cotinine" |
            Compound_name == "4-Chloroaniline" |Compound_name == "Hydrochlorothiazide" |
            Compound_name == "Nicotinamide" |Compound_name == "Quinoxaline" | 
            Compound_name == "Monuron" |Compound_name == "Haloperidol" | 
            Compound_name == "Carbazole" | Compound_name == "Valsartan" | 
            Compound_name == "Alachlor" | Compound_name == "Fipronil" |
            Compound_name == "1-Nitropyrene" | Compound_name == "Abietic acid" |
            Compound_name == "Trifluralin" | Compound_name == "Sudan I" | 
            Compound_name == "Sulfapyridine" | Compound_name == "Metazachlor")

standards <- standards %>%
  group_by(Compound_name, MP, slope,) %>%
  summarise(ret_time = mean(ret_time)) %>%
  ungroup()


knn_st20_1tR_slope_pred <- c()
for(i in 1:391){
  knn_st20_1tR_slope_pred <- c(knn_st20_1tR_slope_pred, assigning_closest_cal_comp_slope(knn_st20_1tR_pH_8_Ac[i,]$ret_time, standards))
}




knn_st20_1tR_pH_8_Ac <- data.frame(knn_st20_1tR_pH_8_Ac,knn_st20_1tR_slope_pred)

knn_st20_1tR_pH_8_Ac <- knn_st20_1tR_pH_8_Ac %>%
  mutate(knn_st20_1tR_c_pred = Peak_Area_IC / knn_st20_1tR_slope_pred)



# New Table

dataset_two = rbind(knn_st20_1tR_pH_2.7, knn_st20_1tR_pH_5_Ac, knn_st20_1tR_pH_8_Ac)%>%
  select(knn_st20_1tR_c_pred,knn_st20_1tR_slope_pred,Compound_name, Peak_Area)


dataset_M <- dataset_M %>%
  left_join(dataset_two)



#Hypothesis five 3 retention times----
# k-nearest neighbour model with euclidean distance. Using the slope in one
# of the mobile phases and the retnetion times from mobile phases with the 
# same organic phase (3 mobile phases).
# k > 1 and 61 standard compounds.

knn_3tR_st61_regressor_pH3_MeCN <- train(`slope_2.7_Formic acid_Acetonitrile` ~  `2.7_Formic acid_Acetonitrile` + `5_Formiate_Acetonitrile` + `8_Bicarbonate_Acetonitrile`, 
                                  data = MPs_separated_columns_model_development,
                                  method = "knn", 
                                  trControl = fitControl, 
                                  verbose = FALSE)




knn_3tR_st61_regressor_pH5_MeCN <- train(`slope_5_Formiate_Acetonitrile` ~ `2.7_Formic acid_Acetonitrile` + `5_Formiate_Acetonitrile` + `8_Bicarbonate_Acetonitrile` , 
                                   data = MPs_separated_columns_model_development, 
                                   method = "knn", 
                                   trControl = fitControl, 
                                   verbose = FALSE)


knn_3tR_st61_regressor_pH8_MeCN <- train(`slope_8_Bicarbonate_Acetonitrile` ~ `2.7_Formic acid_Acetonitrile` + `5_Formiate_Acetonitrile` + `8_Bicarbonate_Acetonitrile`, 
                                    data = MPs_separated_columns_model_development, 
                                    method = "knn", 
                                    trControl = fitControl, 
                                    verbose = FALSE)




knn_3tR_st61_Take_all_MP <- MPs_separated_columns %>% 
  mutate(split = split, 
         knn_3tR_st61_slope_pred_pH3_MeCN = predict(knn_3tR_st61_regressor_pH3_MeCN, newdata = MPs_separated_columns), 
         knn_3tR_st61_slope_pred_pH5_MeCN = predict(knn_3tR_st61_regressor_pH5_MeCN, newdata = MPs_separated_columns), 
         knn_3tR_st61_slope_pred_pH8_MeCN = predict(knn_3tR_st61_regressor_pH8_MeCN, newdata = MPs_separated_columns))

knn_3tR_st61_Take_all_MP <- knn_3tR_st61_Take_all_MP %>%
  select(Compound_name, knn_3tR_st61_slope_pred_pH3_MeCN, knn_3tR_st61_slope_pred_pH5_MeCN, knn_3tR_st61_slope_pred_pH8_MeCN)

knn_3tR_st61_Take_all_MP <- gather(knn_3tR_st61_Take_all_MP, key = MP, value = "knn_3tR_st61_slope_pred", knn_3tR_st61_slope_pred_pH3_MeCN, knn_3tR_st61_slope_pred_pH5_MeCN, knn_3tR_st61_slope_pred_pH8_MeCN)

knn_3tR_st61_Take_all_MP <- knn_3tR_st61_Take_all_MP %>%
  mutate(MP = case_when(
    MP == "knn_3tR_st61_slope_pred_pH3_MeCN" ~ "2.7_Formic acid_Acetonitrile",
    MP == "knn_3tR_st61_slope_pred_pH5_MeCN" ~ "5_Formiate_Acetonitrile",
    MP == "knn_3tR_st61_slope_pred_pH8_MeCN" ~ "8_Bicarbonate_Acetonitrile"
  ))



dataset_M <- dataset_M %>%
  left_join(knn_3tR_st61_Take_all_MP) %>%
  mutate(knn_3tR_st61_c_pred = Peak_Area_IC / knn_3tR_st61_slope_pred)

#Hypothesis 5 2 MPs 60 standards----
# k-nearest neighbour model with euclidean distance. Using the slope in one
# of the mobile phases and the retnetion times from higheas and lowes pH vaue in
# mobile phases with the same organic phase (2 mobile phases).
# k > 1 and 61 standard compounds.

knn_2tR_st61_regressor_pH3_MeCN <- train(`slope_2.7_Formic acid_Acetonitrile` ~  `2.7_Formic acid_Acetonitrile` +  `8_Bicarbonate_Acetonitrile`, 
                                   data = MPs_separated_columns_model_development,
                                   method = "knn", 
                                   trControl = fitControl, 
                                   verbose = FALSE)




knn_2tR_st61_regressor_pH5_MeCN <- train(`slope_5_Formiate_Acetonitrile` ~ `2.7_Formic acid_Acetonitrile` +  `8_Bicarbonate_Acetonitrile` , 
                                   data = MPs_separated_columns_model_development, 
                                   method = "knn", 
                                   trControl = fitControl, 
                                   verbose = FALSE)




knn_2tR_st61_regressor_pH8_MeCN <- train(`slope_8_Bicarbonate_Acetonitrile` ~ `2.7_Formic acid_Acetonitrile` + `8_Bicarbonate_Acetonitrile`, 
                                    data = MPs_separated_columns_model_development, 
                                    method = "knn", 
                                    trControl = fitControl, 
                                    verbose = FALSE)


knn_2tR_st61_Take_all_MP <- MPs_separated_columns %>% 
  mutate(split = split, 
         knn_2tR_st61_slope_pred_pH3_MeCN = predict(knn_2tR_st61_regressor_pH3_MeCN, newdata = MPs_separated_columns), 
         knn_2tR_st61_slope_pred_pH5_MeCN = predict(knn_2tR_st61_regressor_pH5_MeCN, newdata = MPs_separated_columns), 
         knn_2tR_st61_slope_pred_pH8_MeCN = predict(knn_2tR_st61_regressor_pH8_MeCN, newdata = MPs_separated_columns))
         
knn_2tR_st61_Take_all_MP <- knn_2tR_st61_Take_all_MP %>%
  select(Compound_name, knn_2tR_st61_slope_pred_pH3_MeCN, knn_2tR_st61_slope_pred_pH5_MeCN, knn_2tR_st61_slope_pred_pH8_MeCN)

knn_2tR_st61_Take_all_MP <- gather(knn_2tR_st61_Take_all_MP, key = MP, value = "knn_2tR_st61_slope_pred", knn_2tR_st61_slope_pred_pH3_MeCN, knn_2tR_st61_slope_pred_pH5_MeCN, knn_2tR_st61_slope_pred_pH8_MeCN)

knn_2tR_st61_Take_all_MP <- knn_2tR_st61_Take_all_MP %>%
  mutate(MP = case_when(
    MP == "knn_2tR_st61_slope_pred_pH3_MeCN" ~ "2.7_Formic acid_Acetonitrile",
    MP == "knn_2tR_st61_slope_pred_pH5_MeCN" ~ "5_Formiate_Acetonitrile",
    MP == "knn_2tR_st61_slope_pred_pH8_MeCN" ~ "8_Bicarbonate_Acetonitrile"
  ))



dataset_M <- dataset_M %>%
  left_join(knn_2tR_st61_Take_all_MP) %>%
  mutate(knn_2tR_st61_c_pred = Peak_Area_IC / knn_2tR_st61_slope_pred)


#Hypothesis 6----
# k-nearest neighbour model with cosine similarity. Using slopes from one mobile
# phase and retention times for the highest and lowest pH mobile phases with the
# same organic phase.
# k = 1 and 61 standard compounds.




fn_cosine_knn_k <- function(trainingset_all_variables, compound_of_interest, predicted_variable, list_of_input_variables, k) {
  trainingset <- trainingset_all_variables %>% select(list_of_input_variables)
  scaled.trainingset <- scale(trainingset)
  length_of_training_set <- dim(trainingset)[1]
  Matrix <- bind_rows(trainingset, compound_of_interest %>% select(list_of_input_variables))
  Matrix <- scale(Matrix, attr(scaled.trainingset, "scaled:center"), attr(scaled.trainingset, "scaled:scale"))
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  sim <- as.matrix(sim)
  sim <- sim %*% t(sim)
  D_sim <- (1 - sim)
  similarity_to_compound_of_interest <- D_sim[,(length_of_training_set+1)] #only the dissimilarities to the compound of interest
  similarity_to_compound_of_interest <- similarity_to_compound_of_interest[1:length_of_training_set] #remove similarity to itself
  similarity_to_compound_of_interest <- data.frame(similarity_to_compound_of_interest)
  trainingset_with_similarity <- bind_cols(trainingset_all_variables, similarity_to_compound_of_interest)
  slopes <- c()
  for (i in 1:k) {
    output<- trainingset_with_similarity %>% slice(which.min(similarity_to_compound_of_interest)) %>%select(predicted_variable)
    slopes <- c(slopes, unlist(output))
    #print(slopes)
    trainingset_with_similarity <- trainingset_with_similarity %>% slice(-which.min(similarity_to_compound_of_interest))
  }
  return(mean(slopes))
}



Cosine_slope_pred_Ac_2.7 <- c()
for(i in 1:67){
  Cosine_slope_pred_Ac_2.7 <- c(Cosine_slope_pred_Ac_2.7, fn_cosine_knn_k(MPs_separated_columns_model_development, MPs_separated_columns[i,], 'slope_2.7_Formic acid_Acetonitrile', c('2.7_Formic acid_Acetonitrile', '8_Bicarbonate_Acetonitrile'), 1))
}

Cosine_slope_pred_Ac_5 <- c()
for(i in 1:67){
  Cosine_slope_pred_Ac_5 <- c(Cosine_slope_pred_Ac_5, fn_cosine_knn_k(MPs_separated_columns_model_development, MPs_separated_columns[i,], 'slope_5_Formiate_Acetonitrile', c('2.7_Formic acid_Acetonitrile', '8_Bicarbonate_Acetonitrile'), 1))
}

Cosine_slope_pred_Ac_8 <- c()
for(i in 1:67){
  Cosine_slope_pred_Ac_8 <- c(Cosine_slope_pred_Ac_8, fn_cosine_knn_k(MPs_separated_columns_model_development, MPs_separated_columns[i,], 'slope_8_Bicarbonate_Acetonitrile', c('2.7_Formic acid_Acetonitrile', '8_Bicarbonate_Acetonitrile'), 1))
}


MPs_separated_columns <- data.frame(MPs_separated_columns, Cosine_slope_pred_Ac_2.7, Cosine_slope_pred_Ac_5, Cosine_slope_pred_Ac_8)

Cosine_Take_all_MP_pred <- MPs_separated_columns %>%
  select(Compound_name, Cosine_slope_pred_Ac_2.7, Cosine_slope_pred_Ac_5, Cosine_slope_pred_Ac_8)

Cosine_Take_all_MP_pred <- gather(Cosine_Take_all_MP_pred, key = MP, value = "Cosine_slope_pred", Cosine_slope_pred_Ac_2.7, Cosine_slope_pred_Ac_5, Cosine_slope_pred_Ac_8)

Cosine_Take_all_MP_pred <- Cosine_Take_all_MP_pred %>%
  mutate(MP = case_when(
    MP == "Cosine_slope_pred_Ac_2.7" ~ "2.7_Formic acid_Acetonitrile",
    MP == "Cosine_slope_pred_Ac_5" ~ "5_Formiate_Acetonitrile",
    MP == "Cosine_slope_pred_Ac_8" ~ "8_Bicarbonate_Acetonitrile"
  ))

dataset_M <- dataset_M %>%
  left_join(Cosine_Take_all_MP_pred) %>%
  mutate(Cosine_c_pred = Peak_Area_IC / Cosine_slope_pred)



dataset_M_2 <- dataset_M %>%
  select(Compound_name, MP, ret_time, slope) %>%
  group_by(Compound_name, MP) %>%
  summarise(ret_time = mean(ret_time),
            slope = mean(slope)) %>%
  ungroup()

#Creating the H_two_Take_all_MP dataset
dataset_M_2_A <- dataset_M_2 %>%
  select(-slope)
dataset_M_2_B <- dataset_M_2 %>%
  select(-ret_time)

dataset_M_2_a <-  spread(dataset_M_2_A, key = "MP", value = "ret_time")
dataset_M_2_c <- dataset_M_2_B %>%
  mutate(MP = str_c("slope_", MP, sep = ""))
dataset_M_2_b <- spread(dataset_M_2_c, key = "MP", value = "slope")

dataset_M_2 <- dataset_M_2_a %>%
  left_join(dataset_M_2_b)%>%
  drop_na() %>%
  unique() 

dataset_M <- dataset_M %>%
  filter(Linear_6tR_st_61_c_pred > 0,
         MLR_3tR_st61_c_pred > 0,
         MLR_2tR_st61_c_pred > 0)


###adding now the comparisons to the dataset
dataset_M <- dataset_M %>%
  mutate(mean_c_comp = case_when(mean_c < conc_M ~ conc_M / mean_c, 
                                 mean_c > conc_M ~ mean_c /conc_M ),
         median_c_comp = case_when(median_c < conc_M ~ conc_M / median_c, 
                                   median_c > conc_M ~ median_c /conc_M ),
         Linear_6tR_st_61_c_comp = case_when(Linear_6tR_st_61_c_pred < conc_M ~ conc_M / Linear_6tR_st_61_c_pred, 
                                  Linear_6tR_st_61_c_pred > conc_M ~ Linear_6tR_st_61_c_pred /conc_M),
         MLR_3tR_st61_c_comp = case_when(MLR_3tR_st61_c_pred < conc_M ~ conc_M / MLR_3tR_st61_c_pred, 
                                         MLR_3tR_st61_c_pred > conc_M ~ MLR_3tR_st61_c_pred /conc_M ),
         MLR_2tR_st61_c_comp = case_when(MLR_2tR_st61_c_pred < conc_M ~ conc_M / MLR_2tR_st61_c_pred, 
                                         MLR_2tR_st61_c_pred > conc_M ~ MLR_2tR_st61_c_pred /conc_M ),
         knn_3tR_st61_c_comp = case_when(knn_3tR_st61_c_pred < conc_M ~ conc_M / knn_3tR_st61_c_pred, 
                                        knn_3tR_st61_c_pred > conc_M ~ knn_3tR_st61_c_pred /conc_M ),
         knn_2tR_st61_c_comp = case_when(knn_2tR_st61_c_pred < conc_M ~ conc_M / knn_2tR_st61_c_pred, 
                                        knn_2tR_st61_c_pred > conc_M ~ knn_2tR_st61_c_pred /conc_M ),
         knn_st61_1tR_c_comp = case_when(knn_st61_1tR_c_pred < conc_M ~ conc_M / knn_st61_1tR_c_pred, 
                                      knn_st61_1tR_c_pred > conc_M ~ knn_st61_1tR_c_pred /conc_M ),
         knn_st20_1tR_c_comp = case_when(knn_st20_1tR_c_pred < conc_M ~ conc_M / knn_st20_1tR_c_pred, 
                                   knn_st20_1tR_c_pred > conc_M ~ knn_st20_1tR_c_pred /conc_M ),
         Cosine_c_comp = case_when(Cosine_c_pred < conc_M ~ conc_M / Cosine_c_pred, 
                                  Cosine_c_pred > conc_M ~ Cosine_c_pred /conc_M ))

all_concentrations <- dataset_M %>%
  select(Compound_name, mean_c, median_c, Linear_6tR_st_61_c_pred, MLR_3tR_st61_c_pred, MLR_2tR_st61_c_pred, knn_st61_1tR_c_pred, knn_st20_1tR_c_pred, knn_3tR_st61_c_pred, knn_2tR_st61_c_pred, Cosine_c_pred, conc_M)

write_delim(all_concentrations, "all concentrations mix A.csv", delim = ";")

dataset_M <- dataset_M %>%
  mutate(split_two = case_when(  Compound_name == "Creatinine" ~TRUE,
                                 Compound_name == "Theophylline"~TRUE,
                                 Compound_name == "Gabapentin" ~TRUE, 
                                 Compound_name == "Cotinine" ~TRUE,
                                 Compound_name == "4-Chloroaniline" ~TRUE,
                                 Compound_name == "Hydrochlorothiazide" ~TRUE,
                                 Compound_name == "Nicotinamide" ~TRUE,
                                 Compound_name == "Quinoxaline" ~TRUE,
                                 Compound_name == "Monuron" ~TRUE,
                                 Compound_name == "Haloperidol" ~TRUE, 
                                 Compound_name == "Carbazole" ~TRUE,
                                 Compound_name == "Valsartan" ~TRUE,
                                 Compound_name == "Alachlor" ~TRUE,
                                 Compound_name == "Fipronil" ~TRUE,
                                 Compound_name == "1-Nitropyrene" ~TRUE,
                                 Compound_name == "Abietic acid" ~TRUE,
                                 Compound_name == "Trifluralin" ~TRUE,
                                 Compound_name == "Sudan I" ~TRUE,
                                 Compound_name == "Sulfapyridine" ~TRUE,
                                 Compound_name == "Metazachlor" ~TRUE, 
                                 TRUE ~ FALSE))




all_c_comp <- dataset_M %>%
  select(Compound_name, mean_c_comp, median_c_comp, Linear_6tR_st_61_c_comp, MLR_3tR_st61_c_comp, MLR_2tR_st61_c_comp, knn_3tR_st61_c_comp, knn_2tR_st61_c_comp, knn_st61_1tR_c_comp,Cosine_c_comp, MP, split, slope)

all_c_comp <- gather(data = all_c_comp, key = "models", value = error, mean_c_comp, median_c_comp, Linear_6tR_st_61_c_comp, MLR_3tR_st61_c_comp, MLR_2tR_st61_c_comp, knn_3tR_st61_c_comp, knn_2tR_st61_c_comp,  knn_st61_1tR_c_comp, Cosine_c_comp)
all_c_comp_c <- all_c_comp %>%
  filter(split == FALSE)

all_c_comp_c <- all_c_comp_c %>%
  mutate(model = case_when(models == "median_c_comp" ~ "Median slope",
                           models == "mean_c_comp" ~ "Mean slope",
                           models == "Linear_6tR_st_61_c_comp" ~ "Linear regression",
                           models == "MLR_2tR_st61_c_comp" ~ "Multilinear regression (2 MPs)",
                           models == "MLR_3tR_st61_c_comp" ~ "Multilinear regression (3 MPs)",
                           models == "Cosine_c_comp" ~ "Cosine similarity",
                           models == "knn_st61_1tR_c_comp" ~ "k-nearest neighbour (k = 1, 1 MP)",
                           models == "knn_2tR_st61_c_comp" ~ "k-nearest neighbour (k > 1, 2 MPs)",
                           models == "knn_3tR_st61_c_comp" ~ "k-nearest neighbour (k > 1, 3 MPs)"))

all_c_comp_c<- all_c_comp_c %>%
  mutate(mobile_phase = case_when(MP == "2.7_Formic acid_Acetonitrile" ~ "pH 2.7",
                                  MP == "5_Formiate_Acetonitrile" ~ "pH 5",
                                  MP == "8_Bicarbonate_Acetonitrile" ~ "pH 8"))
all_c_comp_d <- all_c_comp_c %>%
  filter(model != "k-nearest neighbour (k > 1, 3 MPs)",
         model != "Multilinear regression (3 MPs)")

#write_delim(all_c_comp_d, "Mix_M_all_c_comp.csv", delim = ";")


less_than_ten <- all_c_comp %>%
  drop_na()%>%
  group_by(models, MP) %>%
  summarise(
    mean_c_comp =mean(error),
    median_c_comp = median(error),
    max_c_comp = max(error),
    quantile_c_comp = quantile(error,probs = c(0.95)),
    n_dp = length(error),
    n_dp_less_than_ten = length(error[error<10 & error>1]),
    percentage_less_than_ten =(length(error[error<10 & error>1]))/(length(error)))%>%
  ungroup()



less_than_ten_FALSE <- all_c_comp_c %>%
  drop_na()%>%
  group_by(models, MP) %>%
  summarise(
    mean_c_comp =mean(error),
    median_c_comp = median(error),
    max_c_comp = max(error),
    quantile_c_comp = quantile(error,probs = c(0.95)),
    n_dp = length(error),
    n_dp_less_than_ten = length(error[error<10 & error>1]),
    percentage_less_than_ten =(length(error[error<10 & error>1]))/(length(error)))%>%
  ungroup()


all_c_comp_a <- dataset_M %>%
  select(Compound_name, knn_st20_1tR_c_comp, MP, split_two)

all_c_comp_a <- gather(data = all_c_comp_a, key = "models", value = error, knn_st20_1tR_c_comp)


all_c_comp_e <- all_c_comp_d%>%
  filter(MP == "2.7_Formic acid_Acetonitrile")%>%
  mutate(mobile_phase = case_when(mobile_phase == "pH 2.7" ~ "pH 2.7 Mix A" ))

less_than_ten_a = all_c_comp_a %>%
  drop_na()%>%
  group_by(models, MP) %>%
  summarise(
    mean_c_comp =mean(error),
    median_c_comp = median(error),
    max_c_comp = max(error),
    quantile_c_comp = quantile(error,probs = c(0.95)),
    n_dp = length(error),
    n_dp_less_than_ten = length(error[error<10 & error>1]),
    percentage_less_than_ten =(length(error[error<10 & error>1]))/(length(error)))%>%
  ungroup()

all_c_comp_b <- all_c_comp_a %>%
  filter(split_two == FALSE)

less_than_ten_a_FALSE = all_c_comp_b %>%
  drop_na()%>%
  group_by(models, MP) %>%
  summarise(
    mean_c_comp =mean(error),
    median_c_comp = median(error),
    max_c_comp = max(error),
    quantile_c_comp = quantile(error,probs = c(0.95)),
    n_dp = length(error),
    n_dp_less_than_ten = length(error[error<10 & error>1]),
    percentage_less_than_ten =(length(error[error<10 & error>1]))/(length(error)))%>%
  ungroup()


#---------------------------THEME---------------------------------------------------------------
font <- choose_font("calibri")
fontsize <- 12
basecolor <- "#717171" #your base color for both font and lines


my_theme <-   theme(
  #remove the background of the plot
  plot.background = element_blank(),
  #and from the panel as well
  panel.background = element_blank(),
  #define the width and color of the axis on the plot
  axis.line = element_line(size = 0.5,
                           color = basecolor),
  #if you use plot title you can specify parameters here
  #PS! use plot title only if you send or show the plot on its own 
  #for plots on the slide/thesis use slide title and figure caption 
  plot.title = element_text(color = basecolor,
                            size = 14,
                            face = "bold"),
  #specify the size and style of the text on the plot, e.g. axis title
  text = element_text(family = font,
                      size = fontsize,
                      color = basecolor),
  legend.key = element_blank(),
  #to remove or adjust the position of the legend
  #"none" - is no legend; "top" "bottom", "right", "left";
  #or by coordinates. 
  #c(0, 0) corresponds to the "bottom left" 
  #and c(1, 1) corresponds to the "top right" position.
  legend.position = c(0.8,0.2),
  #if you have a legend title and text you can specify font size here
  #here it indicates no legend title
  legend.title = element_text(family = font,
                              size = fontsize,
                              color = basecolor),
  legend.text = element_text(family = font,
                             size = fontsize,
                             color = basecolor),
  #specify axis marks text
  axis.text = element_text(family = font,
                           size = fontsize,
                           color = basecolor),
  #remove tick marks
  axis.ticks = element_blank(),
  #define the ratio of x and y axis
  #PS! for scatter plots it needs to be 1!
  #for predicted - measured plots also adjust the ranges!
  aspect.ratio = 1,
  #adjust the position of the axis title
  axis.title.x = element_text(hjust = c(1), vjust = c(0)),
  axis.title.y = element_text(hjust = c(1), vjust = c(1))
)



#--------------------------------------------------------------------------------------------------




ggplot(data = all_c_comp) +
  geom_boxplot(mapping = (aes(x = models, y = error)))+
  scale_color_manual(values = c("#222d50" , "#007d51"))+
  scale_y_log10()+
  coord_flip()+
  facet_wrap(~MP)





box <- ggplot(data = all_c_comp_d) +
  geom_boxplot(mapping = (aes(x = model, y = error)))+
  scale_color_manual(values = c("#222d50" , "#007d51"))+
  scale_y_log10()+
  coord_flip()+
  xlab("")+
  ylab("Error")+
  facet_wrap(~mobile_phase)

box+
  my_theme+
  theme(
    axis.text.x = element_text(angle = 90,
                               hjust = 1))+
  theme(strip.background = element_rect(fill = "#2A8CB7"),
        aspect.ratio = 1.2)


#ggsave("Box_Plot_M.svg", width = 16, height = 8, units = "cm")


box_3 <- ggplot(data = all_c_comp_e) +
  geom_boxplot(mapping = (aes(x = model, y = error)))+
  scale_color_manual(values = c("#222d50" , "#007d51"))+
  scale_y_log10()+
  coord_flip()+
  xlab("")+
  ylab("")+
  facet_wrap(~mobile_phase)


box_3+
  my_theme+
  theme(
    axis.text.x = element_text(angle = 90,
                               hjust = 1))+
  theme(strip.background = element_rect(fill = "#2A8CB7"),
        aspect.ratio = 1.2)


#ggsave("Box_Plot_3_SB.svg", width = 8, height = 8, units = "cm")


box_2 <- ggplot(data = all_c_comp) +
  geom_boxplot(mapping = (aes(x = model, y = error)))+
  scale_color_manual(values = c("#222d50" , "#007d51"))+
  scale_y_log10()+
  coord_flip()+
  xlab("")+
  ylab("Error")+
  facet_wrap(~mobile_phase)

box_2+
  my_theme+
  theme(
    axis.text.x = element_text(angle = 90,
                               hjust = 1))+
  theme(strip.background = element_rect(fill = "#2A8CB7"),
        aspect.ratio = 1.5)


#ggsave("Box_Plot_M_2.svg", width = 16, height = 8, units = "cm")



mid <- mean(log10(dataset_M_2$`slope_2.7_Formic acid_Acetonitrile`)%>% na.omit())
ggplotly(
  geom_point(mapping = aes(x = `2.7_Formic acid_Acetonitrile`, y = `8_Bicarbonate_Acetonitrile`, colour = log10(`slope_2.7_Formic acid_Acetonitrile`), text = Compound_name), size = 2, alpha = 0.7) +
    scale_color_gradient2(midpoint = mid, low = "#2A8CB7", mid = "#EBEBEB", high = "#F62DAE", space = "rgb")+
    xlim(0,25)+
    ylim(0,25)+
    labs(color = "Slope",
         x = "retention time pH 2.7 (min)",
         y = "retention time pH 8 (min)") +
    theme(aspect.ratio = 0.1,)
)
rt2 <- ggplot(data = dataset_M_2)+
  geom_point(mapping = aes(x = `2.7_Formic acid_Acetonitrile`, y = `8_Bicarbonate_Acetonitrile`, colour = log10(`slope_2.7_Formic acid_Acetonitrile`), text = Compound_name), size = 2, alpha = 0.7) +
  scale_color_gradient2(midpoint = mid, low = "#2A8CB7", mid = "#EBEBEB", high = "#F62DAE", space = "rgb")+
  xlim(0,25)+
  ylim(0,25)+
  labs(color = "Slope",
       x = "retention time pH 2.7 (min)",
       y = "retention time pH 8 (min)") +
  theme(aspect.ratio = 0.1,)

rt2+
  my_theme+
  theme(
    #legend.key.size = unit(0.4, "cm"),
    # legend.key.height = unit(0.3, "cm"),
    # legend.key.width = unit(0.5, "cm"),
    legend.direction = "horizontal",
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    legend.position = "bottom")+
  ggtitle("Mix A")


#ggsave("rt2_M.svg", width = 8, height = 8, units = "cm")

histogram <- ggplot(data = dataset_M_2)+
  geom_histogram(mapping = aes(x = log(dataset_M_2$`slope_2.7_Formic acid_Acetonitrile`)))+
  xlab("slope")
histogram+
  my_theme

#ggsave("histogram_M.svg")

ggplotly(
ggplot(data = all_c_comp_c)+
  geom_point(mapping = aes(x = slope, y = error, text = Compound_name)))
error <- ggplot(data = all_c_comp_c)+
  geom_point(mapping = aes(x = slope, y = error), colour = "#2A8CB7")+
  xlab("Slope")+
  ylab("Error")
error+
  my_theme+
  ggtitle("Mix A")+
  theme(
    axis.text.x = element_text(angle = 90,
                               hjust = 1))
#ggsave("error_M.svg", width = 8, height = 8, units = "cm")

ggplotly(
  ggplot(data = MPs_separated_columns)+
    geom_point(mapping = aes(x = MPs_separated_columns$X2.7_Formic.acid_Acetonitrile, y = MPs_separated_columns$slope_2.7_Formic.acid_Acetonitrile, text = Compound_name))+
    geom_line(mapping = aes(x = MPs_separated_columns$X2.7_Formic.acid_Acetonitrile, y =  MPs_separated_columns$slope_2.7_Formic.acid_Acetonitrile))+
    xlab("retention time")+
    ylab("slope")+
    scale_y_log10()
)

rt_slope <-ggplot(data = MPs_separated_columns)+
  #geom_point(mapping = aes(x = MPs_separated_columns$X2.7_Formic.acid_Acetonitrile, y =  MPs_separated_columns$slope_2.7_Formic.acid_Acetonitrile, text = Compound_name))+
  geom_line(mapping = aes(x = MPs_separated_columns$X2.7_Formic.acid_Acetonitrile, y =  MPs_separated_columns$slope_2.7_Formic.acid_Acetonitrile), colour = "#2A8CB7", size = 0.5)+
              xlab("retention time (min)")+
              ylab("slope")+
              scale_y_log10()
            

rt_slope+
  my_theme+
  ggtitle("Mix A")
#ggsave("rt_slope_M.svg", width = 5, height = 5, units = "cm")






#--------------------------Tables for thesis---------------------------------------------

all_c_comp_2.7_Ac<- all_c_comp_c %>%
  filter(MP == "2.7_Formic acid_Acetonitrile")
  #filter(models != "MLR_3tR_st61_c_comp", models != "knn_3tR_st61_c_comp", models != "knn_2tR_st61_c_comp")

Model_comparison_2.7_Ac = all_c_comp_2.7_Ac %>%
  drop_na()%>%
  group_by(models, MP) %>%
  summarise(
    mean_c_comp =mean(error),
    median_c_comp = median(error),
    max_c_comp = max(error),
    quantile_c_comp = quantile(error,probs = c(0.95)),
    percentage_less_than_ten =(length(error[error<10 & error>1]))/(length(error)))%>%
  ungroup()

#write_delim(Model_comparison_2.7_Ac, "Mix_M_evaluation_parameters.csv", delim = ";")

all_c_comp_H_five_M1 <- all_c_comp %>%
  filter(MP == "2.7_Formic acid_Acetonitrile")%>%
  filter(models == "H_five_M1_c_comp" | models == "knn_3tR_st61_c_comp" |models == "knn_2tR_st61_c_comp")

Model_comparison_H_five_M1 = all_c_comp_H_five_M1 %>%
  drop_na()%>%
  group_by(models, MP) %>%
  summarise(
    mean_c_comp =mean(error),
    median_c_comp = median(error),
    max_c_comp = max(error),
    quantile_c_comp = quantile(error,probs = c(0.95)),
    percentage_less_than_ten =(length(error[error<10 & error>1]))/(length(error)))%>%
  ungroup()

all_c_comp_H_three_M2 <- all_c_comp %>%
  filter(models == "knn_2tR_st61_c_comp")

Model_comparison_H_three_M2_all_MP = all_c_comp_H_three_M2 %>%
  drop_na()%>%
  group_by(models, MP) %>%
  summarise(
    mean_c_comp =mean(error),
    median_c_comp = median(error),
    max_c_comp = max(error),
    quantile_c_comp = quantile(error,probs = c(0.95)),
    percentage_less_than_ten =(length(error[error<10 & error>1]))/(length(error)))%>%
  ungroup()

