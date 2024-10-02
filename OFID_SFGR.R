library(tidyverse)
library(ggplot2)
library(lubridate)
library(ranger)
library(pROC)
library(ROCR)
library(gtsummary)
library(flextable)
library(vip)
library(randomForest)
library(tuneRanger)
library(missRanger)
library(gtsummary)
library(skimr)
library(caret)
library(zoo)
library(slider)
library(purrr)  # For mapping and working with lists
library(pROC)   
library(viridis)
library(fuzzyjoin)
library(permimp)
setwd("/Users/rjwilliams/Documents/Rickettsia")

modis <- read.csv('MODIS_External_validation.csv')
modis <- modis[,-c(1)]
modis$Date <- as.Date(modis$Date)
modis$Next_Date<- as.Date(modis$Next_Date)

# Load the necessary library
library(lubridate)

IS<-read.csv('OFID_SFGR_nomal.csv')
colnames(IS)
IS <- IS[,-c(1)]
IS[,1:6] <- lapply(IS[,1:6], factor)
IS$heartrate <- as.numeric(IS$heartrate)
IS$SBP <- as.numeric(IS$SBP)
IS$DBP <- as.numeric(IS$DBP)
IS[,24:25] <- lapply(IS[,24:25], factor)
IS[,27:35] <- lapply(IS[,27:35], factor)

# Convert the dateadmission column to character
IS$dateadmission <- as.character(IS$dateadmission)

# Add leading zeros to dates that are less than 6 digits
IS$dateadmission <- ifelse(nchar(IS$dateadmission) < 6, 
                            paste0("0", IS$dateadmission), 
                            IS$dateadmission)

# Convert the dateadmission column to Date format
IS$dateadmission <- mdy(paste0(substr(IS$dateadmission, 1, 2), "-", 
                                substr(IS$dateadmission, 3, 4), "-", 
                                substr(IS$dateadmission, 5, 8)))




#Lag 1month
IS$lag1 <- IS$dateadmission%m-% months(1) 
#Lag 2month
IS$lag2 <- IS$dateadmission%m-% months(2)
#Lag 3month
IS$lag3 <- IS$dateadmission%m-% months(3)



modis <- modis %>%
  mutate(PrevDate = lead(Date, order_by = Date) - 1)

modis <- modis %>% dplyr::select(-c(Next_Date))

modis$Date <- as.Date(modis$Date)


modis <- modis %>%
  mutate(PrevDate = lead(Date, order_by = Date) - 1)



result <- fuzzy_left_join(
  IS, modis,
  by = c(
    'lag1'= 'Date',
    'lag1' = 'PrevDate'
  ),
  match_fun = list(`>=`, `<=`))

result <- result %>% dplyr::select(-c('PrevDate', 'Date'))

result <- result %>%
  rename_with(~ paste0(., "_lag1"), .cols = setdiff(names(.), names(IS)))


modis_lag2 <- modis %>% rename_with(~ paste0(., "_lag2"))

result <- fuzzy_left_join(
  result, modis_lag2,
  by = c(
    'lag2'= 'Date_lag2',
    'lag2' = 'PrevDate_lag2'
  ),
  match_fun = list(`>=`, `<=`))

result <- result %>% dplyr::select(-c('PrevDate_lag2', 'Date_lag2'))

modis_lag3 <- modis %>% rename_with(~ paste0(., "_lag3"))

result <- fuzzy_left_join(
  result, modis_lag3,
  by = c(
    'lag3'= 'Date_lag3',
    'lag3' = 'PrevDate_lag3'
  ),
  match_fun = list(`>=`, `<=`))

result <- result %>% dplyr::select(-c('PrevDate_lag3', 'Date_lag3', 'lag1', 'lag2', 'lag3'))

logical_columns <- sapply(result, is.logical)  # Find logical columns
result[logical_columns] <- lapply(result[logical_columns], function(x) as.factor(as.integer(x)))



### TOP 10 variables from satellite

#JUST ENVIRONMENTAL PREDICTORS
sf_new <- result %>% dplyr::select(c('sfgr', 'Mean__1_km_16_days_EVI_lag1':'Daytime_Minimum_lag3'))
result <- result %>% dplyr::select(-c(dateadmission))
rf_model_new <- randomForest(sfgr~., data = sf_new, ntree = 1000, keep.forest = TRUE, keep.inbag = TRUE)
perm_imp <- permimp(rf_model_new, conditional=T,
                    data = sf_new,AUC=T)

pid <- data.frame(perm_imp$values)
names <- rownames(pid)
pidf <- data.frame(cbind(names, pid))
pidf <- pidf %>% arrange(desc(perm_imp.values))
satvip<- pidf[1:10,1]

satvip

### TOP 10 variables from clinical
colnames(result)


###JUST CLINICAL PREDICTORS
sf_new <- result %>% dplyr::select(c('sfgr':'BMI'))

### SUBSETTING if IMBALANCED
sf_new <- sf_new %>% dplyr::select(-c(malres))
colnames(sf_new)
sf_new <- sf_new %>% dplyr::select(-c(hematochezia,fever1,meningeal))

rf_model_new <- randomForest(sfgr ~ ., data = sf_new, ntree = 1000, keep.forest = TRUE, keep.inbag = TRUE)
perm_imp <- permimp(rf_model_new, conditional=T,
                    data = sf_new,AUC=T)

skim(sf_new)
pid <- data.frame(perm_imp$values)
names <- rownames(pid)
pidf <- data.frame(cbind(names, pid))
pidf <- pidf %>% arrange(desc(perm_imp.values))
clinvip<- pidf[1:10,1]



#########COMNBINED VIP MODEL

colnames(result)
sf_new <- result[,c('sfgr',clinvip, satvip)]



rf_model_new <- randomForest(sfgr ~ ., data = sf_new, ntree = 1000, keep.forest = TRUE, keep.inbag = TRUE)
perm_imp <- permimp(rf_model_new, conditional=T,
                    data = sf_new,AUC=T)

pid <- data.frame(perm_imp$values)
names <- rownames(pid)
pidf <- data.frame(cbind(names, pid))
pidf <- pidf %>% arrange(desc(perm_imp.values))

pidf$perm_imp.values <- round(pidf$perm_imp.values, digits=5)


############PERFORMANCE WITH 10/10##############
#############################################
############PERFORMANCE WITH 10/10##############
#############################################
############PERFORMANCE WITH 10/10##############
#############################################
############PERFORMANCE WITH 10/10##############
#############################################
############PERFORMANCE WITH 10/10##############
#############################################

sf_new2 <- sf_new %>% dplyr::select(-c(Minimum_NDWI_lag2:Maximum__1_km_16_days_NDVI_lag1))



N = dim(sf_c)[1]
str(sf_new)

results_df <- data.frame(
  i = integer(0),
  CImin = double(0),
  AUC = double(0),
  CImax = double(0)
)

# Create a list of data frames to store CIs for each 'i'
results_list <- list()
i_values <- c( 20,15,10,9,8,7,6,5,4,3,2)

results_list <- list()

for (each in 1:10) {
  while (TRUE) {
    trn_idx = sample(1:N, .8*N)
    train=sf_new[trn_idx,]
    test=sf_new[-trn_idx,]
    
    if (sum(test$sfgr == 1) > 0) {
      # At least one instance with qf == 1 found in the test set
      break
    }
  }
  
  rf_model_original <- randomForest(sfgr ~ ., data = train, ntree = 1000, keep.forest = TRUE, keep.inbag = TRUE)
  rf_model_new <- randomForest(sfgr ~ ., data = train, ntree = 1000, keep.forest = TRUE, keep.inbag = TRUE)
  
  
  # Your code that calls permimp
  perm_imp <- permimp(rf_model_new, forest = rf_model_original$forest, inbag = rf_model_original$inbag, 
                      data = train, do_check = FALSE)
  
  pid <- data.frame(perm_imp$values)
  names <- rownames(pid)
  pidf <- data.frame(cbind(names, pid))
  pidf <- pidf %>% arrange(desc(perm_imp.values))
  
  for (i in i_values) {
    # Extract the desired features based on 'i'
    pidf_i <- pidf[1:i, ]
    train1_i <- train[, c(pidf_i$names, 'sfgr')]
    test1_i <- test[, c(pidf_i$names, 'sfgr')]
    
    # Fit the model and calculate CI for this 'i'
    mod <- ranger(sfgr ~ ., data = train1_i, num.trees = 1000, probability = TRUE)
    pred_train1 <- predict(mod, train1_i)$predictions[, 2][which(train1_i$sfgr == 1)]
    preds_test_i <- predict(mod, test1_i)$predictions[, 2]
    CI_i <- as.numeric(ci(roc(test1_i$sfgr, preds_test_i, plot = FALSE)))
    
    # Store the CI values for this 'i' in the data frame
    results_df <- rbind(results_df, data.frame(i = i, CImin = CI_i[1], AUC = CI_i[2], CImax = CI_i[3]))
  }
}

averages <- results_df %>%
  group_by(i) %>%
  summarize(
    Avg_CImin = mean(CImin),
    Avg_CImax = mean(CImax),
    Avg_AUC = mean(AUC)
  )


##### RF MODEL WITHOUT VIPS
##### RF MODEL WITHOUT VIPS
##### RF MODEL WITHOUT VIPS
##### RF MODEL WITHOUT VIPS
##### RF MODEL WITHOUT VIPS
##### RF MODEL WITHOUT VIPS

sf_new2 <- result[,c('sfgr',clinvip)]
                    
N=dim(sf_new2)[1]
AUC=c()
CI=c()
for (each in 1:100){
  print(each)
  trn_idx = sample(1:N,.8*N)
  train=sf_new2[trn_idx,]
  test=sf_new2[-trn_idx,]
  mod= ranger(sfgr~.,data=train, num.trees = 1000,probability=T)
  pred_train1=predict(mod,train)$predictions[,2][which(train$sfgr==1)]
  pred_train0=predict(mod,train)$predictions[,2][which(train$sfgr==0)]
  preds_test=predict(mod,test)$predictions[,2]
  AUC =rbind(AUC,c(as.numeric(roc(test$sfgr,preds_test)$auc)))
  CI =rbind(CI, c(as.numeric(ci(roc(test$sfgr,preds_test)))))
}



CImin=mean(CI[,1])
AUC=mean(CI[,2])
CImax=mean(CI[,3])




#### ROC TEST #####
#### ROC TEST #####
#### ROC TEST #####
#### ROC TEST #####
#### ROC TEST #####
#### ROC TEST #####
#### ROC TEST #####
#### ROC TEST #####
#### ROC TEST #####
#### ROC TEST #####
#### ROC TEST #####

N=dim(sf_new)[1]
p.value=c()
for (each in 1:100){
  print(each)
  trn_idx = sample(1:N,.8*N)
  train=sf_new[trn_idx,]
  test=sf_new[-trn_idx,]
  mod= ranger(sfgr~.,data=train, num.trees = 1000,probability=T)
  pred_train1=predict(mod,train)$predictions[,2][which(train$sfgr==1)]
  pred_train0=predict(mod,train)$predictions[,2][which(train$sfgr==0)]
  preds_test=predict(mod,test)$predictions[,2]
  roc1 = roc(test$sfgr,preds_test)
  
  train2=sf_new2[trn_idx,]
  test2=sf_new2[-trn_idx,]
  mod2= ranger(sfgr~.,data=train2, num.trees = 1000,probability=T)
  pred_train12=predict(mod2,train2)$predictions[,2][which(train2$sfgr==1)]
  pred_train02=predict(mod2,train2)$predictions[,2][which(train2$sfgr==0)]
  preds_test2=predict(mod2,test2)$predictions[,2]
  roc2 = roc(test2$sfgr,preds_test2)
  
  
  
  rtest <- roc.test(roc1, roc2, method = 'bootstrap')
  p.value = rbind(p.value, c(as.numeric(rtest$p.value)))
}

mean(p.value)
median(p.value)

sum(p.value < .05)



############ MODEL PERFORMANCE
############ MODEL PERFORMANCE
############ MODEL PERFORMANCE
############ MODEL PERFORMANCE
############ MODEL PERFORMANCE
############ MODEL PERFORMANCE
############ MODEL PERFORMANCE
############ MODEL PERFORMANCE
library(epiR)

sensitivity_vec <- numeric(100)
specificity_vec <- numeric(100)
ppv_vec <- numeric(100)
npv_vec <- numeric(100)

N = dim(sf_new2)[1]
CI = c()
for (each in 1:100) {
  trn_idx = sample(1:N, .8 * N)
  train = sf_new2[trn_idx, ]
  test = sf_new2[-trn_idx, ]
  mod = ranger(sfgr ~ ., data = train, num.trees = 1000, probability = TRUE)
  predictions <- predict(mod, data = test)
  risk_scores <- predictions$predictions[, 2]
  roc_model <- roc(test$sfgr, risk_scores)
  youdbest <- coords(roc_model, "best", ret = "threshold", best.method = "youden")
  predictions <- ifelse(risk_scores > youdbest[[1]], 1, 0)
  
  # Create a confusion matrix
  cM <- table(predictions, test$sfgr)
  
  
  p <- epi.tests(cM)
  x<- p$detail
  
  sensitivity_vec[each] <- x[x$statistic == 'se', 'est']
  specificity_vec[each] <- x[x$statistic == 'sp', 'est']
  ppv_vec[each] <- x[x$statistic == 'pv.pos', 'est']
  npv_vec[each] <- x[x$statistic == 'pv.neg', 'est']
}

avg_sensitivity <- mean(sensitivity_vec)
avg_specificity <- mean(specificity_vec)
avg_ppv <- mean(ppv_vec)
avg_npv <- mean(npv_vec)

#######################
######################
#####################
#### RESULTS WITH SET SENS AT 70,80,90 AND VICE VERSA
#########


tten<- pidf[1:10,1]

sf_new2 <- result[,c('sfgr',tten)]
results_all = c()
s70 <- c()
s80<- c()
s90 <- c()

for (each in 1:100) {
  trn_idx <- sample(1:N, 0.8 * N)
  train <- sf_new2[trn_idx, ]
  test <- sf_new2[-trn_idx, ]
  
  mod <- ranger(sfgr ~ ., data = train, num.trees = 1000, probability = TRUE)
  predictions <- predict(mod, data = test)
  risk_scores <- predictions$predictions[, 2]
  
  # Fit the ROC model
  roc_model <- roc(test$sfgr, risk_scores)
  
  yy <- data.frame(sens = roc_model$sensitivities, spec = roc_model$specificities) 
  
  target_sensitivities <- c(0.70, 0.80, 0.90)
  
  
  sens_indices <- sapply(target_sensitivities, function(target_sens) which.min(abs(yy$sens - target_sens)))
  
  # Extract the sensitivity and specificity values at the identified indices
  sens_values <- yy$sens[sens_indices]
  
  t7 <- sens_values[1]
  t8 <- sens_values[2]
  t9 <- sens_values[3]
  
  
  s7 <- yy %>% filter(yy$sens == t7)
  spec7 <- mean(s7$spec)
  s8 <- yy %>% filter(yy$sens == t8)
  spec8 <- mean(s8$spec)
  s9 <- yy %>% filter(yy$sens == t9)
  spec9 <- mean(s9$spec)
  
  s70 <- rbind(spec7, s70)
  s80 <- rbind(spec8, s80)
  s90 <- rbind(spec9, s90)
}
mean(s70)
mean(s80)
mean(s90)



#### SET SPECIFICITIES

results_all = c()
#### RESULTS WITH SET SENS AT 70,80,90 AND VICE VERSA
#########

results_all = c()
s70 <- c()
s80<- c()
s90 <- c()

for (each in 1:100) {
  trn_idx <- sample(1:N, 0.8 * N)
  train <- sf_new[trn_idx, ]
  test <- sf_new[-trn_idx, ]
  
  mod <- ranger(sfgr ~ ., data = train, num.trees = 1000, probability = TRUE)
  predictions <- predict(mod, data = test)
  risk_scores <- predictions$predictions[, 1]
  
  # Fit the ROC model
  roc_model <- roc(test$sfgr, risk_scores)
  
  yy <- data.frame(sens = roc_model$sensitivities, spec = roc_model$specificities) 
  
  target_specificities <- c(0.70, 0.80, 0.90)
  
  
  spec_indices <- sapply(target_specificities, function(target_spec) which.min(abs(yy$spec - target_spec)))
  
  # Extract the sensitivity and specificity values at the identified indices
  spec_values <- yy$spec[spec_indices]
  
  t7 <- spec_values[1]
  t8 <- spec_values[2]
  t9 <- spec_values[3]
  
  
  s7 <- yy %>% filter(yy$spec == t7)
  sens7 <- mean(s7$sens)
  s8 <- yy %>% filter(yy$spec == t8)
  sens8 <- mean(s8$sens)
  s9 <- yy %>% filter(yy$spec == t9)
  sens9 <- mean(s9$sens)
  
  s70 <- rbind(sens7, s70)
  s80 <- rbind(sens8, s80)
  s90 <- rbind(sens9, s90)
}
mean(s70)
mean(s80)
mean(s90)



##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS
##### JUST CLIMATE PREDICTORS

sf_c <- result %>% dplyr::select(c('sfgr', Minimum_NDWI_lag2:Maximum__1_km_16_days_NDVI_lag1))

N=dim(sf_c)[1]
AUC=c()
CI=c()
for (each in 1:100){
  print(each)
  trn_idx = sample(1:N,.8*N)
  train=sf_c[trn_idx,]
  test=sf_c[-trn_idx,]
  mod= ranger(sfgr~.,data=train, num.trees = 1000,probability=T)
  pred_train1=predict(mod,train)$predictions[,2][which(train$sfgr==1)]
  pred_train0=predict(mod,train)$predictions[,2][which(train$sfgr==0)]
  preds_test=predict(mod,test)$predictions[,2]
  AUC =rbind(AUC,c(as.numeric(roc(test$sfgr,preds_test)$auc)))
  CI =rbind(CI, c(as.numeric(ci(roc(test$sfgr,preds_test)))))
}



CImin=mean(CI[,1])
AUC=mean(CI[,2])
CImax=mean(CI[,3])




pred_fun = function(X.model, newdata) {
  predict(X.model, newdata)$predictions[,2]
}




CI_df <- data.frame(
  i = integer(0),
  CImin = double(0),
  AUC = double(0),
  CImax = double(0)
)

i_values <- c(20,15,10,9,8,7,6,5,4,3,2)

N=dim(sf_c)[1]
CI=c()
for(each in 1:100){
  trn_idx = sample(1:N, .8*N)
  train=sf_c[trn_idx,]
  test=sf_c[-trn_idx,]
  mod=ranger(sfgr~., data=train, num.trees = 1000, probability=T)
  VIP <- vip(object = mod, method = 'permute', target = 'sfgr', metric = 'roc_auc', train= train,
             pred_wrapper= pred_fun, event_level = 'second', num_features = 25)
  VIP <- VIP$data
  
  for(i in i_values) {
    selected_variables <- VIP %>%
      slice(1:i) %>%
      pull(Variable)
    
    train1 <- train[,c(selected_variables, 'sfgr')]
    test1 <- test[,c(selected_variables, 'sfgr')]
    mod=ranger(sfgr~., data=train1, num.trees = 1000, probability=T)
    pred_train1=predict(mod,train1)$predictions[,2][which(train1$sfgr==1)]
    pred_train0=predict(mod,train1)$predictions[,2][which(train1$sfgr==0)]
    preds_test=predict(mod,test1)$predictions[,2]
    
    CI= c(as.numeric(ci(roc(test$sfgr, preds_test, plot=F))))
    CI_df <- rbind(CI_df, data.frame(i = i, CImin = CI[1], AUC = CI[2], CImax = CI[3]))
  }
}


averages <- CI_df %>%
  group_by(i) %>%
  summarize(
    Avg_CImin = mean(CImin),
    Avg_CImax = mean(CImax),
    Avg_AUC = mean(AUC)
  )


if (sum(test$sfgr == 1) > 0) {
  # At least one instance with qf == 1 found in the test set
  break
}



# Initialize vectors to store probabilities and their corresponding classes
probabilities <- c()
classes <- c()

AUC=c()
CI=c()
N=dim(sf_new)[1]
for (each in 1:100){
  print(each)
  trn_idx = sample(1:N,.8*N)
  train=sf_new[trn_idx,]
  test=sf_new[-trn_idx,]
  mod= ranger(sfgr~.,data=train, num.trees = 1000,probability=T)
  pred_train1=predict(mod,train)$predictions[,2][which(train$sfgr==1)]
  pred_train0=predict(mod,train)$predictions[,2][which(train$sfgr==0)]
  preds_test=predict(mod,test)$predictions[,2]
  AUC =rbind(AUC,c(as.numeric(roc(test$sfgr,preds_test)$auc)))
  CI =rbind(CI, c(as.numeric(ci(roc(test$sfgr,preds_test)))))
  
  # Store probabilities and their corresponding classes
  probabilities <- c(probabilities, preds_test)
  classes <- c(classes, ifelse(test$sfgr == 1, "Positive", "Negative"))
}

CImin=mean(CI[,1])
AUC=mean(CI[,2])
CImax=mean(CI[,3])

# Create a dataframe for probabilities and their corresponding classes
df_prob <- data.frame(Probability = probabilities, Class = classes)

# Plot histogram
library(ggplot2)
p <- ggplot(df_prob, aes(x = Probability, fill = Class)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(title = "Histogram of Probabilities", x = "Probability", y = "Count") +
  scale_fill_manual(values = c("Positive" = "blue", "Negative" = "red")) +
  theme_minimal()



ggsave("hist_probs.pdf", plot = p, device = "pdf", width = 7, height = 7,  dpi = 600)















N = dim(sf_new)[1]
str(sf_new)

# Create a list of data frames to store CIs for each 'i'
results_list <- list()
i_values <- c( 20,15,10,9,8,7,6,5,4,3,2)

results_list <- list()

for (each in 1:100) {
  while (TRUE) {
    trn_idx <- sample(1:N, 0.8 * N)
    train <- sf_new[trn_idx,]
    test <- sf_new[-trn_idx,]
    
    if (sum(test$sfgr == 1) > 0) {
      # At least one instance with qf == 1 found in the test set
      break
    }
  }
  
  rf_model_original <- randomForest(sfgr ~ ., data = train, ntree = 1000, keep.forest = TRUE, keep.inbag = TRUE)
  rf_model_new <- randomForest(sfgr ~ ., data = train, ntree = 1000, keep.forest = TRUE, keep.inbag = TRUE)
  
  
  # Your code that calls permimp
  perm_imp <- permimp(rf_model_new, forest = rf_model_original$forest, inbag = rf_model_original$inbag, 
                      data = train, do_check = FALSE)
  
  pid <- data.frame(perm_imp$values)
  names <- rownames(pid)
  pidf <- data.frame(cbind(names, pid))
  pidf <- pidf %>% arrange(desc(perm_imp.values))
  
  for (i in i_values) {
    # Extract the desired features based on 'i'
    pidf_i <- pidf[1:i, ]
    train1_i <- train[, c(pidf_i$names, 'sfgr')]
    test1_i <- test[, c(pidf_i$names, 'sfgr')]
    
    # Fit the model and calculate CI for this 'i'
    mod <- ranger(sfgr ~ ., data = train1_i, num.trees = 1000, probability = TRUE)
    pred_train1 <- predict(mod, train1_i)$predictions[, 2][which(train1_i$sfgr == 1)]
    preds_test_i <- predict(mod, test1_i)$predictions[, 2]
    CI_i <- data.frame(as.numeric(ci(roc(test1_i$sfgr, preds_test_i, plot = FALSE))))
    
    # Store the CI data frame for this 'i' in the list
    if (is.null(results_list[[as.character(i)]])) {
      results_list[[as.character(i)]] <- list()
    }
    results_list[[as.character(i)]][[each]] <- CI_i
  }
}




results_list[[1]]

dd <- data.frame(unlist(results_list[1]))

colnames(dd)[1] <- "Value"

average1 <- mean(dd$Value[seq(1, nrow(dd), by = 3)])
average2 <- mean(dd$Value[seq(2, nrow(dd), by = 3)])
average3 <- mean(dd$Value[seq(3, nrow(dd), by = 3)])


sfgr_performance <- data.frame(model = '20',
                               CImin = average1,
                               AUC = average2,
                               CImax=average3)

new_row <- c(2, average1, average2, average3)
sfgr_performance <- rbind(sfgr_performance, new_row)


sfgr_performance[, 2:4] <- round(sfgr_performance[, 2:4], 3)

sfgr_performance <- sfgr_performance[-c(6:10),]

done <- sfgr_performance %>% arrange(desc(AUC))



sfgr_performance$combined <- paste(sfgr_performance$CImin, sfgr_performance$CImax, sep = '-')

# Remove the second and fourth columns
sfgr_performance <- sfgr_performance[-c(2, 4)]


rf_model_new <- randomForest(sfgr ~ ., data = sf_new, ntree = 1000, keep.forest = TRUE, keep.inbag = TRUE)
perm_imp <- permimp(rf_model_new, conditional=T,
                    data = sf_new,AUC=T)


pid <- data.frame(perm_imp$values)
names <- rownames(pid)
pidf <- data.frame(cbind(names, pid))
pidf <- pidf %>% arrange(desc(perm_imp.values))





