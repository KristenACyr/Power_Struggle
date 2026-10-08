######## Reference Tag Paper - Position Yield - daily mean counts 

## load libraries 
library(tidyverse)
library(ggplot2)
library(geodata)
library(terra)
library(dplyr)
library(sf)
library(lubridate)
library(ggmap)
library(stringr)
library(boot)
library(stats)
library(multcomp)


## set working directory 
setwd("D:/Documents/Masters/CoLab paper/REFTAG/R Code")

## open datasheet and set name to "highlow" 
highlow = read.csv("ReferenceTag_Positions.csv")

## change date/time column from character to date/time format
highlow$detection_timestamp_AST = as.POSIXct(highlow$detection_timestamp_AST, 
                                             format = "%Y-%m-%d %H:%M")



## calculate the total number of positions per tag power (high vs low) and receivers used (3, 4, 5)
tot_positions = highlow %>% 
  group_by(power, configuration, nRxUsed) %>% 
  summarise(n = n(), .groups = "drop")

## 'n' is produced as an integer, so change from integer to numeric 
tot_positions$n <- as.numeric(tot_positions$n)



## determine the number of days each configuration was deployed using ORDINAL DAY rather than date/time
deployed_days <- highlow %>%
  group_by(configuration, year) %>%
  summarise(start_day = min(ordinal_day, na.rm = TRUE),
    end_day = max(ordinal_day, na.rm = TRUE),
    n_days = end_day - start_day + 1)

## determine the number of detections/day for each configuration/power level/number of receivers
days_detections <- deployed_days %>% 
  left_join(tot_positions, by = "configuration") %>% 
  mutate(detections_per_day = n / n_days)

## convert configuration from numeric to factor
days_detections$configuration <- as.factor(days_detections$configuration)


########## GENERAL LINEAR MODELS TO DETERMINE WHAT AFFECTS DETECTIONS PER DAY #############
## (models estimate how the # of detections per day (n) changes with the tag powers (tag_power), 
#the number of receivers used (nRxUsed), and the configuration (configuration))



## GLM1 - glm with POWER x NUMBER OF RECEIVERS USED and POWER x CONFIGURATION 
highlow_glm1 <- glm(detections_per_day~ (power * nRxUsed) + (power * configuration),
                    data = days_detections, family = gaussian(link = "log"))


summary(highlow_glm1)


## GLM2 - glm with POWER x CONFIGURATION and an additional covariate of NUMBER OF RECEIVERS USED
highlow_glm2 <- glm(detections_per_day~ (power * configuration) + nRxUsed, 
                    data = days_detections, family = gaussian(link = "log"))


summary(highlow_glm2)

## GLM3 - glm with POWER x NUMBER OF RECEIVERS USED and an additional covariate of CONFIGURATION 
highlow_glm3 <- glm(detections_per_day~ (power * nRxUsed) + configuration, 
                    data = days_detections, family = gaussian(link = "log"))


summary(highlow_glm3)


#glm4 - no interaction effects 
highlow_glm4 <- glm(detections_per_day~ power + nRxUsed + configuration, 
                    data = days_detections, family = gaussian(link = "log"))


summary(highlow_glm4)

## checking residuals 
qqnorm(residuals(highlow_glm4))
hist(residuals(highlow_glm4))




################### AIC TO DETERMINE WHICH MODEL BEST FITS THE DATA ################
## load libraries
library(MuMIn)
library(AICcmodavg)

## list the models 
model <- list(highlow_glm1, highlow_glm2, highlow_glm3, highlow_glm4)

## name the models 
model.names <- c('highlow_glm1', 'highlow_glm2', 'highlow_glm3', 'highlow_glm4')


## run the AIC 
aictab(cand.set = model, modnames = model.names)
# higlow_glm5 is the best fit model for this data as it has a Delta_AICc of 0.00





####################TUKEY POST HOC#############################################
#tukey post hoc used after running a glm 

# load library
library(emmeans)

# 1. Calculate estimated marginal means
emm_obj <- emmeans(highlow_glm4, ~ configuration)

# 2. Perform pairwise comparisons (Tukey is the default)
pairs(emm_obj)






