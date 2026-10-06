##required libraries
library(lubridate)
library(tidyverse)
library(dplyr)
library(tidyr)
library(cowplot)
library(lme4)
library(AICcmodavg)
library(ggeffects)
library(emmeans)
############################################################################


###########################IMPORTED DATASETS####################################
#############################################################################
refdetects = read.csv("refdetects.csv")

hltags = read.csv("Highlow_AllYearsNEW.csv",
                  colClasses = c(
                    Latitude = "numeric",
                    Longitude = "numeric"
                  ))

###############################################################################















#############################FORMAT DATASETS FOR ANALYSIS#####################
################################################################################
options(digits = 15)
hltags$Longitude[1]
hltags$Latitude[1]

#Change date/time column to date/time format
hltags$detection_timestamp_AST = 
  as.POSIXct(hltags$detection_timestamp_AST, format = "%Y-%m-%d %H:%M:%S")


hltags = hltags %>% 
  mutate(
    detection_timestamp_AST = as.POSIXct(detection_timestamp_AST),
    ordinal_day = lubridate::yday(detection_timestamp_AST)   # day of year
  )  



refdetects <- refdetects %>%
  mutate(
    detection_timestamp_AST = if_else(
      str_detect(as.character(detection_timestamp_AST), ":"),
      as.character(detection_timestamp_AST),                     # already has time
      paste0(as.character(detection_timestamp_AST), " 00:00:00") # add midnight
    ),
    detection_timestamp_AST = lubridate::ymd_hms(
      detection_timestamp_AST,
      tz = "America/Halifax"
    )
  )


#Distinguish between detections versus positions
refdetects$detection_type = "raw"
hltags$detection_type = "VPS"
####################################################################################






#############CODE FOR POSITION PROBABILITY FIGURES AND ANALYSIS###################


#SUMMARY INFORMATION


#calculate total detections from each receiver for each year
totdets_year_rec = refdetects %>% 
  group_by(receiver_id, year) %>% 
  summarise(total_dets = n(), .groups = "drop")


#calculate total positions that each receiver was apart of estimating across each year
totpos_year_rec = hltags %>%
  mutate(RxUsed = str_split(RxUsed, " ")) %>%   # split by space
  unnest(RxUsed) %>%
  #distinguish which receivers were associated with each station as provided
  #by Innvosea
  mutate(receiver_id = case_when(
    year == 2024 & RxUsed == "St01" ~ 710470,
    year == 2024 & RxUsed == "St02" ~ 710471,
    year == 2024 & RxUsed == "St03" ~ 710472,
    year == 2024 & RxUsed == "St04" ~ 710473,
    year == 2025 & RxUsed == "St01" ~ 710470,
    year == 2025 & RxUsed == "St02" ~ 710471,
    year == 2025 & RxUsed == "St03" ~ 710472,
    year == 2025 & RxUsed == "St04" ~ 710473,
    year == 2025 & RxUsed == "St05" ~ 710534,
    TRUE ~ NA_integer_
  )) %>%
  group_by(receiver_id, year) %>%
  summarise(total_pots = n(), .groups = "drop")






##Summary table of total positions and detections

#calculate total positions and detections for each reference tag (high versus low)
#across all configurations and years
tagPorp <- hltags %>%
  group_by(power, year, configuration) %>%
  summarise(total_positions = n(), .groups = "drop") %>%
  left_join(
    refdetects %>%
      group_by(power, year, configuration) %>%
      summarise(total_detections = n(), .groups = "drop"),
    by = c("power", "year", "configuration")
  ) %>%
  mutate(percent_positions_per_detection = 
           round((total_positions / total_detections)*100, 3))











#Daily probability of positions estimated

#calculate the probability high versus low power tags were positioned from 
#detections per day
tag_percent <- refdetects %>%
  group_by(power, ordinal_day, configuration) %>%
  summarise(total_detections = n(), .groups = "drop") %>%
  left_join(
    hltags %>%
      group_by(power, ordinal_day, configuration) %>%
      summarise(total_positions = n(), .groups = "drop"),
    by = c("power", "ordinal_day", "configuration")
  ) %>%
  mutate(
    total_positions = tidyr::replace_na(total_positions, 0),
    percent_positions_per_detection =
      round((total_positions / total_detections) * 100, 3)
  )





#Figure 5 generated from daily probability of position estimates
p = ggplot(tag_percent,
       aes(x = as.factor(configuration),
           y = percent_positions_per_detection,
           colour = power)) +
  geom_point(alpha = 0.7,size = 5, shape = 3,
             position = position_jitter(width = 0.1, height = 0)) +
  stat_summary(
    aes(group = power),
    fun = mean,
    geom = "line",
    linewidth = 2
  ) +
  stat_summary(
    aes(group = power),
    fun = mean,
    geom = "point",
    size = 4
  ) +
  scale_y_continuous(
    limits = c(0, 30)) +
  labs(
    x = "Configuration",
    y = "Daily Detection Porportion (%)",
    colour = "Tag power"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30))
p
#################################################################################




























##########################GLM MODEL FOR POSITION PROBABILITY#################
#############################################################################

#Create dataframe to compare position probability for glm model
tag_percentMODELB <- hltags %>%
  #using my positions group by tag power, day, and ref tag config
  group_by(power, ordinal_day, configuration) %>%
  #using this grouping get their total positions
  summarise(total_positions = n(), .groups = "drop") %>%
  #join this with my detection dataframe grouping and summarizing the same way
  #but instead of positions get detections
  left_join(
    refdetects %>%
      group_by(power, ordinal_day, configuration) %>%
      summarise(total_detections = n(), .groups = "drop"),
    by = c("power", "ordinal_day", "configuration")
  ) %>%
  # calculate the porportion of detections that made it into positions
  #and calculate the total times a detection was NOT positioned
  mutate(
    prop_positions_per_detection = total_positions / total_detections,
    failures = total_detections - total_positions 
  )


#need to make tag power and config factor for the model
tag_percentMODELB$power = as.factor(tag_percentMODELB$power)
tag_percentMODELB$configuration = as.factor(tag_percentMODELB$configuration)



#test 2 models for with and without interaction of power and configuration


#NOTE: CANNOT ADD receiver ID as factor, because there are multiple receivers 
#used for positions. 

# does probability of successful position depend on 
#interaction between power and configuration
modelA <- glm(
  cbind(total_positions, failures) ~ power * configuration,
  family = binomial,
  data = tag_percentMODELB
)

summary(modelA)


hist(residuals(modelA))
qqnorm(residuals(modelA))
qqline(residuals(modelA))

# does probability of successful position depend on 
#power and configuration
modelB <- glm(
  cbind(total_positions, failures) ~ power + configuration,
  family = binomial,
  data = tag_percentMODELB
)

hist(residuals(modelB))
qqnorm(residuals(modelB))
qqline(residuals(modelB))


#Run AIC for both models
models <- list(
  GLM1 = modelA,
  GLM2 = modelB
)

aictab(models)

#Model with interaction best model
summary(modelA)


#Plot model predictions
pred <- ggpredict(modelA, terms = c("power", "configuration"))
plot(pred)

# Estimated marginal means for power within each configuration
emm <- emmeans(modelA, ~ power | configuration)

# Pairwise comparisons (High vs Low) within each configuration
power_comparisons <- pairs(emm, adjust = "tukey")
power_comparisons

####################################################################################





