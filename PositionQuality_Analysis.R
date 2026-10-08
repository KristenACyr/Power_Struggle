#### Reference Tag Paper - MCP Analysis, Accuracy and Precision

## running packages
library(sf)
library(dplyr)
library(emmeans)
library(tidyverse)
library(MuMIn)
library(AICcmodavg)
library(writexl)
library(adehabitatHR)
library(ggplot2)
library(ggspatial)
library(scales)
library(sp)

#### open datasets
hltags = read.csv("ReferenceTag_Positions.csv")
deploy_points = read_csv("deploy_points.csv")

## keep 4 first columns
deploy_points = deploy_points[,c(1:4)]

# Convert deploy_points dataframe to spatial file using lat/lon columns
deploy_points_sf = st_as_sf(deploy_points, 
                            coords = c("deploy_lon", "deploy_lat"), 
                            crs = 4326,   # WGS84 (lat/lon)
                            remove = FALSE)  # keeps original columns

# Transform to UTM Zone 20N (meters)
deploy_points = st_transform(deploy_points_sf, crs = 32620)

# Extract lat and lon, to keep as coordinate columns, otherwise removes lat lon
coords <- st_coordinates(deploy_points)

#calculate lat lon as meters x and meters y
deploy_points$deploy_meters_x <- coords[,1]
deploy_points$deploy_meters_y <- coords[,2]

# Need to convert spatial file back to dataframe
deploy_points <- as.data.frame(deploy_points)

################################################################################
#turn varaibles into factors
hltags$nRxUsed <- as.factor(hltags$nRxUsed)
hltags$configuration <- as.factor(hltags$configuration)
hltags$power = as.factor(hltags$power)
deploy_points$configuration = as.factor(deploy_points$configuration)
deploy_points$power = as.factor(deploy_points$power)


######################CONVERT METERS COLUMNS NOT DEGREES!!!!############################

################reference tag dataframe
# Convert hltags dataframe to spatial file using lat/lon columns
hltags_sf = st_as_sf(hltags, 
                     coords = c("Longitude", "Latitude"), 
                     crs = 4326,   # WGS84 (lat/lon)
                     remove = FALSE)  # keeps original columns

# Transform to UTM Zone 20N (meters)
hltags = st_transform(hltags_sf, crs = 32620)

# Extract lat and lon, to keep as coordinate columns, otherwise removes lat lon
coords <- st_coordinates(hltags)

#calculate lat lon as meters x and meters y
hltags$meters_x <- coords[,1]
hltags$meters_y <- coords[,2]

# Need to convert spatial file back to dataframe
hltags <- as.data.frame(hltags)




#######################  ACCURACY CALCULATIONS #########################################
# Calculate distance between estimated and actual deployment locations (I.E., DEPLOY_POINTS)
hltags <- hltags %>%
  left_join(
    deploy_points %>%
      dplyr::select(configuration, deploy_meters_x, deploy_meters_y, power),
    by = c("configuration", "power")
  ) %>%
  mutate(
    accuracy_m = sqrt(
      (meters_x - deploy_meters_x)^2 +
        (meters_y - deploy_meters_y)^2
    )
  )

######################################################################################


########################### PRECISION CALCULATIONS ###################################
# Calculate mean position for each configuration and power
mean_positions <- hltags %>%
  group_by(configuration, power) %>%
  summarise(
    mean_x = mean(meters_x, na.rm = TRUE),
    mean_y = mean(meters_y, na.rm = TRUE),
    .groups = "drop"
  )

# Add mean positions to each observation
hltags_precision <- hltags %>%
  left_join(mean_positions, by = c("configuration", "power"))

# Calculate distance from mean position (precision) to the actual position
hltags <- hltags_precision %>%
  mutate(
    precision_m = sqrt(
      (meters_x - mean_x)^2 +
        (meters_y - mean_y)^2
    )
  )
############################################################################################

















#########################################################################################
###########################Linear Regress models for precision###########################

# Log and check distribution
hist(log10(hltags$precision_m))


## log and check distribution
hist(log10(hltags$accuracy_m))

#create logged column so we can compare later with true value and make sure it's no below 0
hltags$logged_precision = (log10(hltags$precision_m))

################################## ACCURACY MODELS #########################################
#### shifting accuracy 
summary(hltags$accuracy_m)
hltags$accuracy_adj <- hltags$accuracy_m + 0.01


#glm1 - glm with power interacting with number of receivers used and power interacting with configuration 
#**power because we are comparing high vs low 
#**added interactions between main effects and number of receivers and main effects and positions (configs)
##**Model estimates how the number of detections (n) changes with the tag power 
##levels (tag_power), the number of receivers used (RxUsed) and the configuration 
##(position)


highlow_glm1acc <- glm(accuracy_adj
                       ~ (power * nRxUsed) 
                       + (power * configuration),
                       family = Gamma(link = "log"),
                       data = hltags)

summary(highlow_glm1acc)

# CHECKING RESIDUALS #
qqnorm(residuals(highlow_glm1acc))
hist(residuals(highlow_glm1acc))

# model coefficients and p-values
summary(highlow_glm1acc)$coefficients
# Calculate confidence intervals
confint(highlow_glm1acc)


#glm2 - glm with power interacting with configuration and number of receivers used

highlow_glm2acc <- glm(accuracy_adj ~ 
                      (power * configuration) 
                       + nRxUsed,
                       family = Gamma(link = "log"),
                       data = hltags)


summary(highlow_glm2acc)

##checking residuals
qqnorm(residuals(highlow_glm2acc))
hist(residuals(highlow_glm2acc))


#glm3 - power interacting with number of recievers used and configuration
highlow_glm3acc <- glm(accuracy_adj~ (power * nRxUsed) + configuration, 
                       family = Gamma(link = "log"),
                       data = hltags)

summary(highlow_glm3acc)

### Residuals
qqnorm(residuals(highlow_glm3acc))
hist(residuals(highlow_glm3acc))


#glm4 - no interaction effects 
highlow_glm4acc <- glm(accuracy_adj~ 
                      power 
                    + nRxUsed 
                    + configuration, 
                    family = Gamma(link = "log"),
                    data = hltags)

summary(highlow_glm4acc)

####################################################################
##################AIC TO DETERMINE WHICH MODEL BEST FITS DATA################


#Ran diff GLMs and AIC showed which was the best for this data frame
#list models 
model <- list(highlow_glm1acc, highlow_glm2acc, highlow_glm3acc, highlow_glm4acc)

#add column model names
model.names <- c('highlow_glm1', 'highlow_glm2', 'highlow_glm3', 'highlow_glm4')

## run AIC
aictab(cand.set = model, modnames = model.names)

#######################################################################################

######################### Tukey post-hoc test for GLM 1 ###############################

# Set maximum number of printed results
options(max.print = 10000)

# Calculate emmeans for nRxUsed × power
emmnRxusedacc <- emmeans(
  highlow_glm1acc,
  ~ power * nRxUsed
)
# post-hoc comparisons
pairs(emmnRxusedacc)

# Calculate emmeans for power × configuration
emmpoweracc <- emmeans(
  highlow_glm1acc,
  ~ power * configuration
)
# Pairwise post-hoc comparisons
pairs(emmpoweracc)

###############################################################################



#####################################Linear Regress models - PRECISION  ####################################
#glm1 - glm with power interacting with number of receivers used and power interacting with configuration
highlow_glm1prec <- glm(logged_precision 
                        ~ (power * nRxUsed) 
                        + (power * configuration), # (nRxUsed * configuration), 
                        data = hltags)

summary(highlow_glm1prec)$coefficients

## Residuals ##
qqnorm(residuals(highlow_glm1prec))
hist(residuals(highlow_glm1prec))


#glm2 - glm with power interacting with configuration and number of receivers used 
highlow_glm2prec <- glm(logged_precision 
                        ~ (power * configuration) 
                        + nRxUsed, 
                        data = hltags)

summary(highlow_glm2prec)

##checking residuals
qqnorm(residuals(highlow_glm2prec))
hist(residuals(highlow_glm2prec))




#glm3 - power interacting with number of recievers used and configuration
highlow_glm3prec <- glm(logged_precision
                        ~ (power * nRxUsed) 
                        + configuration, 
                        data = hltags)

summary(highlow_glm3prec)

## Residuals
qqnorm(residuals(highlow_glm3prec))
hist(residuals(highlow_glm3prec))


#glm4 - no interaction effects 
highlow_glm4prec <- glm(logged_precision
                        ~ power 
                        + nRxUsed 
                        + configuration, 
                        data = hltags)

qqnorm(residuals(highlow_glm4prec))
hist(residuals(highlow_glm4prec))
summary(highlow_glm4prec)


####################################################################
##################AIC TO DETERMINE WHICH MODEL BEST FITS DATA################

## list models
model <- list(highlow_glm1prec, highlow_glm2prec, highlow_glm3prec, highlow_glm4prec)


## name columns
model.names <- c('highlow_glm1', 'highlow_glm2', 'highlow_glm3', 'highlow_glm4')

## run aic
aictab(cand.set = model, modnames = model.names)


#######################################################################################
######################### Tukey post-hoc test for GLM 1 

# Calculate emmeans for nRxUsed × configuration
emmnRxusedprec <- emmeans(
  highlow_glm1prec,
  ~ power * nRxUsed
)

# post-hoc comparisons
pairs(emmnRxusedprec)


# Calculate emmeans for power × configuration
emmpowerprec <- emmeans(
  highlow_glm1prec,
  ~ power * configuration
)

# Pairwise post-hoc comparisons
pairs(emmpowerprec)


################### Finding the means of Accuracy and Precision ############

# Select accuracy, precision, power, and configuration columns
exportaccprec <- hltags[, c("accuracy_m", "precision_m", "power", "configuration")]

# Calculate mean accuracy and precision by power AND configuration
meansaccprec <- exportaccprec %>%
  group_by(power, configuration) %>%
  summarise(
    mean_accuracy = mean(accuracy_m, na.rm = TRUE),
    mean_precision = mean(precision_m, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

####### Exporting
# Write results to an Excel file
write_xlsx(meansaccprec, "accuracy_precision_export.xlsx")





#####################################################################################
##################################### MCP CODE ######################################

################## Housekeeping prior to making MCP - datasets ######################

# Create dataframe for each configuration
highlowconfig1 <- hltags %>% filter(configuration == 1)
highlowconfig2 <- hltags %>% filter(configuration == 2)
highlowconfig3 <- hltags %>% filter(configuration == 3)
highlowconfig4 <- hltags %>% filter(configuration == 4)
highlowconfig5 <- hltags %>% filter(configuration == 5)
highlowconfig6 <- hltags %>% filter(configuration == 6)
highlowconfig7 <- hltags %>% filter(configuration == 7)

################### Making MCP ######################


### NOTE: For each configuration - change out the number with the associated config#
### Ex. If you want to create MCP for Config 2 - all code below should read ,"highlowconfig2" etc. 

# Convert configuration 1 to spatial format
if (inherits(highlowconfig1, "sf")) {
  highlowconfig1_sp <- as(highlowconfig1, "Spatial")
} else {
  highlowconfig1_sp <- highlowconfig1
  coordinates(highlowconfig1_sp) <- ~Longitude + Latitude
  proj4string(highlowconfig1_sp) <- CRS("+proj=longlat +datum=WGS84")
}

# Calculate 95% MCP
mcp_res1 <- adehabitatHR::mcp(
  highlowconfig1_sp[, "Transmitter"],
  percent = 95
)

# Convert MCP to sf
mcp_sf1 <- st_as_sf(mcp_res1)

# Label power groups
mcp_sf1$power <- factor(
  mcp_sf1$id,
  levels = c(1000, 2000),
  labels = c("Low power", "High power")
)

# Set UTM CRS
utm_crs <- 32620

# Transform MCP to UTM
mcp_sf1 <- st_transform(mcp_sf1, utm_crs)

# Find MCP centroids
centroids1 <- st_centroid(mcp_sf1)

# Add power labels
centroids1$power <- mcp_sf1$power


#####################################################
### Reference tags for configuration 1

# Convert high ref tag to spatial data
high_sf1 <- st_as_sf(
  deploy_points[
    deploy_points$configuration == 1 &
      deploy_points$power == "high",
  ],
  coords = c("deploy_lon", "deploy_lat"),
  crs = 4326
) |>
  st_transform(utm_crs)

# Convert low ref tag to spatial data
low_sf1 <- st_as_sf(
  deploy_points[
    deploy_points$configuration == 1 &
      deploy_points$power == "low",
  ],
  coords = c("deploy_lon", "deploy_lat"),
  crs = 4326
) |>
  st_transform(utm_crs)


#####################################################
### Get plot boundaries

bbox <- sf::st_bbox(mcp_sf1)

x_limits <- c(
  bbox[["xmin"]],
  bbox[["xmax"]]
)

y_limits <- c(
  bbox[["ymin"]],
  bbox[["ymax"]]
)


#####################################################
### Create plot

ggplot() +
  
  # Plot MCP areas
  geom_sf(
    data = mcp_sf1,
    aes(fill = power),
    alpha = 0.3,
    colour = "black",
    linewidth = 1.2
  ) +
  
  # Plot MCP centroids
  geom_sf(
    data = centroids1,
    aes(colour = power),
    shape = 23,
    size = 3,
    show.legend = FALSE
  ) +
  
  # Plot high reference tag
  geom_sf(
    data = high_sf1,
    aes(colour = "High Reference Tag"),
    size = 5
  ) +
  
  # Plot low reference tag
  geom_sf(
    data = low_sf1,
    aes(colour = "Low Reference Tag"),
    size = 5
  ) +
  
  # Set axis breaks
  scale_x_continuous(
    breaks = x_limits,
    labels = scales::label_number(big.mark = "")
  ) +
  
  scale_y_continuous(
    breaks = y_limits,
    labels = scales::label_number(big.mark = "")
  ) +
  
  # Set plot theme
  theme_classic() +
  theme(
    plot.margin = margin(
      t = 20, r = 20, b = 10, l = 25
    ),
    
    axis.line = element_line(
      linewidth = 1.2,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      size = 16,
      face = "bold",
      margin = margin(r = 15)
    ),
    
    axis.text.x = element_text(
      size = 28,
      face = "bold",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    axis.text.y = element_text(
      size = 28,
      face = "bold"
    )
  ) +
  
  # Set MCP colours
  scale_fill_manual(
    name = "Tag power",
    values = c(
      "Low power" = "#00BFC4",
      "High power" = "#F8766D"
    )
  ) +
  
  # Set reference tag colours
  scale_color_manual(
    name = NULL,
    values = c(
      "High Reference Tag" = "#F8766D",
      "Low Reference Tag" = "#00BFC4",
      "Low power" = "#00BFC4",
      "High power" = "#F8766D"
    )
  ) +
  
  # Set map boundaries and CRS
  coord_sf(
    xlim = x_limits,
    ylim = y_limits,
    expand = FALSE,
    datum = utm_crs
  )

#######################################################





