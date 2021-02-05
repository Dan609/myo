# Cell movement trajectory analysis by Dan Bobkov, 2019 # dan.bobkov@gmail.com
# Script desined to analyze tracks collected in ImageJ Manual tracking
# plugin from CQ1 cytometer derived images, recorded during 24 h experiment
# with MSCWJ-1 cells conditioned by long-term cultivation.
# Images were collected from passages 9, 15, 28, 36 to model the classical
# replicative senescense effects.
# If you use this script in your publications, please cite:
# Bobkov D, Polyanskaya A, Musorina A, Lomert E, Shabelnikov S, Poljanskaya G. 
# Replicative senescence in MSCWJ-1 human umbilical cord mesenchymal stem cells is marked by characteristic changes in motility, 
# cytoskeletal organization, and RhoA localization. Mol Biol Rep. 2020; 47:3867–3883. https://doi.org/10.1007/s11033-020-05476-6
# This script based on trajr package
# If you use trajr in your publications, please cite McLean DJ, Skowron Volponi MA. trajr:
# An R package for characterisation of animal trajectories. Ethology. 2018; 12739.
# https://doi.org/10.1111/eth.12739.
#
#################################################
################## Set workspace  ###############
#################################################
# Load libraries-------------
library(PerformanceAnalytics)
library(GGally)
library(trajr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(plyr)
library(xlsx)
library(car)
library(rgl)
library(dunn.test)
library(DescTools)
library(FactoMineR)
library(factoextra)
library(Hmisc)
library(gplots)
library(ggsignif)
library(stringi)
library(PMCMRplus)
library(pca3d)
library(xtable)
library(DAAG)
library(cluster)

#################################################
################## Initialization ###############
#################################################
# Set names---------------
alltracks <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                      c("track",
                         "length",
                         "distance",
                         "straight",
                         "square_displacement",
                         "mean_speed",
                         "sd_speed",
                         "max_speed",
                         "min_speed",
                         "sinuosity",
                         "emax",
                         "DC",
                         "SDDC",
                         "mean_angle",
                         "probe"))
tracks_p09 <- alltracks
tracks_p15 <- alltracks
tracks_p28 <- alltracks
tracks_p36 <- alltracks
# Remove outliers function------
outliers.rm <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75))
  q <- unname(q)
  for (i in x)
  {
    if ( i < q[1] - 1.5*IQR(x) ||
         i > q[2] + 1.5*IQR(x)) {
      x <- x[x != i]
      return(x)
    }
    else {
      return(x)
    }
  } }
#################################################
##### Trajectory Analisys Functions #############
#################################################
# Load Track analysis functions
# p09-------------
traj_analysis_p09 <- function(input) {

  data <- read.csv(input)

  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))

  for (i in unique(data$Track)) {

    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")

    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 95/3600/24)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 1.315789 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 1.3158 / 1, "micrometer") # 0.6579 for p28
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)

    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj), # D/L ratio
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           probe = 'p09')

    head(traj_params)

  }

  # print(traj_params)
  write.csv(traj_params, file = 'traj_p09.csv')
  tracks <<- traj_params
  return(traj_params) } # 96 frames
# p15-----------
traj_analysis_p15 <- function(input) {

  data <- read.csv(input)

  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))

  for (i in unique(data$Track)) {

    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")

    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 95/3600/24)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 1.315789 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 1.3158 / 1, "micrometer")
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)

    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj),
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           probe = 'p15')

    head(traj_params)

  }

  # print(traj_params)
  write.csv(traj_params, file = 'traj_p15.csv')
  tracks <<- traj_params
  return(traj_params)} # 96 frames
# p28-----------
traj_analysis_p28 <- function(input) {

  data <- read.csv(input)

  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))

  for (i in unique(data$Track)) {

    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")

    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 95/3600/24)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.6579 / 1, "micrometer")
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)

    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj),
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           probe = 'p28')

    head(traj_params)

  }

  # print(traj_params)
  write.csv(traj_params, file = 'traj_p28.csv')
  tracks <<- traj_params
  return(traj_params)} 

# p28-----------
traj_analysis_p36 <- function(input) {

  data <- read.csv(input)

  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))

  for (i in unique(data$Track)) {

    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")

    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 95/3600/24)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.6579 / 1, "micrometer")
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)

    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj),
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           probe = 'p28')

    head(traj_params)

  }

  # print(traj_params)
  write.csv(traj_params, file = 'traj_p28.csv')
  tracks <<- traj_params
  return(traj_params)}
#################################################
##### Tracks importing from local dirs ##########
#################################################
# Choose dir ---------p09
file_list_p09 <- list.files(path = , choose.dir(default = "",
                                            caption = "Select folder"),
                        pattern = "csv",
                        all.files = FALSE,
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = FALSE, include.dirs = FALSE,
                        no.. = FALSE)
# Choose dir ---------p15
file_list_p15 <- list.files(path = , choose.dir(default = "",
                                            caption = "Select folder"),
                        pattern = "csv",
                        all.files = FALSE,
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = FALSE, include.dirs = FALSE,
                        no.. = FALSE)
# Choose dir ---------p28
file_list_p28 <- list.files(path = , choose.dir(default = "",
                        caption = "Select folder"),
                        pattern = "csv",
                        all.files = FALSE,
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = FALSE, include.dirs = FALSE,
                        no.. = FALSE)
# Choose dir ---------p36
file_list_p36 <- list.files(path = , choose.dir(default = "",
                                            caption = "Select folder"),
                        pattern = "csv",
                        all.files = FALSE,
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = FALSE, include.dirs = FALSE,
                        no.. = FALSE)
#################################################
############# Call to function ##################
#################################################

# Start scan for p09
for (file_name in file_list_p09) {
  traj_analysis_p09(file_name)
  tracks_p09 <- rbind(tracks_p09, tracks)
}
# Start scan for p15
for (file_name in file_list_p15) {
  traj_analysis_p15(file_name)
  tracks_p15 <- rbind(tracks_p15, tracks)
}
# Start scan for p28
for (file_name in file_list_p28) {
  traj_analysis_p28(file_name)
  tracks_p28 <- rbind(tracks_p28, tracks)
}
# Start scan for p36
for (file_name in file_list_p36) {
  traj_analysis_p36(file_name)
  tracks_p36 <- rbind(tracks_p36, tracks)
}

#################################################
############# Prepare data ######################
#################################################
# Merge all tracks-----------------------
alltracks <- rbind(tracks_p09, tracks_p15, tracks_p28, tracks_p36) # collect all tracks
summary(alltracks)
write.csv(alltracks, file = 'alltracks.csv') #save results
# Order probe levels--------------------
alltracks$probe <- as.factor(alltracks$probe)
alltracks$probe <- ordered(alltracks$probe,
                      levels = c("p09", "p15", "p28", "p36"))
# Set time to hours, remove tracks and mean_angle columns----------------
all.h <- alltracks
head(all.h)
all.h <- cbind(all.h[,c(6,7,8,9)]*3600, all.h[,c(2,3,4,5,10,11,12,13,15)])
data <- all.h
head(data)

#################################################
######## Plot all correlations ##################
#################################################

# Plot all-in-one-------------------------
# plot(data)
# png()
# ggpairs(data)
chart.Correlation(data[,1:12], histogram=TRUE, pch=19)
ggcorr(data, palette = "RdBu", label = TRUE)
# dev.off()

#################################################
################ Plot 3D images #################
#################################################

## Plot data in three factor space
scatter3d(x = data$length, y = data$sinuosity, z = data$straight, 
          groups = data$probe,
          xlab=deparse(substitute(Length)), 
          ylab=deparse(substitute(Sinuosity)),
          zlab=deparse(substitute(Straightness)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE, 
          ellipsoid = TRUE,
          level=0.8, ellipsoid.alpha=0.5, id=FALSE,
          model.summary=FALSE)

#par3d("windowRect")
#r3dDefaults$windowRect <- c(0,100, 1600, 1600) 
#plot3d(data, col="black", size = 10)  
## Plot data in three factor space
scatter3d(x = data$SDDC, y = data$sinuosity, z = data$mean_speed, 
          groups = data$probe,
          xlab=deparse(substitute(SDDC)), 
          ylab=deparse(substitute(Sinuosity)),
          zlab=deparse(substitute(mean_speed)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE, 
          ellipsoid = TRUE,
          level=0.8, ellipsoid.alpha=0.5, id=FALSE,
          model.summary=FALSE)

## Plot data in three factor space
scatter3d(x = data$length, y = data$distance, z = data$DC, 
          groups = data1$probe,
          xlab=deparse(substitute(length)), 
          ylab=deparse(substitute(distance)),
          zlab=deparse(substitute(DC)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE, 
          ellipsoid = TRUE,
          level=0.8, ellipsoid.alpha=0.5, id=FALSE,
          model.summary=FALSE)

## Plot data in 3D factor space-----------------------------
scatter3d(x = data$mean_speed, y = data$sinuosity, z = data$straight, 
          groups = data$probe,
          xlab=deparse(substitute(mean_speed)), 
          ylab=deparse(substitute(sinuosity)),
          zlab=deparse(substitute(straight)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE,
          ellipsoid = TRUE,
          level=0.7, ellipsoid.alpha=0.6, id=FALSE,
          model.summary=FALSE)

scatter3d(x = data$mean_speed, y = data$sinuosity, z = data$straight, 
          groups = data1$probe,
          xlab=deparse(substitute(mean_speed)), 
          ylab=deparse(substitute(sinuosity)),
          zlab=deparse(substitute(straight)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE,
          ellipsoid = TRUE,
          level=0.5, ellipsoid.alpha=0.5, id=FALSE,
          model.summary=FALSE)

## Plot data in three factor space
scatter3d(x = data$SDDC, y = data$sinuosity, z = data$mean_speed, 
          groups = data1$probe,
          xlab=deparse(substitute(SDDC)), 
          ylab=deparse(substitute(Sinuosity)),
          zlab=deparse(substitute(mean_speed)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE, 
          ellipsoid = TRUE,
          level=0.8, ellipsoid.alpha=0.5, id=FALSE,
          model.summary=FALSE)

## Plot data in three factor space
scatter3d(x = data[data$probe != 'p36',]$sinuosity, 
          y = data[data$probe != 'p36',]$mean_speed, 
          z = data[data$probe != 'p36',]$straight, 
          groups = data[data$probe != 'p36',]$probe,
          xlab=deparse(substitute(Sinuosity)), 
          ylab=deparse(substitute(Mean_spees)),
          zlab=deparse(substitute(Straightness)), 
          axis.scales=TRUE, axis.ticks=TRUE,
          surface=FALSE, 
          ellipsoid = TRUE,
          level=0.7, ellipsoid.alpha=0.7, id=FALSE,
          model.summary=TRUE)

#library(lattice)
#cloud(data$mean_speed ~ data$straight*data$sinuosity)

#################################################
########### Fit statistical models ##############
#################################################

# fit generalized linear models
fit <- glm(probe ~ mean_speed + sinuosity + straight, data,
           family = "binomial")
anova(fit, test='Chisq')
newobject1 <- anova(fit, test='Chisq')
print(xtable(newobject1, type = "latex"), file = "filename1.tex")

fit <- glm(probe ~ mean_speed + sinuosity + straight, data,
           family = "quasibinomial")
anova(fit, test='Chisq')
newobject2 <- anova(fit, test='Chisq')
print(xtable(newobject2, type = "latex"), file = "filename2.tex")

fit <- glm(probe ~ mean_speed + length + straight + sinuosity + 
             mean_speed * sinuosity + 
             mean_speed * length + 
             straight * sinuosity +
             mean_speed * straight * sinuosity , data,
           family = "binomial")

library(MuMIn)
options(na.action = "na.fail") # change the default "na.omit" to prevent models
# from being fitted to different datasets in
# case of missing values.
globalmodel <- lm(probe ~ length + straight + sinuosity, data = data1)
combinations <- dredge(globalmodel)
print(combinations)
coefTable(combinations)

# GLM fit generalized linear models------------------------
fit <- glm(probe ~ mean_speed + sinuosity + straight, data,
           family = "binomial")
anova(fit, test='Chisq')
newobject1 <- anova(fit, test='Chisq')
print(xtable(newobject1, type = "latex"), file = "filename1.tex")

fit <- glm(probe ~ mean_speed + sinuosity + straight + SDDC, data,
           family = "binomial")
anova(fit, test='Chisq')
newobject1 <- anova(fit, test='Chisq')
print(xtable(newobject1, type = "latex"), file = "filename1.tex")

fit <- glm(probe ~ mean_speed + sinuosity + straight + SDDC + SDDC* straight, data,
           family = "binomial")
anova(fit, test='Chisq')
newobject1 <- anova(fit, test='Chisq')
print(xtable(newobject1, type = "latex"), file = "filename1.tex")

fit <- glm(probe ~ mean_speed + sinuosity + straight, data,
           family = "quasibinomial")
anova(fit, test='Chisq')
newobject2 <- anova(fit, test='Chisq')
print(xtable(newobject2, type = "latex"), file = "filename2.tex")

fit <- glm(probe ~ mean_speed + length + straight + sinuosity + 
             mean_speed * sinuosity + 
             mean_speed * length + 
             straight * sinuosity +
             mean_speed * straight * sinuosity , data,
           family = "binomial")

fit3 <- glm(probe ~ straight + sinuosity + SDDC + length, data, family = 'binomial')
summary(fit3)
anova(fit3, test = "Chisq")


# LM -----------

head(data)

fit_1 <- lm(sinuosity ~ ., data)
summary(fit_1)

fit_1 <- lm(sinuosity + length ~ ., data)
summary(fit_1)

fit_1 <- lm(length ~ ., data)
summary(fit_1)

fit_1 <- lm(straightness ~ ., data)
summary(fit_1)

cor.test(~ sinuosity + length, data)

vif(fit_1)

fit_2 <- lm(sinuosity ~ ., select(data, -SDDC))
summary(fit_2)

vif(fit_2)


library(MuMIn)
options(na.action = "na.fail") # change the default "na.omit" to prevent models
# from being fitted to different datasets in
# case of missing values.
globalmodel <- lm(probe ~ length + straight + sinuosity, data = data1)
combinations <- dredge(globalmodel)
print(combinations)
coefTable(combinations)

#################################################
######## Investigate for normality ##############
#################################################

# Density------------------------------------
ggdensity(data[data$probe=='p09',]$mean_speed,
          main = "Density plot of mean_speed in p09",
          xlab = "mean_speed")
ggdensity(data[data$probe=='p15',]$mean_speed,
          main = "Density plot of mean_speed in p15",
          xlab = "mean_speed")
ggdensity(data[data$probe=='p28',]$mean_speed,
          main = "Density plot of mean_speed in p28",
          xlab = "mean_speed")
ggdensity(data[data$probe=='p36',]$mean_speed,
          main = "Density plot of mean_speed in p36",
          xlab = "mean_speed")
# QQ plot--------------------------
# ggqqplot(data$mean_speed, main = 'mean_speed')
ggqqplot(data$mean_speed, main = 'mean_speed')
#ggqqplot(data$sd_speed, main = 'sd_speed')
ggqqplot(data$max_speed, main = 'max_speed')
#ggqqplot(data$min_speed, main = 'min_speed')
ggqqplot(data$length, main = 'length')
ggqqplot(data$distance, main = 'distance')
ggqqplot(data$straight, main = 'straight')
ggqqplot(data$square_displacement, main = 'square_displacement')
ggqqplot(data$sinuosity, main = 'sinuosity')
ggqqplot(data$emax, main = 'emax')
ggqqplot(data$DC, main = 'DC')
ggqqplot(data$SDDC, main = 'SDDC')

ggqqplot(data[data$probe=='p09',]$mean_speed, main = 'mean_speed')
ggqqplot(data[data$probe=='p15',]$mean_speed, main = 'mean_speed')
ggqqplot(data[data$probe=='p28',]$mean_speed, main = 'mean_speed')
ggqqplot(data[data$probe=='p36',]$mean_speed, main = 'mean_speed')

# qq plot
qqnorm(data$mean_speed, pch = 1, frame = FALSE)
qqline(data$mean_speed, col = "steelblue", lwd = 2)
qqplot(data$mean_speed)

# Shapiro test --------------------------
shapiro.test(data$SDDC)
shapiro.test(data[data$probe=='p09',]$SDDC)
shapiro.test(data[data$probe=='p15',]$SDDC)
shapiro.test(data[data$probe=='p28',]$SDDC)
shapiro.test(data[data$probe=='p36',]$SDDC)

# Test for normality----------------------------
shapiro.test(data[data$probe=='p09',]$mean_speed)
shapiro.test(data[data$probe=='p15',]$mean_speed)
shapiro.test(data[data$probe=='p28',]$mean_speed)
shapiro.test(data[data$probe=='p36',]$mean_speed) # data is not normally distributed

#################################################
### Compute mean and sd for track parameters ####
#################################################
# Mean and sd --------------------------
mean(data[data$probe=='p09',]$mean_speed)
sd(data[data$probe=='p09',]$mean_speed)

mean(data[data$probe=='p15',]$mean_speed)
sd(data[data$probe=='p15',]$mean_speed)

mean(data[data$probe=='p28',]$mean_speed)
sd(data[data$probe=='p28',]$mean_speed)

mean(data[data$probe=='p36',]$mean_speed)
sd(data[data$probe=='p36',]$mean_speed)
#
mean(data[data$probe=='p09',]$max_speed)
sd(data[data$probe=='p09',]$max_speed)

mean(data[data$probe=='p15',]$max_speed)
sd(data[data$probe=='p15',]$max_speed)

mean(data[data$probe=='p28',]$max_speed)
sd(data[data$probe=='p28',]$max_speed)

mean(data[data$probe=='p36',]$max_speed)
sd(data[data$probe=='p36',]$max_speed)
#
mean(data[data$probe=='p09',]$length)
sd(data[data$probe=='p09',]$length)

mean(data[data$probe=='p15',]$length)
sd(data[data$probe=='p15',]$length)

mean(data[data$probe=='p28',]$length)
sd(data[data$probe=='p28',]$length)

mean(data[data$probe=='p36',]$length)
sd(data[data$probe=='p36',]$length)
#
mean(data[data$probe=='p09',]$distance)
sd(data[data$probe=='p09',]$distance)

mean(data[data$probe=='p15',]$distance)
sd(data[data$probe=='p15',]$distance)

mean(data[data$probe=='p28',]$distance)
sd(data[data$probe=='p28',]$distance)

mean(data[data$probe=='p36',]$distance)
sd(data[data$probe=='p36',]$distance)
#
mean(data[data$probe=='p09',]$SDDC)
sd(data[data$probe=='p09',]$SDDC)

mean(data[data$probe=='p15',]$SDDC)
sd(data[data$probe=='p15',]$SDDC)

mean(data[data$probe=='p28',]$SDDC)
sd(data[data$probe=='p28',]$SDDC)

mean(data[data$probe=='p36',]$SDDC)
sd(data[data$probe=='p36',]$SDDC)
#
mean(data[data$probe=='p09',]$sinuosity)
sd(data[data$probe=='p09',]$sinuosity)

mean(data[data$probe=='p15',]$sinuosity)
sd(data[data$probe=='p15',]$sinuosity)

mean(data[data$probe=='p28',]$sinuosity)
sd(data[data$probe=='p28',]$sinuosity)

mean(data[data$probe=='p36',]$sinuosity)
sd(data[data$probe=='p36',]$sinuosity)
#

mean(data[data$probe=='p09',]$sinuosity)
sd(data[data$probe=='p09',]$sinuosity)

mean(data[data$probe=='p15',]$sinuosity)
sd(data[data$probe=='p15',]$sinuosity)

mean(data[data$probe=='p28',]$sinuosity)
sd(data[data$probe=='p28',]$sinuosity)

mean(data[data$probe=='p36',]$sinuosity)
sd(data[data$probe=='p36',]$sinuosity)
# Compute the analysis of variance------

res.aov <- aov(mean_speed ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)
# Straightness index
# Compute the analysis of variance
res.aov <- aov(straight ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)

# Compute the analysis of variance------
res.aov <- aov(sinuosity ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)
kruskal.test(data$straight ~ data$probe)


compare_means(straight ~ probe,  data = data, method = "kruskal.test")
write.xlsx(compare_means(mean_speed ~ probe,  data = data, method = "kruskal.test"),
           file = 'kruskal.test.mean_speed.xlsx')

compare_means(straight ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

write.xlsx(compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test"),
           file = 'wilcox.test.mean_speed.xlsx')


# Compute mean and SD-----------------------------
df.summary.sinuosity <- group_by(data, probe) %>%
  summarise(
    sd = sd(sinuosity, na.rm = TRUE),
    sinuosity = mean(sinuosity)
  )

df.summary.sinuosity


# Perform pairwise comparisons -------------------------
compare_means(SDDC ~ probe,  data = data, method = "t.test")
compare_means(SDDC ~ probe,  data = data, method = "wilcox.test")

compare_means(sinuosity ~ probe,  data = data, method = "t.test")
compare_means(sinuosity ~ probe,  data = data, method = "wilcox.test")

write.csv(compare_means(SDDC ~ probe,  data = data, method = "t.test"), file="t.test.SDDC.csv")
write.csv(compare_means(SDDC ~ probe,  data = data, method = "wilcox.test"), file="wilcox.test.SDDC.csv")

write.csv(compare_means(sinuosity ~ probe,  data = data, method = "t.test"), file="t.test.sinuosity.csv")
write.csv(compare_means(sinuosity ~ probe,  data = data, method = "wilcox.test"), file="wilcox.test.sinuosity.csv")

compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

# Visualize: Specify the comparisons you want color---------------------
my_comparisons <- list( c("p15","p09"), c("p28","p09"), c("p15", "p28"))

ggboxplot(data, x = "probe", y = "sinuosity")+ 
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 0)

#################################################
############### Make PCA and MFA ################
#################################################

# scree plot -------------------
Scree.Plot <- function(R,main="Scree Plot",sub=NULL){
  roots <- eigen(R)$values
  x <- 1:dim(R)[1]
  plot(x,roots,type="b",col='blue',ylab="Eigenvalue",
       xlab="Component Number",main=main,sub=sub) 
  abline(h=1,lty=2,col="red") }

Scree.Plot(cor(data[,-13]),main ="Scree Plot (Motility Data)")

# PCA ------------------------
#
res.pca <- PCA(data[,-13], graph = FALSE)
res.pca <- PCA(data[,-c(2, 4, 5,8, 10, 13)], graph = FALSE)
#print(res.pca)
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])

# Default plot
plot(res.pca, choix = "var")
fviz_pca_var(res.pca)
fviz_pca_var(res.pca, col.var="contrib") + theme_classic(base_size=14)
# Change the gradient color
fviz_pca_var(res.pca, col.var="contrib") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=5) + theme_classic(base_size=14)
# Control the transparency of variables using their contribution
# Possible values for the argument alpha.var are :
# "cos2", "contrib", "coord", "x", "y"
fviz_pca_var(res.pca, alpha.var="contrib") + theme_classic(base_size=14)
#Individuals factor map :
fviz_pca_ind(res.pca, label="none", col.ind = data$probe) + theme_classic(base_size=14)
#Change individual colors by groups :
fviz_pca_ind(res.pca, label="none", habillage=data$probe,
             addEllipses=TRUE, ellipse.level=0.95)

# MFA - Multiple Factor Analysis -----------------
res.mfa <- MFA(data, 
               group = c(4, 1, 1, 1, 1, 1, 1, 2, 1), 
               type = c("s", "s", "s","s", "s", "s","s", "s", "n" ),
               name.group = c("speed","length","distance",
                              "straightness", "square_displacement", 
                              "sinuosity", "emax" , "DC_SDDC", "probe"),
               #num.group.sup = c(1, 6),
               graph = FALSE)

eig.val <- get_eigenvalue(res.mfa)
head(eig.val)
fviz_screeplot(res.mfa)

fviz_mfa_var(res.mfa, "quanti.var", palette = "jco", 
             col.var.sup = "violet", repel = TRUE)

group <- get_mfa_var(res.mfa, "group")
group
# Coordinates of groups
head(group$coord)
# Cos2: quality of representation on the factore map
head(group$cos2)
# Contributions to the  dimensions
head(group$contrib)

fviz_mfa_var(res.mfa, "group")

# Contribution to the first dimension
fviz_contrib(res.mfa, "group", axes = 1)
# Contribution to the second dimension
fviz_contrib(res.mfa, "group", axes = 2)
# Contribution to the second dimension
fviz_contrib(res.mfa, "group", axes = 3)

quanti.var <- get_mfa_var(res.mfa, "quanti.var")
quanti.var 

# Coordinates
head(quanti.var$coord)
# Cos2: quality on the factore map
head(quanti.var$cos2)
# Contributions to the dimensions
head(quanti.var$contrib)

fviz_mfa_var(res.mfa, "quanti.var", palette = "jco", 
             col.var.sup = "violet", repel = TRUE)

fviz_mfa_var(res.mfa, "quanti.var", palette = "jco", 
             col.var.sup = "violet", repel = TRUE,
             geom = c("point", "text"), legend = "bottom")

# Contributions to dimension 1
fviz_contrib(res.mfa, choice = "quanti.var", axes = 1, top = 20,
             palette = "jco")
# Contributions to dimension 2
fviz_contrib(res.mfa, choice = "quanti.var", axes = 2, top = 20,
             palette = "jco")

fviz_mfa_var(res.mfa, "quanti.var", col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             col.var.sup = "violet", repel = TRUE,
             geom = c("point", "text"))

# Color by cos2 values: quality on the factor map
fviz_mfa_var(res.mfa, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             col.var.sup = "violet", repel = TRUE)

fviz_cos2(res.mfa, choice = "quanti.var", axes = 1)


ind <- get_mfa_ind(res.mfa)
ind

# Graph of partial axes
fviz_mfa_axes(res.mfa)

#################################################
############### Cluster analysis ################
#################################################

# Cluster--------
KM <- kmeans(data[,-c(2, 4, 5, 8, 10, 13)], 3, iter.max = 1000, algorithm = "Hartigan-Wong")
clusplot(data, KM$cluster, color = TRUE, 
         shade = FALSE, labels=1,
         main = 'Cluster Analysis for Motility Data')

data.stand <- scale(data[,-13])  # To standarize the variables

# K-Means
k.means.fit <- kmeans(data.stand, 5)
clusplot(data, k.means.fit$cluster, main='2D representation of the Cluster solution',
         color=TRUE, shade=TRUE,
         labels=1, lines=0)

attributes(k.means.fit)
# Centroids:
k.means.fit$centers
# Clusters:
k.means.fit$cluster
# Cluster size:
k.means.fit$size

wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot(data.stand, nc=6) 

# table(data[,1],k.means.fit$cluster)


#################################################
## Compute statistics and Plot results ##########
#################################################
# Summary -----------------------
summary(data[data$probe=='p09',])
summary(data[data$probe=='p15',])
summary(data[data$probe=='p28',])
summary(data[data$probe=='p36',])

# Kruskal test ---------

kruskal.test(data$mean_speed ~ data$probe)

# (01) Mean speed ----------------------
kruskal.test(data$mean_speed ~ data$probe)
compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
qplot(probe, mean_speed, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Mean speed") +
  labs(y = 'Micrometers per hour',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(100),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(90),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  
  geom_signif(y_position = c(110),
              xmin = c(1),
              xmax = c(3),
              annotation = "ns",
              tip_length = 0.04) +
  geom_signif(y_position = c(100),
              xmin = c(3),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04) +
  # annotation_custom(grob) +
  
  labs(title = "Mean speed", caption = "Kruskal-Wallis p-value < 2.2e-16")
  # stat_compare_means(label.y = 3)

kruskal.test(data$max_speed ~ data$probe)
compare_means(max_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (02) Max speed ----------------------
qplot(probe, max_speed, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Max speed" ) +
  labs(y = 'Micrometers',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(550),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(500),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  
  # annotation_custom(grob) +
  
  labs(title = "Max speed",
       # subtitle = "",
       caption = "Kruskal-Wallis p-value < 2.2e-16")


kruskal.test(data$length ~ data$probe)
compare_means(length ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (03) Length ----------------------
qplot(probe, length, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Length" ) +
  labs(y = 'Micrometers',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(2500),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(2200),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  geom_signif(y_position = c(2500),
              xmin = c(3),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  geom_signif(y_position = c(2800),
              xmin = c(1),
              xmax = c(3),
              annotation = "ns",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  # annotation_custom(grob) +
  
  labs(title = "Length",
       # subtitle = "",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data$distance ~ data$probe)
compare_means(distance ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (04) Distance ----------------------
qplot(probe, distance, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Distance" ) +
  labs(y = 'Micrometers',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1150),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04) +
  geom_signif(y_position = c(1000),
              xmin = c(1),
              xmax = c(2),
              annotation = "***",
              tip_length = 0.04) +
  geom_signif(y_position = c(1000),
              xmin = c(3),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1450),
              xmin = c(2),
              xmax = c(4),
              annotation = "ns",
              tip_length = 0.04) +
  geom_signif(y_position = c(1300),
              xmin = c(1),
              xmax = c(3),
              annotation = "ns",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  # xmin / xmax positions should match the x-axis labels' positions
  
  # annotation_custom(grob) +
  
  labs(title = "Distance",
       # subtitle = "",
       caption = "Kruskal-Wallis p-value = 2.251e-09")


kruskal.test(data$straight ~ data$probe)
compare_means(straight ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (05) Straightness ----------------------
qplot(probe, straight, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Straightness", ) +
  labs(y = 'Straightness index',
       x = "Cell passage")  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.1),
              xmin = c(1),
              xmax = c(2),
              annotation = "ns",
              tip_length = 0.04) +
  geom_signif(y_position = c(1.22),
              xmin = c(2),
              xmax = c(3),
              annotation = "ns",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.35),
              xmin = c(3),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04)+
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.5),
              xmin = c(1),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  
  # annotation_custom(grob) +
  
  labs(title = "Straightness",
       # subtitle = "",
       caption = "Kruskal-Wallis p-value = 4.142e-09")

kruskal.test(data$sinuosity ~ data$probe)
compare_means(sinuosity ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (6) Sinuosity------------------------
# Add jitter and change fill color by probe
qplot(probe, sinuosity, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Random search path tortuosity", ) +
  labs(y = 'Sinuosity index',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.9),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04)+
  geom_signif(y_position = c(0.9),
              xmin = c(3),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04)+
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.1),
              xmin = c(1),
              xmax = c(3),
              annotation = "*",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  
  # annotation_custom(grob) +
  
  labs(title = "Sinuosity",
       # subtitle = "",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data[data$probe != "p36",]$emax ~ data[data$probe != "p36",]$probe)
compare_means(emax ~ probe,  data = data[data$probe != "p36",], method = "wilcox.test") # pairwise comparisons
# (7) Emax ------------------------
# Emax-a is a dimensionless, scale-independent measure of the maximum possible expected displacement
qplot(probe, emax, data = data[data$probe != "p36",],
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Trajectory straightness index, E-max", ) +
  labs(y = 'emax',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  geom_signif(y_position = c(15),
              xmin = c(1),
              xmax = c(3),
              annotation = "*",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  labs(title = "Straightness index, E-max",
       caption = "Kruskal-Wallis p-value = 0.04776")

# selected passages
qplot(probe, emax, data = data[data$probe != "p36",],
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Trajectory straightness index, E-max", ) +
  labs(y = 'emax',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  labs(title = "Trajectory straightness index, E-max",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data$sd_speed ~ data$probe)
compare_means(sd_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (8) sd_speed ------------------------
qplot(probe, sd_speed, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "sd_speed", ) +
  labs(y = 'sd_speed',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  labs(title = "sd_speed",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data$square_displacement ~ data$probe)
compare_means(square_displacement ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (9) square_displacement ------------------------
qplot(probe, square_displacement, data = data[data$probe != "p36",],
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "square_displacement", ) +
  labs(y = 'square_displacement',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  labs(title = "square_displacement",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data$DC ~ data$probe)
compare_means(DC ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (10) DC ------------------------
qplot(probe, DC, data = data[data$probe != "p36",],
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "DC", ) +
  labs(y = 'DC',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  labs(title = "DC",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data[data$probe != "p36",]$SDDC ~ data[data$probe != "p36",]$probe)
compare_means(SDDC ~ probe,  data = data[data$probe != "p36",], method = "wilcox.test") # pairwise comparisons
# (11) SDDC ------------------------
qplot(probe, SDDC, data = data[data$probe != "p36",],
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "SDDC", ) +
  labs(y = 'SDDC',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  labs(title = "SDDC",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data$min_speed ~ data$probe)
compare_means(min_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (12) min_speed ------------------------
qplot(probe, min_speed, data = data[data$probe != "p36",],
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "min_speed", ) +
  labs(y = 'min_speed',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  labs(title = "min_speed",
       caption = "Kruskal-Wallis p-value < 2.2e-16")


# jitter plot--------------------------------
ggplot(data, aes(probe, straight)) +
  geom_bar(stat = "identity", data = df.summary.straight, size=1.2,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.3),  size=1, color = 'red') +
  geom_errorbar(
    aes(ymin = straight-sd, ymax = straight+sd), color = 'black', size=1.2,
    data = df.summary.straight, width = 0.2) +
  ggtitle('straight')+
  ylim(0, 1) +
  ggtitle("MSC-WJ1, 24h trajecroty straightness") +
  labs(y="straight", x = "Passage")  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.7),
              xmin = c(1),
              xmax = c(2),
              annotation = "ns",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.8),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04)+
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.9),
              xmin = c(1),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04) +
  theme(
    # Change axis lines
    axis.line = element_line(size = 1),
    # Change axis ticks text labels: font color, size and face
    axis.text.x = element_text(face = "bold",
                               size = 12, angle = 0),     # Change x axis tick labels only
    axis.text.y = element_text(face = "bold",
                               size = 12, angle = 0),     # Change y axis tick labels only
    # Change axis ticks line: font color, size, linetype and length
    axis.ticks = element_line(),      # Change ticks line fo all axes
    axis.ticks.x = element_line(),    # Change x axis ticks only
    axis.ticks.y = element_line(),    # Change y axis ticks only
    axis.ticks.length = unit(3, "pt"), # Change the length of tick marks
    legend.title = element_text(color = "black", size = 15),
    legend.text = element_text(color = "black", size = 15))



# Plot means
plotmeans(mean_speed ~ probe, data = data, frame = FALSE, 
          mean.labels=FALSE, connect=TRUE,
          n.label=TRUE, text.n.label="n = ",
          xlab = "Passages", ylab = "Micrometers per hour",
          main="WJMSC-1 cells 24 h track Mean Speed, 
          \nMean Plot with 95% CI")

plotmeans(sinuosity ~ probe, data = data, frame = FALSE, 
          mean.labels=FALSE, connect=TRUE,
          n.label=TRUE, text.n.label="n = ",
          xlab = "Passages", ylab = "Sinuosity index",
          main="WJMSC-1 cells 24 h tracks Sinuosity, 
          \nMean Plot with 95% CI")

plotmeans(straight ~ probe, data = data, frame = FALSE, 
          mean.labels=FALSE, connect=TRUE,
          n.label=TRUE, text.n.label="n=",
          xlab = "Passages", ylab = "Straightness index",
          main="WJMSC-1 24 h tracks Straightness, 
          \nMean Plot with 95% CI")

#  ------------

#################################################
################ Beta-galactozidase #############
#################################################

BinomCI(224, 3724, conf.level = 0.95)
BinomCI(494, 2404, conf.level = 0.95)
BinomCI(475, 1103, conf.level = 0.95)

#################################################
################ Mean Comparison ################
#################################################

mean(data[data$probe=='p09',]$sinuosity)
mean(data[data$probe=='p15',]$sinuosity)
mean(data[data$probe=='p28',]$sinuosity)
mean(data[data$probe=='p36',]$sinuosity)

mean(data[data$probe=='p09',]$mean_speed)
mean(data[data$probe=='p15',]$mean_speed)
mean(data[data$probe=='p28',]$mean_speed)
mean(data[data$probe=='p36',]$mean_speed)


#################################################
################## Cumulative plots #############
#################################################
# Aggregate keypoint data -----------------------
all <- data.frame(probe = c(9, 15, 28, 36), 
                  colocM9act = c(0.6359818, 0.4870214, 0.7352667, 0.6310091), 
                  colocA4act = c(NA, NA, NA, NA), 
                  colocRhoAnuc = c(0.2182182, 0.5042364, 0.0989375, -0.01028462), 
                  speed = c(38.32615, 25.02812, 38.05941, 18.25367), 
                  sinuosity = c(0.2649515, 0.3288146, 0.2509206, 0.2192818), 
                  bgal = c(0.06, 0.21, 0.43, 0.53) )
# all
# ggplots ----------------------------
ggplot(data = all, aes(x = probe)) +
  geom_line(aes(y = speed, colour = "Mean_speed"), lwd=2, alpha = .99, lty=2) +
  scale_colour_manual("", 
                      breaks = c("Mean_speed"),
                      values = c( 
                        "Mean_speed"="red"
                      )) +
  xlab("Passages") +
  scale_y_continuous("Micrometers per hour", limits = c(0,50)) + 
  scale_x_continuous(breaks = c(9,15,28,36), labels = c("9", "15", "28", "36")) +
  labs(title="") + 
  geom_vline(xintercept = 15, lwd=1, lty=3) +
  geom_vline(xintercept = 28, lwd=1, lty=3) +
  theme_classic2()

ggplot(data = all, aes(x = probe)) +
  geom_line(aes(y = sinuosity, colour = "Sinuosity"), lwd=2, alpha = .99, lty=2) +
  scale_colour_manual("", 
                      breaks = c("Sinuosity"),
                      values = c( 
                        "Sinuosity"="blue"
                      )) +
  xlab("Passages") +
  scale_y_continuous("Sinuosity", limits = c(0.2,.35)) + 
  scale_x_continuous(breaks = c(9,15,28,36), labels = c("9", "15", "28", "36")) +
  labs(title="") + 
  geom_vline(xintercept = 15, lwd=1, lty=3) +
  geom_vline(xintercept = 28, lwd=1, lty=3) +
  theme_classic2()

ggplot(data = all, aes(x = probe)) +
  geom_line(aes(y = bgal, colour = "Beta-galactosidase"), lwd=2, alpha = .99, lty=2) +
  scale_colour_manual("", 
                      breaks = c("Beta-galactosidase"),
                      values = c( 
                        "Beta-galactosidase"="violetred4"
                      )) +
  xlab("Passages") +
  scale_y_continuous("", limits = c(0,.6)) + 
  scale_x_continuous(breaks = c(9,15,28,36), labels = c("9", "15", "28", "36")) +
  labs(title="") + 
  geom_vline(xintercept = 15, lwd=1, lty=3) +
  geom_vline(xintercept = 28, lwd=1, lty=3) +
  theme_classic2()



ggplot(data = all, aes(x = probe)) +
  geom_line(aes(y = colocRhoAnuc, colour = "RhoA_nucleus"), lwd=2, alpha = .9) +
  geom_line(aes(y = colocM9act, colour = "Myosin-9_F-actin"), lwd=2, alpha = .77) +
  scale_colour_manual("", 
                      breaks = c( 
                        "Myosin-9_F-actin",
                        "RhoA_nucleus"),
                      values = c( 
                        "RhoA_nucleus" = "deeppink",
                        "Myosin-9_F-actin" = "springgreen4"
                      )) +
  xlab("Passages") +
  scale_y_continuous("", limits = c(-0.2,1)) + 
  scale_x_continuous(breaks = c(9,15,28,36), labels = c("9", "15", "28", "36")) +
  labs(title="") + 
  geom_vline(xintercept = 15, lwd=1, lty=3) +
  geom_vline(xintercept = 28, lwd=1, lty=3) +
  theme_classic2()


#################################################
## Function includes time rediscretization ######
#################################################
traj_analysis_rediscret <- function(input) {  # modified for 24 frames-data

  data <- read.csv(input)

  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))

  for (i in unique(data$Track)) {

    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         #times = c(1:24))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")

    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 23/3600/24)

    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)

    trj <- TrajScale(trj, 1.3158 / 1, "micrometer")

    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory

    par(mar=c(5,5,5,5))

    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    trj <- TrajResampleTime(trj, 901)

    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)

    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)


    # Plot rediscretized trajectory in red
    lines(trj, col = "#FF0000A0", lwd = 2)
    points(trj, type = 'p', col = "#FF0000A0", pch = 16)

    legend("topright", c("Original", "Resampled"), col = c("black", "red"),
           lwd = 2, inset = c(0.01, 0.02))


    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    TrajMeanVectorOfTurningAngles(trj) # eturns the mean vector of the turning angles
    TrajAngles(trj) # Returns the turning angles (radians) of a trajectory
    TrajMeanVelocity(trj) # Returns the mean velocity vector of the trajectory (or a portion)
    TrajGetTimeUnits(trj) # Returns the temporal units of a trajectory

    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)

    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj),
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           probe = 'p36')

    head(traj_params)

  }

  # print(traj_params)
  write.csv(traj_params, file = 'traj_p36.csv')
  tracks <<- traj_params
  return(traj_params) }
# 24 frames, TrajRediscretize

# Radarchart -----------
# help(fmsb)
library(fmsb)

data2 <- data.frame(data)

radarchart(data2[,2], axistype=1, seg=5, plty=1, vlabels=c(), 
           title="(axis=1with specified vlabels)", vlcex=0.5)

radarchart(dat, axistype=2, pcol=topo.colors(3), plty=1, pdensity=c(5, 10, 30), 
           pangle=c(10, 45, 120), pfcol=topo.colors(3), 
           title="(topo.colors, fill, axis=2)")

radarchart(dat, axistype=3, pty=32, plty=1, axislabcol="grey", na.itp=FALSE,
           title="(no points, axis=3, na.itp=FALSE)")

radarchart(dat, axistype=1, plwd=1:5, pcol=1, centerzero=TRUE, 
           seg=4, caxislabels=c("worst", "", "", "", "best"),
           title="(use lty and lwd but b/w, axis=1,\n centerzero=TRUE, with centerlabels)")
par(op)
