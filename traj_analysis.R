# Cell movement trajectory analysis by Dan Bobkov, 2019 # dan.bobkov@gmail.com

# Load libraries-------------
library(GGally)
library(trajr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(xlsx)
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
  }
}
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
  return(traj_params)
  
} # 96 frames
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
  return(traj_params)
  
} # 96 frames
# p36----------
# modified for 24 frames-data, includes time rediscritization:
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
  return(traj_params)
  
} # 24 frames, TrajRediscretize
# Choose dir ---------
file_list_p09 <- list.files(path = , choose.dir(default = "", 
                                            caption = "Select folder"),
                        pattern = "csv", 
                        all.files = FALSE,
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = FALSE, include.dirs = FALSE,
                        no.. = FALSE)
# Choose dir
file_list_p15 <- list.files(path = , choose.dir(default = "", 
                                            caption = "Select folder"),
                        pattern = "csv", 
                        all.files = FALSE,
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = FALSE, include.dirs = FALSE,
                        no.. = FALSE)
# Choose dir
file_list_p36 <- list.files(path = , choose.dir(default = "", 
                                            caption = "Select folder"),
                        pattern = "csv", 
                        all.files = FALSE,
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = FALSE, include.dirs = FALSE,
                        no.. = FALSE)
# Call to function--------------
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
# Start scan for p36
for (file_name in file_list_p36) {
  traj_analysis_p36(file_name)
  tracks_p36 <- rbind(tracks_p36, tracks)
}
# Merge all tracks-----------------------
alltracks <- rbind(tracks_p09, tracks_p15, tracks_p36) # collect all tracks
summary(alltracks)
write.csv(alltracks, file = 'alltracks.csv') #save results
# Order probe levels--------------------
alltracks$probe <- as.factor(alltracks$probe)
alltracks$probe <- ordered(alltracks$probe,
                      levels = c("p09", "p15", "p36"))
# Set time to hours, remove tracks and mean_angle columns----------------
all.h <- alltracks
head(all.h)
all.h <- cbind(all.h[,c(6,7,8,9)]*3600, all.h[,c(2,3,4,5,10,11,12,13,15)])
data <- all.h
head(data)
# Plot all-in-one-------------------------
# plot(data)
png()
ggpairs(data)
ggcorr(data, palette = "RdBu", label = TRUE)
dev.off()


# Summary -----------------------
summary(data[data$probe=='p09',])
summary(data[data$probe=='p15',])
summary(data[data$probe=='p36',])
### Stat analysis ###


# (1) Mean speed------------------------
# Compute the analysis of variance------
res.aov <- aov(mean_speed ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)
# Compute mean and SD-----------------------------
df.summary.mean_speed <- group_by(data, probe) %>%
  summarise(
    sd = sd(mean_speed, na.rm = TRUE),
    mean_speed = mean(mean_speed)
  )

df.summary.mean_speed

mean(data[data$probe=='p09',]$mean_speed)
sd(data[data$probe=='p09',]$mean_speed)

mean(data[data$probe=='p15',]$mean_speed)
sd(data[data$probe=='p15',]$mean_speed)

mean(data[data$probe=='p36',]$mean_speed)
sd(data[data$probe=='p36',]$mean_speed)
# Density------------------------------------
ggdensity(data[data$probe=='p09',]$mean_speed, 
          main = "Density plot of mean_speed in p09",
          xlab = "mean_speed")
ggdensity(data[data$probe=='p15',]$mean_speed, 
          main = "Density plot of mean_speed in p15",
          xlab = "mean_speed")
ggdensity(data[data$probe=='p36',]$mean_speed, 
          main = "Density plot of mean_speed in p36",
          xlab = "mean_speed")
# ggqqplot(data$mean_speed, main = 'mean_speed')------------------
ggqqplot(data[data$probe=='p09',]$mean_speed, main = 'mean_speed')
ggqqplot(data[data$probe=='p15',]$mean_speed, main = 'mean_speed')
ggqqplot(data[data$probe=='p36',]$mean_speed, main = 'mean_speed')
# Test for normality
shapiro.test(data[data$probe=='p09',]$mean_speed)
shapiro.test(data[data$probe=='p15',]$mean_speed)
shapiro.test(data[data$probe=='p36',]$mean_speed) # data is not normally distributed

kruskal.test(data$mean_speed ~ data$probe)

compare_means(mean_speed ~ probe,  data = data, method = "kruskal.test")
write.xlsx(compare_means(mean_speed ~ probe,  data = data, method = "kruskal.test"), 
           file = 'kruskal.test.mean_speed.xlsx')

compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

write.xlsx(compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test"),
           file = 'wilcox.test.mean_speed.xlsx')


# Bar plot with signifiers ----------------------------
ggplot(df.summary.mean_speed, aes(probe, mean_speed)) +
  geom_bar(stat = "identity", fill = 'gray', 
           color = "black", size= 1, show.legend=TRUE) +
  geom_errorbar(aes(ymin = mean_speed-sd, ymax = mean_speed+sd), width = 0.2, size=1) +
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
    axis.ticks.length = unit(3, "pt") # Change the length of tick marks
  ) +
  geom_point() +
  ylim(0, 70) + 
  ggtitle("Mean speed, MSCWJ1, 24h") + 
  labs(y="Mean speed, micrometers per hour", x = "Passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(60),
              xmin = c(1),
              xmax = c(2),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(45),
              xmin = c(2),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04)
# jitter plot--------------------------------
ggplot(data, aes(probe, mean_speed)) +
  geom_bar(stat = "identity", data = df.summary.mean_speed, size=1.2,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.3),  size=1, color = 'red') + 
  geom_errorbar(
    aes(ymin = mean_speed-sd, ymax = mean_speed+sd), color = 'black', size=1.2,
    data = df.summary.mean_speed, width = 0.2) + 
    ggtitle('mean_speed')+
  ylim(0, 110) +
  ggtitle("MSC-WJ1, 24h Mean speed") + 
  labs(y="Micrometers per hour", x = "Passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(100),
              xmin = c(1),
              xmax = c(2),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(80),
              xmin = c(2),
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
  
# Add jitter and change fill color by probe----------------------
qplot(probe, mean_speed, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Mean speed, MSCWJ1, 24h", ) + 
  labs(y = 'Micrometers per hour',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(100),
              xmin = c(1),
              xmax = c(2),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(80),
              xmin = c(2),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04) + theme_bw()

# (2) Max speed------------------------
# Compute the analysis of variance-----
res.aov <- aov(max_speed ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)
# Compute mean and SD----
df.summary.max_speed <- group_by(data, probe) %>%
  summarise(
    sd = sd(max_speed, na.rm = TRUE),
    max_speed = mean(max_speed)
  )

df.summary.max_speed


mean(data[data$probe=='p09',]$max_speed)
sd(data[data$probe=='p09',]$max_speed)

mean(data[data$probe=='p15',]$max_speed)
sd(data[data$probe=='p15',]$max_speed)

mean(data[data$probe=='p36',]$max_speed)
sd(data[data$probe=='p36',]$max_speed)
# Density------------------------------------
ggdensity(data[data$probe=='p09',]$max_speed, 
          main = "Density plot of max_speed in p09",
          xlab = "max_speed")
ggdensity(data[data$probe=='p15',]$max_speed, 
          main = "Density plot of max_speed in p15",
          xlab = "max_speed")
ggdensity(data[data$probe=='p36',]$max_speed, 
          main = "Density plot of max_speed in p36",
          xlab = "max_speed")
#---
ggqqplot(data[data$probe=='p09',]$max_speed, main = 'max_speed')
ggqqplot(data[data$probe=='p15',]$max_speed, main = 'max_speed')
ggqqplot(data[data$probe=='p36',]$max_speed, main = 'max_speed')

# Test for normality
shapiro.test(data[data$probe=='p09',]$max_speed)
shapiro.test(data[data$probe=='p15',]$max_speed)
shapiro.test(data[data$probe=='p36',]$max_speed) # data is not normally distributed

kruskal.test(data$max_speed ~ data$probe)

compare_means(max_speed ~ probe,  data = data, method = "kruskal.test")
write.xlsx(compare_means(max_speed ~ probe,  data = data, method = "kruskal.test"), 
           file = 'kruskal.test.max_speed.xlsx')

compare_means(max_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

write.xlsx(compare_means(max_speed ~ probe,  data = data, method = "wilcox.test"),
           file = 'wilcox.test.max_speed.xlsx')


# Bar plot with signifiers ------------------
ggplot(df.summary.max_speed, aes(probe, max_speed)) +
  geom_bar(stat = "identity", fill = 'gray', 
           color = "black", size= 1, show.legend=TRUE) +
  geom_errorbar(aes(ymin = max_speed-sd, ymax = max_speed+sd), width = 0.2, size=1) +
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
    axis.ticks.length = unit(3, "pt") # Change the length of tick marks
  ) +
  geom_point() +
  ylim(0, 300) + 
  ggtitle("MSCWJ1 24h-trajectory max_speed") + 
  labs(y="max_speed, micrometers per hour", x = "Passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(300),
              xmin = c(1),
              xmax = c(2),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(250),
              xmin = c(2),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04)
# jitter plot--------------------------------
ggplot(data, aes(probe, max_speed)) +
  geom_bar(stat = "identity", data = df.summary.max_speed, size=1.2,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.3),  size=1, color = 'blue') + 
  geom_errorbar(
    aes(ymin = max_speed-sd, ymax = max_speed+sd), color = 'black', size=1.2,
    data = df.summary.max_speed, width = 0.2) + 
  ggtitle('Length')+
  ylim(0, 400) +
  ggtitle("MSC-WJ1, 24h trajectory max_speed") + 
  labs(y="Micrometers", x = "Passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(400),
              xmin = c(1),
              xmax = c(2),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(300),
              xmin = c(2),
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

# Add jitter and change fill color by probe----------------------
qplot(probe, max_speed, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Max speed on 24h-trajectory, MSCWJ1" ) + 
  labs(y = 'Micrometers',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(500),
              xmin = c(1),
              xmax = c(2),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(450),
              xmin = c(2),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04) + theme_bw()


# (3) Length   ------------------------
# Compute the analysis of variance-----
res.aov <- aov(length ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)
# Compute mean and SD----
df.summary.length <- group_by(data, probe) %>%
  summarise(
    sd = sd(length, na.rm = TRUE),
    length = mean(length)
  )

df.summary.length


mean(data[data$probe=='p09',]$length)
sd(data[data$probe=='p09',]$length)

mean(data[data$probe=='p15',]$length)
sd(data[data$probe=='p15',]$length)

mean(data[data$probe=='p36',]$length)
sd(data[data$probe=='p36',]$length)
# Density------------------------------------
ggdensity(data[data$probe=='p09',]$length, 
          main = "Density plot of length in p09",
          xlab = "length")
ggdensity(data[data$probe=='p15',]$length, 
          main = "Density plot of length in p15",
          xlab = "length")
ggdensity(data[data$probe=='p36',]$length, 
          main = "Density plot of length in p36",
          xlab = "length")
#---
ggqqplot(data[data$probe=='p09',]$length, main = 'length')
ggqqplot(data[data$probe=='p15',]$length, main = 'length')
ggqqplot(data[data$probe=='p36',]$length, main = 'length')

# Test for normality
shapiro.test(data[data$probe=='p09',]$length)
shapiro.test(data[data$probe=='p15',]$length)
shapiro.test(data[data$probe=='p36',]$length) # data is not normally distributed

kruskal.test(data$length ~ data$probe)

compare_means(length ~ probe,  data = data, method = "kruskal.test")
write.xlsx(compare_means(length ~ probe,  data = data, method = "kruskal.test"), 
           file = 'kruskal.test.length.xlsx')

compare_means(length ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

write.xlsx(compare_means(length ~ probe,  data = data, method = "wilcox.test"),
           file = 'wilcox.test.length.xlsx')


# Bar plot with signifiers ------------------
ggplot(df.summary.length, aes(probe, length)) +
  geom_bar(stat = "identity", fill = 'gray', 
           color = "black", size= 1, show.legend=TRUE) +
  geom_errorbar(aes(ymin = length-sd, ymax = length+sd), width = 0.2, size=1) +
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
    axis.ticks.length = unit(3, "pt") # Change the length of tick marks
  ) +
  geom_point() +
  ylim(0, 1500) + 
  ggtitle("MSCWJ1 24h-trajectory length") + 
  labs(y="Length, micrometers", x = "Passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1400),
              xmin = c(1),
              xmax = c(2),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1300),
              xmin = c(2),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04)
# jitter plot--------------------------------
ggplot(data, aes(probe, length)) +
  geom_bar(stat = "identity", data = df.summary.length, size=1.2,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.3),  size=1, color = 'blue') + 
  geom_errorbar(
    aes(ymin = length-sd, ymax = length+sd), color = 'black', size=1.2,
    data = df.summary.length, width = 0.2) + 
  ggtitle('Length')+
  ylim(0, 2700) +
  ggtitle("MSC-WJ1, 24h trajectory length") + 
  labs(y="Micrometers", x = "Passage") +
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

# Add jitter and change fill color by probe----------------------
qplot(probe, length, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Length of 24h-trajectory, MSCWJ1" ) + 
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
              tip_length = 0.04) + theme_bw()

# (4) Distance ------------------------
# Compute the analysis of variance
res.aov <- aov(distance ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)

# Compute mean and SD
df.summary.distance <- group_by(data, probe) %>%
  summarise(
    sd = sd(distance, na.rm = TRUE),
    distance = mean(distance)
  )

df.summary.distance


mean(data[data$probe=='p09',]$distance)
sd(data[data$probe=='p09',]$distance)

mean(data[data$probe=='p15',]$distance)
sd(data[data$probe=='p15',]$distance)

mean(data[data$probe=='p36',]$distance)
sd(data[data$probe=='p36',]$iu)

# Density
ggdensity(data[data$probe=='p09',]$distance, 
          main = "Density plot of distance in p09",
          xlab = "distance")
ggdensity(data[data$probe=='p15',]$distance, 
          main = "Density plot of distance in p15",
          xlab = "distance")
ggdensity(data[data$probe=='p36',]$distance, 
          main = "Density plot of distance in p36",
          xlab = "distance")

ggqqplot(data[data$probe=='p09',]$distance, main = 'distance')
ggqqplot(data[data$probe=='p15',]$distance, main = 'distance')
ggqqplot(data[data$probe=='p36',]$distance, main = 'distance')

# Test for normality
shapiro.test(data[data$probe=='p09',]$distance)
shapiro.test(data[data$probe=='p15',]$distance)
shapiro.test(data[data$probe=='p36',]$distance) # data is not normally distributed

kruskal.test(data$distance ~ data$probe)

compare_means(distance ~ probe,  data = data, method = "kruskal.test")
write.xlsx(compare_means(distance ~ probe,  data = data, method = "kruskal.test"), 
           file = 'kruskal.test.distance.xlsx')

compare_means(distance ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

write.xlsx(compare_means(distance ~ probe,  data = data, method = "wilcox.test"),
           file = 'wilcox.test.distance.xlsx')


# Bar plot with signifiers ------------------
ggplot(df.summary.distance, aes(probe, distance)) +
  geom_bar(stat = "identity", fill = 'gray', 
           color = "black", size= 1, show.legend=TRUE) +
  geom_errorbar(aes(ymin = distance-sd, ymax = distance+sd), width = 0.2, size=1) +
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
    axis.ticks.length = unit(3, "pt") # Change the length of tick marks
  ) +
  geom_point() +
  ylim(0, 600) + 
  ggtitle("Distance") + 
  labs(y="Distance, micrometers", x = "Passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(600),
              xmin = c(1),
              xmax = c(2),
              annotation = "***", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(500),
              xmin = c(1),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(400),
              xmin = c(2),
              xmax = c(3),
              annotation = "ns", 
              tip_length = 0.04)
# jitter plot--------------------------------
ggplot(data, aes(probe, distance)) +
  geom_bar(stat = "identity", data = df.summary.distance, size=1.2,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.3),  size=1, color = 'blue') + 
  geom_errorbar(
    aes(ymin = distance-sd, ymax = distance+sd), color = 'black', size=1.2,
    data = df.summary.distance, width = 0.2) + 
  ggtitle('Distance')+
  ylim(0, 700) +
  ggtitle("Distance") + 
  labs(y="Micrometers", x = "Passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(600),
              xmin = c(1),
              xmax = c(2),
              annotation = "***", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(500),
              xmin = c(1),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(400),
              xmin = c(2),
              xmax = c(3),
              annotation = "ns", 
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

# Add jitter and change fill color by probe----------------------
qplot(probe, distance, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Distance, 24h, MSCWJ1" ) + 
  labs(y = 'Micrometers',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1000),
              xmin = c(1),
              xmax = c(2),
              annotation = "***", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(900),
              xmin = c(1),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(800),
              xmin = c(2),
              xmax = c(3),
              annotation = "ns", 
              tip_length = 0.04) + theme_bw()



# (5) Straightness---------------------
# Straightness index
# Compute the analysis of variance
res.aov <- aov(straight ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)
# Compute mean and SD

df.summary.straight <- group_by(data, probe) %>%
  summarise(
    sd = sd(straight, na.rm = TRUE),
    straight = mean(straight)
  )

df.summary.straight


mean(data[data$probe=='p09',]$length)
sd(data[data$probe=='p09',]$length)

mean(data[data$probe=='p15',]$length)
sd(data[data$probe=='p15',]$length)

mean(data[data$probe=='p36',]$length)
sd(data[data$probe=='p36',]$length)


# Density
ggdensity(data[data$probe=='p09',]$mean_speed, 
          main = "Density plot of mean_speed in p09",
          xlab = "mean_speed")
ggdensity(data[data$probe=='p15',]$mean_speed, 
          main = "Density plot of mean_speed in p15",
          xlab = "mean_speed")
ggdensity(data[data$probe=='p36',]$mean_speed, 
          main = "Density plot of mean_speed in p36",
          xlab = "mean_speed")

# ggqqplot(data$mean_speed, main = 'mean_speed')
ggqqplot(data[data$probe=='p09',]$mean_speed, main = 'mean_speed')
ggqqplot(data[data$probe=='p15',]$mean_speed, main = 'mean_speed')
ggqqplot(data[data$probe=='p36',]$mean_speed, main = 'mean_speed')

# Test for normality
shapiro.test(data[data$probe=='p09',]$mean_speed)
shapiro.test(data[data$probe=='p15',]$mean_speed)
shapiro.test(data[data$probe=='p36',]$mean_speed) # data is not normally distributed

kruskal.test(data$straight ~ data$probe)

compare_means(straight ~ probe,  data = data, method = "kruskal.test")
write.xlsx(compare_means(mean_speed ~ probe,  data = data, method = "kruskal.test"), 
           file = 'kruskal.test.mean_speed.xlsx')

compare_means(straight ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

write.xlsx(compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test"),
           file = 'wilcox.test.mean_speed.xlsx')



# Compute the analysis of variance------
res.aov <- aov(straight ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)
# Compute mean and SD-----------------------------
df.summary.straight <- group_by(data, probe) %>%
  summarise(
    sd = sd(straight, na.rm = TRUE),
    straight = mean(straight)
  )

df.summary.straight

mean(data[data$probe=='p09',]$straight)
sd(data[data$probe=='p09',]$straight)

mean(data[data$probe=='p15',]$straight)
sd(data[data$probe=='p15',]$straight)

mean(data[data$probe=='p36',]$straight)
sd(data[data$probe=='p36',]$straight)
# Density------------------------------------
ggdensity(data[data$probe=='p09',]$straight, 
          main = "Density plot of straight in p09",
          xlab = "straight")
ggdensity(data[data$probe=='p15',]$straight, 
          main = "Density plot of straight in p15",
          xlab = "straight")
ggdensity(data[data$probe=='p36',]$straight, 
          main = "Density plot of straight in p36",
          xlab = "straight")

ggqqplot(data[data$probe=='p09',]$straight, main = 'straight')
ggqqplot(data[data$probe=='p15',]$straight, main = 'straight')
ggqqplot(data[data$probe=='p36',]$straight, main = 'straight')
# Test for normality
shapiro.test(data[data$probe=='p09',]$straight)
shapiro.test(data[data$probe=='p15',]$straight)
shapiro.test(data[data$probe=='p36',]$straight) # data is not normally distributed

kruskal.test(data$straight ~ data$probe)

compare_means(straight ~ probe,  data = data, method = "kruskal.test")
write.xlsx(compare_means(straight ~ probe,  data = data, method = "kruskal.test"), 
           file = 'kruskal.test.straight.xlsx')

compare_means(straight ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

write.xlsx(compare_means(straight ~ probe,  data = data, method = "wilcox.test"),
           file = 'wilcox.test.straight.xlsx')


# Bar plot with signifiers ----------------------------
ggplot(df.summary.straight, aes(probe, straight)) +
  geom_bar(stat = "identity", fill = 'gray', 
           color = "black", size= 1, show.legend=TRUE) +
  geom_errorbar(aes(ymin = straight-sd, ymax = straight+sd), width = 0.2, size=1) +
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
    axis.ticks.length = unit(3, "pt") # Change the length of tick marks
  ) +
  geom_point() +
  ylim(0, 0.9) + 
  ggtitle("straightness index, MSCWJ1, 24h") + 
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
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.9),
              xmin = c(1),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04)
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

# Add jitter and change fill color by probe----------------------
qplot(probe, straight, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Straightness index, MSCWJ1, 24h", ) + 
  labs(y = 'straight',
       x = "Cell passage")  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.1),
              xmin = c(1),
              xmax = c(2),
              annotation = "ns", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.25),
              xmin = c(2),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04)+
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.4),
              xmin = c(1),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04) + theme_bw()



# (6) Sinuosity------------------------
# Compute the analysis of variance------
res.aov <- aov(sinuosity ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)
# Compute mean and SD-----------------------------
df.summary.sinuosity <- group_by(data, probe) %>%
  summarise(
    sd = sd(sinuosity, na.rm = TRUE),
    sinuosity = mean(sinuosity)
  )

df.summary.sinuosity

mean(data[data$probe=='p09',]$sinuosity)
sd(data[data$probe=='p09',]$sinuosity)

mean(data[data$probe=='p15',]$sinuosity)
sd(data[data$probe=='p15',]$sinuosity)

mean(data[data$probe=='p36',]$sinuosity)
sd(data[data$probe=='p36',]$sinuosity)
# Density------------------------------------
ggdensity(data[data$probe=='p09',]$sinuosity, 
          main = "Density plot of sinuosity in p09",
          xlab = "sinuosity")
ggdensity(data[data$probe=='p15',]$sinuosity, 
          main = "Density plot of sinuosity in p15",
          xlab = "sinuosity")
ggdensity(data[data$probe=='p36',]$sinuosity, 
          main = "Density plot of sinuosity in p36",
          xlab = "sinuosity")

ggqqplot(data[data$probe=='p09',]$sinuosity, main = 'sinuosity')
ggqqplot(data[data$probe=='p15',]$sinuosity, main = 'sinuosity')
ggqqplot(data[data$probe=='p36',]$sinuosity, main = 'sinuosity')
# Test for normality
shapiro.test(data[data$probe=='p09',]$sinuosity)
shapiro.test(data[data$probe=='p15',]$sinuosity)
shapiro.test(data[data$probe=='p36',]$sinuosity) # data is not normally distributed

kruskal.test(data$sinuosity ~ data$probe)

compare_means(sinuosity ~ probe,  data = data, method = "kruskal.test")
write.xlsx(compare_means(sinuosity ~ probe,  data = data, method = "kruskal.test"), 
           file = 'kruskal.test.sinuosity.xlsx')

compare_means(sinuosity ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

write.xlsx(compare_means(sinuosity ~ probe,  data = data, method = "wilcox.test"),
           file = 'wilcox.test.sinuosity.xlsx')


# Bar plot with signifiers ----------------------------
ggplot(df.summary.sinuosity, aes(probe, sinuosity)) +
  geom_bar(stat = "identity", fill = 'gray', 
           color = "black", size= 1, show.legend=TRUE) +
  geom_errorbar(aes(ymin = sinuosity-sd, ymax = sinuosity+sd), width = 0.2, size=1) +
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
    axis.ticks.length = unit(3, "pt") # Change the length of tick marks
  ) +
  geom_point() +
  ylim(0, 0.9) + 
  ggtitle("Sinuosity, MSCWJ1, 24h") + 
  labs(y="Sinuosity", x = "Passage")  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.82),
              xmin = c(1),
              xmax = c(2),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.75),
              xmin = c(2),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.9),
              xmin = c(1),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04)
# jitter plot--------------------------------
ggplot(data, aes(probe, sinuosity)) +
  geom_bar(stat = "identity", data = df.summary.sinuosity, size=1.2,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.3),  size=1, color = 'red') + 
  geom_errorbar(
    aes(ymin = sinuosity-sd, ymax = sinuosity+sd), color = 'black', size=1.2,
    data = df.summary.sinuosity, width = 0.2) + 
  ggtitle('sinuosity')+
  ylim(0, 1) +
  ggtitle("MSC-WJ1, 24h trajecroty sinuosity") + 
  labs(y="sinuosity", x = "Passage")  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.9),
              xmin = c(1),
              xmax = c(2),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.8),
              xmin = c(2),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04)+
    # xmin / xmax positions should match the x-axis labels' positions
    geom_signif(y_position = c(0.7),
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

# Add jitter and change fill color by probe----------------------
qplot(probe, sinuosity, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Sinuosity, MSCWJ1, 24h", ) + 
  labs(y = 'Sinuosity',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.8),
              xmin = c(1),
              xmax = c(2),
              annotation = "****", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.9),
              xmin = c(2),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04)+
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1),
              xmin = c(1),
              xmax = c(3),
              annotation = "****", 
              tip_length = 0.04) + theme_bw()




#--------
# Perform pairwise comparisons

compare_means(distance ~ probe,  data = data, method = "anova")
compare_means(distance ~ probe,  data = data, method = "kruskal.test")
compare_means(distance ~ probe,  data = data, method = "t.test")

compare_means(distance ~ probe,  data = data, method = "wilcox.test")
compare_means(straight ~ probe,  data = data, method = "wilcox.test")
compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test")

write.xlsx(compare_means(ros ~ probe,  data = data, method = "kruskal.test"), 
           file = 'kruskal.test.xlsx')

write.xlsx(compare_means(ros ~ probe, data = data, method = "anova"),
           file = 'anova.xlsx')

write.xlsx(compare_means(ros ~ probe,  data = data, method = "t.test"),
           file = 't.test.xlsx')







ggqqplot(data$length, main = 'lenght')
ggqqplot(data$distance, main = 'distance')
ggqqplot(data$square_displacement, main = 'square_displacement')

ggqqplot(data$sd_speed, main = 'sd_speed')
ggqqplot(data$max_speed, main = 'max_speed')
ggqqplot(data$min_speed, main = 'min_speed')
ggqqplot(data$sinuosity, main = 'sinuosity')
ggqqplot(data$emax, main = 'emax')
ggqqplot(data$DC, main = 'DC')
ggqqplot(data$SDDC, main = 'SDDC')

colnames(data)



#--------

df.summary.length <- group_by(data, probe) %>%
  summarise(
    sd = sd(length, na.rm = TRUE),
    length = mean(length)
  )

df.summary.length

df.length <- data

ggplot(df.length, aes(probe, length)) +
  geom_bar(stat = "identity", data = df.summary.length,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") + 
  geom_errorbar(
    aes(ymin = length-sd, ymax = length+sd),
    data = df.summary.length, width = 0.2) + ggtitle('length')

#--------

df.summary.distance <- group_by(data, probe) %>%
  summarise(
    sd = sd(distance, na.rm = TRUE),
    distance = mean(distance)
  )

df.summary.distance

df.distance <- data

ggplot(df.distance, aes(probe, distance)) +
  geom_bar(stat = "identity", data = df.summary.distance,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") + 
  geom_errorbar(
    aes(ymin = distance-sd, ymax = distance+sd),
    data = df.summary.distance, width = 0.2) #+ ggtitle('distance')+
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(580),
              xmin = c(1),
              xmax = c(2),
              annotation = "**", 
              tip_length = 0.04)# +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(500),
              xmin = c(1),
              xmax = c(3),
              annotation = "***", 
              tip_length = 0.04)

#--------

df.summary.square_displacement <- group_by(data, probe) %>%
  summarise(
    sd = sd(square_displacement, na.rm = TRUE),
    square_displacement = mean(square_displacement)
  )

df.summary.square_displacement

df.square_displacement <- data

ggplot(df.square_displacement, aes(probe, square_displacement)) +
  geom_bar(stat = "identity", data = df.summary.square_displacement,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") + 
  geom_errorbar(
    aes(ymin = square_displacement-sd, ymax = square_displacement+sd),
    data = df.summary.square_displacement, width = 0.2) + ggtitle('square_displacement') 


#--------
df.summary.sd_speed <- group_by(data, probe) %>%
  summarise(
    sd = sd(sd_speed, na.rm = TRUE),
    sd_speed = mean(sd_speed)
  )

df.summary.sd_speed

df.sd_speed <- data

ggplot(df.sd_speed, aes(probe, sd_speed)) +
  geom_bar(stat = "identity", data = df.summary.sd_speed,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") + 
  geom_errorbar(
    aes(ymin = sd_speed-sd, ymax = sd_speed+sd),
    data = df.summary.sd_speed, width = 0.2) + ggtitle('sd_speed')

#--------

df.summary.max_speed <- group_by(data, probe) %>%
  summarise(
    sd = sd(max_speed, na.rm = TRUE),
    max_speed = mean(max_speed)
  )

df.summary.max_speed

df.max_speed <- data

ggplot(df.max_speed, aes(probe, max_speed)) +
  geom_bar(stat = "identity", data = df.summary.max_speed,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") + 
  geom_errorbar(
    aes(ymin = max_speed-sd, ymax = max_speed+sd),
    data = df.summary.max_speed, width = 0.2) + ggtitle('max_speed')

#--------

df.summary.min_speed <- group_by(data, probe) %>%
  summarise(
    sd = sd(min_speed, na.rm = TRUE),
    min_speed = mean(min_speed)
  )

df.summary.min_speed

df.min_speed <- data

ggplot(df.min_speed, aes(probe, min_speed)) +
  geom_bar(stat = "identity", data = df.summary.min_speed,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") + 
  geom_errorbar(
    aes(ymin = min_speed-sd, ymax = min_speed+sd),
    data = df.summary.min_speed, width = 0.2) + ggtitle('min_speed')


#--------

df.summary.sinuosity <- group_by(data, probe) %>%
  summarise(
    sd = sd(sinuosity, na.rm = TRUE),
    sinuosity = mean(sinuosity)
  )

df.summary.sinuosity

df.sinuosity <- data


ggplot(df.sinuosity, aes(probe, sinuosity)) +
  geom_bar(stat = "identity", data = df.summary.sinuosity,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") + 
  geom_errorbar(
    aes(ymin = sinuosity-sd, ymax = sinuosity+sd),
    data = df.summary.sinuosity, width = 0.2) + ggtitle('sinuosity')

#--------





df.summary.emax <- group_by(data, probe) %>%
  summarise(
    sd = sd(emax, na.rm = TRUE),
    emax = mean(emax)
  )

df.summary.emax

df.emax <- data

ggplot(df.emax, aes(probe, emax)) +
  geom_bar(stat = "identity", data = df.summary.emax,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") + 
  geom_errorbar(
    aes(ymin = emax-sd, ymax = emax+sd),
    data = df.summary.emax, width = 0.2) + ggtitle('emax')

#--------

df.summary.DC <- group_by(data, probe) %>%
  summarise(
    sd = sd(DC, na.rm = TRUE),
    DC = mean(DC)
  )

df.summary.DC

df.DC <- data

ggplot(df.DC, aes(probe, DC)) +
  geom_bar(stat = "identity", data = df.summary.DC,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") + 
  geom_errorbar(
    aes(ymin = DC-sd, ymax = DC+sd),
    data = df.summary.DC, width = 0.2) + ggtitle('DC')

#--------

df.summary.SDDC <- group_by(data, probe) %>%
  summarise(
    sd = sd(SDDC, na.rm = TRUE),
    SDDC = mean(SDDC)
  )

df.summary.SDDC

df.SDDC <- data

ggplot(df.SDDC, aes(probe, SDDC)) +
  geom_bar(stat = "identity", data = df.summary.SDDC,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") + 
  geom_errorbar(
    aes(ymin = SDDC-sd, ymax = SDDC+sd),
    data = df.summary.SDDC, width = 0.2) + ggtitle('SDDC')

#--------


#--------
# par(mfrow = c(2,2))
png()

qplot(probe, length, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), #log = "y",
      main = "MSCWJ1, 24h Path length") + theme_classic2()

qplot(probe, distance, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), #log = "y",
      main = "WJ1, 24h")

qplot(probe, square_displacement, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), log = "y",
      main = "WJ1, 24h")

qplot(probe, sd_speed, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), log = "y",
      main = "WJ1, 24h")

qplot(probe, max_speed, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), log = "y",
      main = "WJ1, 24h")

qplot(probe, min_speed, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), log = "y",
      main = "WJ1, 24h")

qplot(probe, sinuosity, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), log = "y",
      main = "WJ1, 24h")

qplot(probe, DC, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), log = "y",
      main = "WJ1, 24h")

qplot(probe, SDDC, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), log = "y",
      main = "WJ1, 24h")

qplot(probe, emax, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), log = "y",
      main = "WJ1, 24h")

qplot(probe, straight, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.6), #log = "y",
      main = "Straightness, WJ1, 24h")

# dev.off()




#--------
all <- read.csv('alltracks.csv')
all <- data

b <- ggplot(all, aes(x = probe, y = emax))

b + geom_point()
b + geom_jitter()

b + geom_boxplot(aes(color = probe))
b + geom_boxplot(aes(fill = probe))
b + geom_boxplot(aes(fill = probe)) + scale_fill_grey()

b + geom_boxplot() + coord_flip()
b + geom_boxplot(notch = TRUE)

b + stat_boxplot()

b + geom_violin()
b + geom_line()

b + geom_dotplot(binaxis = "y", stackdir = "center",
                 stackratio = 1, dotsize = 0.2)




