# cell shape analisys

# Actin cytoskeleton fractality analysis
# If you use this script in your publications, please cite:
# Bobkov D, Polyanskaya A, Musorina A, Lomert E, Shabelnikov S, Poljanskaya G. 
# Replicative senescence in MSCWJ-1 human umbilical cord mesenchymal stem cells is marked by characteristic changes in motility, 
# cytoskeletal organization, and RhoA localization. Mol Biol Rep. 2020;:1–17. DOI: 10.1007/s11033-020-05476-6
# mailto: bobkov@incras.ru
# Statistical analysis 
# Import libraries
# Load libraries-------------
library(GGally)
library(trajr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(xlsx)
library(plyr)
library(dplyr) 
library(ggpubr)
library(car)
library(stringi)
library(Hmisc)
library(gplots)
library(PMCMRplus)
library(dunn.test)
library(DescTools)
library(ggsignif)
library(grid)
library(xtable)
library(parallel)
library(nortest)
library(KScorrect)
################## Load data from csv file:
data <- read.csv('wj1shape.csv', sep = ';', dec = '.',
                 header=TRUE, stringsAsFactors=FALSE)

datafrac <- read.csv('wj1fracmod.csv', sep = ';', dec = ',')

data <- read.csv('../frac/lmwj1/wj1frac.csv', sep = ';', dec = ',')
data <- read.csv('../frac/lmwj1/wj1fracmod.csv', sep = ';', dec = ',')
#data <- read.csv('adhfrac.csv', sep = ';', dec = ',')
data$Dhat <- remove_outliers(data$Dhat)

data$Area <- remove_outliers(data$Area)

plotmeans(data$Dhat ~ data$passage,
          xlab = c('Passage number'),
          ylab = c('F-actin LCFD'),
          use.t = TRUE)
###
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

##
## png()
#par(mfrow = c(1, 2))
boxplot(data$Area)
## dev.off()
data$Area <- remove_outliers(data$Area)
data$Circ <- remove_outliers(data$Circ)
data$Perim <- remove_outliers(data$Perim)

data$Dhat <- remove_outliers(data$Dhat)
datafrac$Dhat <- remove_outliers(datafrac$Dhat)

boxplot(data$Area ~ data$passage)

boxplot(data$Dhat ~ data$passage)


data$passage <- as.factor(data$passage)

glimpse(data)

###################
ggpairs(data)

ggcorr(data, palette = "RdBu", label = TRUE)
##################

plotmeans(data$Dhat ~ data$passage,
          xlab = c('Passage number'),
          ylab = c('F-actin LCFD'),
          ylim = c(1.7, 1.95),
          use.t = TRUE) 

plotmeans(data$Area ~ data$passage,
          xlab = c('Passage number'),
          ylab = c('Cell Area, square microns'),
          ylim = c(3000, 10000),
          use.t = TRUE)


plotmeans(data$Circ ~ data$passage,
          xlab = c('Passage number'),
          ylab = c('Cell circularity'),
          ylim = c(.3, .5),
          use.t = TRUE)

plotmeans(data$Perim ~ data$passage,
          xlab = c('Passage number'),
          ylab = c('Cell Perimeter, microns'),
          ylim = c(300, 700),
          use.t = TRUE)

plotmeans(data$AR ~ data$passage,
          xlab = c('Passage number'),
          ylab = c('Cell AR'),
          use.t = TRUE)

plotmeans(data$Round ~ data$passage,
          xlab = c('Passage number'),
          ylab = c('Cell Round'),
          use.t = TRUE)

plotmeans(data$Solidity ~ data$passage,
          xlab = c('Passage number'),
          ylab = c('Cell Solidity'),
          use.t = TRUE)


mean(data$Area, na.rm = TRUE)


install.packages("plotrix")

library("plotrix")

std.error(c(1,2,3,4))

# Mean and 95% CI calculation
t.test(data[data$passage=='9',]$Area)$"conf.int"
t.test(data[data$passage=='9',]$Area)
(4297.736 - 5202.035)/-2  #452.1495

t.test(data[data$passage=='15',]$Area)$"conf.int"
t.test(data[data$passage=='15',]$Area)
(4247.510 - 5217.284)/-2 #484.887

t.test(data[data$passage=='28',]$Area)$"conf.int"
t.test(data[data$passage=='28',]$Area)
(4835.616 - 5701.752)/-2 #433.068

t.test(data[data$passage=='36',]$Area)$"conf.int"
t.test(data[data$passage=='36',]$Area)
(7522.276 - 8741.374)/-2 #609.549

t.test(data[data$passage=='38',]$Area)$"conf.int"
t.test(data[data$passage=='38',]$Area)
(8359.333 - 9801.028)/-2 #720.8475


# Mean and 95% CI calculation
t.test(data[data$passage=='9',]$Circ)$"conf.int"
t.test(data[data$passage=='9',]$Circ)
(0.4219042 - 0.4841546)/-2  #0.0311252

t.test(data[data$passage=='15',]$Circ)$"conf.int"
t.test(data[data$passage=='15',]$Circ)
(0.3831267 - 0.4435255)/-2 #0.0301994

t.test(data[data$passage=='28',]$Circ)$"conf.int"
t.test(data[data$passage=='28',]$Circ)
(0.3528737 - 0.4091055)/-2 #0.0281159

t.test(data[data$passage=='36',]$Circ)$"conf.int"
t.test(data[data$passage=='36',]$Circ)
(0.3661626 -0.4228374)/-2 #0.0283374

t.test(data[data$passage=='38',]$Circ)$"conf.int"
t.test(data[data$passage=='38',]$Circ)
(0.3410044 - 0.3934738)/-2 #0.0262347

# Mean and 95% CI calculation
t.test(data[data$passage=='9',]$Dhat)$"conf.int"
t.test(data[data$passage=='9',]$Dhat)
(1.792890- 1.843771)/-2  #0.0254405

t.test(data[data$passage=='15',]$Dhat)$"conf.int"
t.test(data[data$passage=='15',]$Dhat)
(1.730628 -1.793550)/-2 #0.031461

t.test(data[data$passage=='28',]$Dhat)$"conf.int"
t.test(data[data$passage=='28',]$Dhat)
(1.812878 - 1.868111)/-2 #0.0276165

t.test(data[data$passage=='36',]$Dhat)$"conf.int"
t.test(data[data$passage=='36',]$Dhat)
(1.841569   -1.900131)/-2 #0.029281

t.test(data[data$passage=='38',]$Dhat)$"conf.int"
t.test(data[data$passage=='38',]$Dhat)
(1.830497 -1.906781)/-2 #0.038142

mean(data[data$passage=='9',]$Area, na.rm = TRUE)
sd(data[data$passage=='9',]$Area, na.rm = TRUE)


mean(data[data$passage=='15',]$Area, na.rm = TRUE)
sd(data[data$passage=='15',]$Area, na.rm = TRUE)

mean(data[data$passage=='28',]$Area, na.rm = TRUE)
sd(data[data$passage=='28',]$Area, na.rm = TRUE)

mean(data[data$passage=='36',]$Area, na.rm = TRUE)
sd(data[data$passage=='36',]$Area, na.rm = TRUE)

mean(data[data$passage=='38',]$Area, na.rm = TRUE)
sd(data[data$passage=='38',]$Area, na.rm = TRUE)



mean(data[data$passage=='9',]$Circ, na.rm = TRUE)
sd(data[data$passage=='9',]$Circ, na.rm = TRUE)

mean(data[data$passage=='15',]$Circ, na.rm = TRUE)
sd(data[data$passage=='15',]$Circ, na.rm = TRUE)

mean(data[data$passage=='28',]$Circ, na.rm = TRUE)
sd(data[data$passage=='28',]$Circ, na.rm = TRUE)

mean(data[data$passage=='36',]$Circ, na.rm = TRUE)
sd(data[data$passage=='36',]$Circ, na.rm = TRUE)

mean(data[data$passage=='38',]$Circ, na.rm = TRUE)
sd(data[data$passage=='38',]$Circ, na.rm = TRUE)


mean(data[data$passage=='9',]$Dhat, na.rm = TRUE)
sd(data[data$passage=='9',]$Dhat, na.rm = TRUE)
std.error(data[data$passage=='9',]$Dhat, na.rm = TRUE)

mean(data[data$passage=='15',]$Dhat, na.rm = TRUE)
sd(data[data$passage=='15',]$Dhat, na.rm = TRUE)

mean(data[data$passage=='28',]$Dhat, na.rm = TRUE)
sd(data[data$passage=='28',]$Dhat, na.rm = TRUE)

mean(data[data$passage=='36',]$Dhat, na.rm = TRUE)
sd(data[data$passage=='36',]$Dhat, na.rm = TRUE)

mean(data[data$passage=='38',]$Dhat, na.rm = TRUE)
sd(data[data$passage=='38',]$Dhat, na.rm = TRUE)
######## Investigate for normality ##############
#################################################
shapiro.test(data[data$passage=='9',]$Area)
shapiro.test(data[data$passage=='15',]$Area)
shapiro.test(data[data$passage=='28',]$Area)
shapiro.test(data[data$passage=='36',]$Area)
shapiro.test(data[data$passage=='38',]$Area)

# Density------------------------------------
ggdensity(data[data$passage=='9',]$Area,
          main = "Density plot of 9 in W1",
          xlab = "Area")

ggdensity(data[data$passage=='15',]$Area,
          main = "Density plot of 15 in W1",
          xlab = "Area")

ggdensity(data[data$passage=='28',]$Area,
          main = "Density plot of 28 in W1",
          xlab = "Area")

ggdensity(data[data$passage=='36',]$Area,
          main = "Density plot of 36 in W1",
          xlab = "Area")

ggdensity(data[data$passage=='38',]$Area,
          main = "Density plot of 38 in W1",
          xlab = "Area")


# Density------------------------------------
ggdensity(data[data$passage=='9',]$Dhat,
          main = "Density plot of 9 in W1",
          xlab = "Area")

ggdensity(data[data$passage=='15',]$Dhat,
          main = "Density plot of 15 in W1",
          xlab = "Area")

ggdensity(data[data$passage=='28',]$Dhat,
          main = "Density plot of 28 in W1",
          xlab = "Area")

ggdensity(data[data$passage=='36',]$Dhat,
          main = "Density plot of 36 in W1",
          xlab = "Area")

ggdensity(data[data$passage=='38',]$Dhat,
          main = "Density plot of 38 in W1",
          xlab = "Area")



ggplot() + 
  geom_density(data=data, 
               aes(x=passage, 
                   group=passage, 
                   fill=passage),
               alpha=0.5, adjust=1) + 
  xlab("Passage") +
  ylab("Area Density")


# Use custom palette
ggdensity(data, x = "Area",
          add = "mean", rug = TRUE,
          color = "passage", fill = "passage")


# Q-Q plot: Q-Q plot (or quantile-quantile plot) draws the correlation between a given sample and the normal distribution. 
ggqqplot(data$Area, main = 'Area')
ggqqplot(data$Circ, main = 'Circ')

shapiro.test(data[data$passage=='9',]$Dhat)
shapiro.test(data[data$passage=='15',]$Dhat)
shapiro.test(data[data$passage=='28',]$Dhat)
shapiro.test(data[data$passage=='36',]$Dhat)
shapiro.test(data[data$passage=='38',]$Dhat)

############

ggplot(data, aes(x = passage, y = Area)) + 
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  # theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  # geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  # add zero level
  
  labs(title = "Area in WJ1",
       caption = "cell area") +
  
  geom_smooth(method='lm', formula= y ~ x) +# add linear regression
  
  scale_x_continuous(name = "Cell passage",
                     breaks = c(seq(7, 40, 1)), 
                     limits = c(7, 40))  
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   size = 9)) +
  #  
  scale_y_continuous(name = "Dlc", 
                     breaks = c(seq(0, 1000, 5)), 
                     limits = c(0, 1000))


fit <- lm(Area ~ passage, data = data)

summary(fit)


# Perform pairwise comparisons
compare_means(Area ~ passage,  data = data, method = "anova")
compare_means(Area ~ passage,  data = data, method = "kruskal.test")
compare_means(Area ~ passage,  data = data, method = "t.test")
compare_means(Area ~ passage,  data = data, method = "wilcox.test")

###################

# Perform pairwise comparisons
compare_means(Dhat ~ passage,  data = data, method = "anova")
compare_means(Dhat ~ passage,  data = data, method = "kruskal.test")
compare_means(Dhat ~ passage,  data = data, method = "t.test")
compare_means(Dhat ~ passage,  data = data, method = "wilcox.test")

############

ggplot(data, aes(x = passage, y = Circ)) + 
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  # theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  # geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  # add zero level
  
  labs(title = "WJ1",
       caption = "cell circularity") +
  
  geom_smooth(method='lm', formula= y ~ x) +# add linear regression
  
  scale_x_continuous(name = "Cell passage",
                     breaks = c(seq(7, 40, 1)), 
                     limits = c(7, 40)) 

theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                 size = 9)) +
  #  
  scale_y_continuous(name = "Dlc", 
                     breaks = c(seq(0, 1000, 5)), 
                     limits = c(0, 1000))


fit <- lm(Dhat ~ passage, data = data)

summary(fit)


shapiro.test(data[data$passage=='9',]$Circ)
shapiro.test(data[data$passage=='15',]$Circ)
shapiro.test(data[data$passage=='28',]$Circ)
shapiro.test(data[data$passage=='36',]$Circ)
shapiro.test(data[data$passage=='38',]$Circ)


shapiro.test(data[data$passage=='9',]$Dhat)
shapiro.test(data[data$passage=='15',]$Dhat)
shapiro.test(data[data$passage=='28',]$Dhat)
shapiro.test(data[data$passage=='36',]$Dhat)
shapiro.test(data[data$passage=='38',]$Dhat)


# Perform pairwise comparisons
compare_means(Circ ~ passage,  data = data, method = "anova")
compare_means(Circ ~ passage,  data = data, method = "kruskal.test")
compare_means(Circ ~ passage,  data = data, method = "t.test")
compare_means(Circ ~ passage,  data = data, method = "wilcox.test")


# function to produce summary statistics (mean and +/- sd), as required for ggplot2
data_summary <- function(x) {
  mu <- mean(x)
  sigma1 <- mu-sd(x)
  sigma2 <- mu+sd(x)
  return(c(y=mu,ymin=sigma1,ymax=sigma2))
}

# function for computing mean, DS, max and min values
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}



summary(dunnettTest(data$Area, data$passage, alternative = "two.sided"))
# plot 1 Mean +- SD
#ggsave("plot1.jpg", width = 12, height = 12, units = "cm")
data$passage <- as.factor(data$passage)


ggplot(data, aes(x = passage, y = Area)) +
  
  
  ylim(c(0, 18050)) +
  
  theme_classic(base_size=14) +
  
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", 
               aes(width=0.7),
               position=position_dodge(width=1.5), size=0.8) +
  
  geom_jitter(position = position_jitter(.15),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
 # scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Cell area, mkm^2',
       x = "Cell passage") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 1,
                                   size = 12, face="bold",
                                   colour="black" ), 
        axis.text.y = element_text(angle = 0, hjust = 1,
                                   size = 11, face = 'bold',
                                   colour="black" )) +
  
  geom_signif(y_position = c(15000),
              xmin = c(1),
              xmax = c(2),
              annotation = "NS", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold") +
  
  geom_signif(y_position = c(17700),
              xmin = c(1),
              xmax = c(3),
              annotation = "*", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold") +
  
  geom_signif(y_position = c(16300),
              xmin = c(2),
              xmax = c(3),
              annotation = "*", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold") +
  
  geom_signif(y_position = c(15000),
              xmin = c(3),
              xmax = c(4),
              annotation = "****", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold")

# Perform pairwise comparisons
compare_means(Area ~ passage,  data = data, method = "anova")
compare_means(Area ~ passage,  data = data, method = "kruskal.test")
compare_means(Area ~ passage,  data = data, method = "t.test")
compare_means(Area ~ passage,  data = data, method = "wilcox.test")


ggplot(data, aes(x = passage, y = Circ)) +
  
  
  ylim(c(0.1, .9)) +
  
  theme_classic(base_size=14) +
  
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", 
               aes(width=0.7),
               position=position_dodge(width=1.5), size=0.8) +
  
  geom_jitter(position = position_jitter(.15),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  # scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Cell circularity',
       x = "Cell passage") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 1,
                                   size = 12, face="bold",
                                   colour="black" ), 
        axis.text.y = element_text(angle = 0, hjust = 1,
                                   size = 11, face = 'bold',
                                   colour="black" ))
  
  geom_signif(y_position = c(15000),
              xmin = c(1),
              xmax = c(2),
              annotation = "NS", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold") +
  
  geom_signif(y_position = c(17700),
              xmin = c(1),
              xmax = c(3),
              annotation = "*", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold") +
  
  geom_signif(y_position = c(16300),
              xmin = c(2),
              xmax = c(3),
              annotation = "*", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold") +
  
  geom_signif(y_position = c(15000),
              xmin = c(3),
              xmax = c(4),
              annotation = "****", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold")

# Perform pairwise comparisons
compare_means(Circ ~ passage,  data = data, method = "anova")
compare_means(Circ ~ passage,  data = data, method = "kruskal.test")
compare_means(Circ ~ passage,  data = data, method = "t.test")
compare_means(Circ ~ passage,  data = data, method = "wilcox.test")



# Perform pairwise comparisons
compare_means(Dhat ~ passage,  data = datafrac, method = "anova")
compare_means(Dhat ~ passage,  data = datafrac, method = "kruskal.test")
compare_means(Dhat ~ passage,  data = datafrac, method = "t.test")
compare_means(Dhat ~ passage,  data = datafrac, method = "wilcox.test")
#################################################
################ Plot 3D images #################
#################################################

## Plot data in three factor space
scatter3d(x = data$Area, y = data$Circ, z = data$Round, 
          groups = data$passage,
          xlab=deparse(substitute(Area)), 
          ylab=deparse(substitute(Circ)),
          zlab=deparse(substitute(Round)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE, 
          ellipsoid = TRUE,
          level=0.8, ellipsoid.alpha=0.5, id=FALSE,
          model.summary=FALSE)











ggplot(datafrac, aes(x = passage, y = Dhat, group = passage)) +
  
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", 
               aes(width=0.7),
               position=position_dodge(width=1.5), size=0.8) +
  
  geom_jitter(position = position_jitter(.15),
              cex = .9,
              shape = 16)
  
  
  ylim(c(0, 2)) +
  
  theme_classic(base_size=14) +
  
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", 
               aes(width=0.7),
               position=position_dodge(width=1.5), size=0.8) +
  
  geom_jitter(position = position_jitter(.15),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  # scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Cell circularity',
       x = "Cell passage") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 1,
                                   size = 12, face="bold",
                                   colour="black" ), 
        axis.text.y = element_text(angle = 0, hjust = 1,
                                   size = 11, face = 'bold',
                                   colour="black" ))

geom_signif(y_position = c(15000),
            xmin = c(1),
            xmax = c(2),
            annotation = "NS", textsize = 4, family="serif",
            tip_length = 0.04, size = 1, face="bold") +
  
  geom_signif(y_position = c(17700),
              xmin = c(1),
              xmax = c(3),
              annotation = "*", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold") +
  
  geom_signif(y_position = c(16300),
              xmin = c(2),
              xmax = c(3),
              annotation = "*", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold") +
  
  geom_signif(y_position = c(15000),
              xmin = c(3),
              xmax = c(4),
              annotation = "****", textsize = 4, family="serif",
              tip_length = 0.04, size = 1, face="bold")
