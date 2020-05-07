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

######################################
data <- read.csv('wj1frac.csv', sep = ';', dec = ',')

glimpse(data)
# Q-Q plot: Q-Q plot (or quantile-quantile plot) draws the correlation between a given sample and the normal distribution. 
ggqqplot(data$Dhat, main = 'Dlc')
ggqqplot(data[data$age=='young',]$Dhat, main = 'Dlc')
ggqqplot(data[data$age=='middle',]$Dhat, main = 'Dlc')
ggqqplot(data[data$age=='old',]$Dhat, main = 'Dlc')
# normality test
# perform the Shapiro-Wilk test of normality for one variable (univariate):
shapiro.test(data$Dhat)

mean(data[data$age=='young',]$Dhat)
sd(data[data$age=='young',]$Dhat)

mean(data[data$age=='middle',]$Dhat)
sd(data[data$age=='middle',]$Dhat)

mean(data[data$age=='old',]$Dhat)
sd(data[data$age=='old',]$Dhat)

mean(data$Dhat)
sd(data$Dhat)

# Density plot: the density plot provides a visual judgment about whether the distribution is bell shaped.
ggdensity(data$Dhat, 
          main = "Density plot of Dlc",
          xlab = "Dlc") + #
  
  scale_x_continuous(name = "D",
                     breaks = c(seq(1.4, 2, .2)), 
                     limits = c(1.4, 2))

ggdensity(data[data$age=='young',]$Dhat, 
          main = "Density plot of Dlc young",
          xlab = "Dlc")+ #
  
  scale_x_continuous(name = "D",
                     breaks = c(seq(1.4, 2, .2)), 
                     limits = c(1.4, 2))

ggdensity(data[data$age=='middle',]$Dhat, 
          main = "Density plot of Dlc middle",
          xlab = "Dlc")+ #
  
  scale_x_continuous(name = "D",
                     breaks = c(seq(1.4, 2, .2)), 
                     limits = c(1.4, 2))

ggdensity(data[data$age=='old',]$Dhat, 
          main = "Density plot of Dlc old",
          xlab = "Dlc")+ #
  
  scale_x_continuous(name = "D",
                     breaks = c(seq(1.4, 2, .2)), 
                     limits = c(1.4, 2))
# Perform pairwise comparisons
compare_means(Dhat ~ passage,  data = data, method = "anova")
compare_means(Dhat ~ passage,  data = data, method = "kruskal.test")
compare_means(Dhat ~ passage,  data = data, method = "t.test")
compare_means(Dhat ~ passage,  data = data, method = "wilcox.test")

# Kruskal Wallis Test One Way Anova by Ranks
kruskal.test(data$Dhat ~ data$passage)
#--boxplot---

qplot(passage, Dhat, data = data, aes(group = passage),
      geom = c("boxplot"), alpha = I(0.3),
      main = "F-actin local fractal dimension" ) + 
  labs(y = 'Dlc',
       x = "Cell passage") +   theme_classic(base_size=14) + # add linear regression
  
  scale_x_continuous(name = "Cell passage",
                     breaks = c(seq(7, 18, 1)), 
                     limits = c(7, 20))
# xmin / xmax positions should match the x-axis labels' positions
#geom_signif(y_position = c(1.3),
 #           xmin = c(1),
  #          xmax = c(4),
   #         annotation = "0.0091", 
    #        tip_length = 0.04)

# Load data 1
data$age <- as.factor(data$age)

# order levels
data$age <- ordered(data$age,
                        levels = c("young", "middle", "old"))
# xmin / xmax positions should match the x-axis labels' positions
ggplot(data, aes(x = passage, y = Dhat)) + 
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  # theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  # geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  # add zero level
  
  labs(title = "F-actin local fractal dimension in WJ1",
       caption = "") +
  
  geom_smooth(method='lm', formula= y ~ x) + # add linear regression
  
  scale_x_continuous(name = "Cell passage",
                     breaks = c(seq(7, 36, 1)), 
                     limits = c(7, 36))  +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   size = 9)) +
  #  
  scale_y_continuous(name = "Dlc", 
                     breaks = c(seq(1.5, 2, .1)), 
                     limits = c(1.5, 2))

fit <- lm(Dhat ~ passage, data = data)

summary(fit)

###### divide in stages
ggplot(data, aes(x = passage, y = Dhat,
                 color = age)) + 
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +

 # theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  # geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  # add zero level
  
  labs(title = "F-actin in MSCWJ1",
       subtitle = "Replicative senescense stages",
       caption = "") +
  
  geom_smooth(method='lm', formula = y ~ x) + # add linear regression
  
 scale_x_continuous(name = "Cell passage",
                   breaks = c(7,9, 12, 15,
                              18, 21, 25,
                              28, 36), 
                   limits = c(7, 36))  +
  
 theme(axis.text.x = element_text(angle = 0, hjust = 1,
                                  size = 10)) +
#  
 scale_y_continuous(name = "Fractal dimension", 
                  breaks = c(seq(1.5, 2, .1)), 
                    limits = c(1.5, 2)) +
  
  scale_colour_manual(values = c('darkgreen','blue','red'))

# LM

fit_young <- lm(Dhat ~ passage, data = data[data$age == 'young',])

summary(fit_young)

fit_middle<- lm(Dhat ~ passage, data = data[data$age == 'middle',])

summary(fit_middle)

fit_old<- lm(Dhat ~ passage, data = data[data$age == 'old',])

summary(fit_old)

#########################

# Compute linear regression
data1 <- data[data$passage == c("7", "9", "12", "15"),]

# Set X - axis names
CellScilabs <- c("7", "9", "12", "15", "18",
                 "21", "25", "27", "28", "35",  "36")


# Bar plot with signifiers 

df.summary <- group_by(data, passage) %>%
  summarise(
    sd = sd(Dhat, na.rm = TRUE),
    Dhat = mean(Dhat)
  )

df.summary

## 

ggplot(df.summary, aes(data$passage, data$Dhat)) +
  geom_bar(stat = "identity", fill = 'gray', 
           color = "black", size= 1, show.legend=TRUE) +
  geom_errorbar(aes(ymin = Dhat-sd, ymax = Dhat+sd), width = 0.2, size=1) +
  theme(
    # Change axis lines
    axis.line = element_line(size = 1),
    # Change axis ticks text labels: font color, size and face
    axis.text.x = element_text(face = "bold",
                               size = 12, angle = 90),     # Change x axis tick labels only
    axis.text.y = element_text(face = "bold", 
                               size = 12, angle = 0),     # Change y axis tick labels only
    # Change axis ticks line: font color, size, linetype and length
    axis.ticks = element_line(),      # Change ticks line fo all axes
    axis.ticks.x = element_line(),    # Change x axis ticks only
    axis.ticks.y = element_line(),    # Change y axis ticks only
    axis.ticks.length = unit(3, "pt") # Change the length of tick marks
  ) +
  geom_point() +
  ylim(0, 2) + 
  ggtitle("F-actin local fractal dimension") + 
  labs(y="Dlc", x = "Passage number")
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.95),
              xmin = c(5),
              xmax = c(9),
              annotation = "***", 
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.85),
              xmin = c(2),
              xmax = c(4),
              annotation = "*", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.9),
              xmin = c(9),
              xmax = c(11),
              annotation = "*", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.77),
              xmin = c(4),
              xmax = c(6),
              annotation = "**", 
              tip_length = 0.04)


# where plots saved:

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)

# Now, you can copy these files to your desired directory, as follows:
# file.copy(from=plots.png.paths, to="path_to_folder")
savehistory(file='myscript.R')
#####################################




compare_means(bTau ~ probe,  data = data1, method = "t.test")
compare_means(bTau ~ probe,  data = data1, method = "wilcox.test")

write.csv(compare_means(bTau ~ probe,  data = data1, method = "t.test"), file="t.test.bTau.csv")
write.csv(compare_means(bTau ~ probe,  data = data1, method = "wilcox.test"), file="wilcox.test.bTau.csv")

#
t.test <- compare_means(bTau ~ probe,  data = data1, method = "t.test", ref.group = 'p07')
wilcox.test <- compare_means(bTau ~ probe,  data = data1, method = "wilcox.test", ref.group = 'p07')


# Test for normality


kruskal.test(data$max_speed ~ data$probe)


write.xlsx(compare_means(max_speed ~ probe,  data = data, method = "kruskal.test"), 
           file = 'kruskal.test.max_speed.xlsx')

compare_means(bTau ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

compare_means(Rs ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

compare_means(Rval ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

write.xlsx(compare_means(max_speed ~ probe,  data = data, method = "wilcox.test"),
           file = 'wilcox.test.max_speed.xlsx')


#--------RhoA-----------

# Load data

data2 <- read.csv('colocRhoAF.csv')
data3 <- read.csv('colocRhoAH.csv')
head(data3)
data3$probe

# Name dependent variables
data2$probe <- as.factor(data2$probe)

# order levels
data2$probe <- ordered(data2$probe,
                       levels = c( "p09", "p15", 
                                   "p28", "p36"))
data3$probe <- ordered(data3$probe,
                       levels = c( "p09", "p15", 
                                   "p28", "p36"))
data3[data2$probe=='p15',]
#
# plot(data)
# png()
ggpairs(data3)
ggcorr(data3, palette = "RdBu", label = TRUE)
# dev.off()


#--aA4 Rvak---

qplot(Probe, Rval, data = data4,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = Probe,
      main = "Pearson's R value" ) + 
  labs(y = 'Rval',
       x = "Cell passage") +   theme_classic(base_size=14)+
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.2),
              xmin = c(1),
              xmax = c(2),
              annotation = "0.02", 
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.1),
              xmin = c(4),
              xmax = c(5),
              annotation = "0.0195", 
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1),
              xmin = c(2),
              xmax = c(3),
              annotation = "0.0045", 
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(.9),
              xmin = c(3),
              xmax = c(4),
              annotation = "0.004", 
              tip_length = 0.04)

compare_means(Dhat ~ passage,  data = data, method = "wilcox.test") # pairwise comparisons
kruskal.test(data4$Rval ~ data4$Probe)


compare_means(tM1 ~ Probe,  data = data4, method = "wilcox.test") # pairwise comparisons
kruskal.test(data4$tM1 ~ data4$Probe)

scatter3d(x = data$bTau, y = data$Rval, z = data$tM1, 
          groups = data$probe,
          xlab=deparse(substitute(bTau)), 
          ylab=deparse(substitute(Rval)),
          zlab=deparse(substitute(tM1)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE,
          ellipsoid = TRUE,
          level=0.5, ellipsoid.alpha=0.1, id=FALSE,
          model.summary=FALSE)

# Load data 1
data$probe <- as.factor(data$probe)

# order levels
data$probe <- ordered(data$probe,
                      levels = c("p07", "p09", "p12", "p15", "p18",
                                 "p21", "p25", "p27", "p28", "p35", "p36"))

# Set X - axis names
CellScilabs <- c("7", "9", "12", "15", "18",
                 "21", "25", "27", "28", "35",  "36")

# One-way ANOVA
# Compute the analysis of variance

res.aov <- aov(Dhat ~ passage, data = data)

# Summary of the analysis

summary(res.aov)
TukeyHSD(res.aov)
par(mar = c(4.5, 8, 4.5, 4.5))
plot(TukeyHSD(res.aov), las = 1)


#######################


ggplot(data[data$age=='young',], aes(x = passage, y = Dhat)) + 
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  # geom_hline(yintercept=0, linetype="dashed", color = "red") + # add zero level
  
  labs(title = "F-actin local fractal dimension in WJ1",
       caption = "lm: 0.000637 R-squared: -0.05105 ") +
  
  geom_smooth(method='lm', formula= y ~ x) + # add linear regression
  
  scale_x_continuous(name = "Cell passage",
                     breaks = c(seq(7, 36, 1)), 
                     limits = c(7, 36))  +
  #  
  scale_y_continuous(name = "Dlc", 
                     breaks = c(seq(1.5, 2, .1)), 
                     limits = c(1.5, 2))



ggplot(data[data$age=='old',], aes(x = passage, y = Dhat)) + 
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  # geom_hline(yintercept=0, linetype="dashed", color = "red") + # add zero level
  
  labs(title = "F-actin local fractal dimension in WJ1",
       caption = "lm: 0.000637 R-squared: -0.05105 ") +
  
  geom_smooth(method='lm', formula= y ~ x) + # add linear regression
  
  scale_x_continuous(name = "Cell passage",
                     breaks = c(seq(7, 36, 1)), 
                     limits = c(7, 36))  +
  #  
  scale_y_continuous(name = "Dlc", 
                     breaks = c(seq(1.5, 2, .1)), 
                     limits = c(1.5, 2))
  


plotmeans(data$Dhat ~ data$age,
          xlab = c('Passage groups'),
          ylab = c('Mean F-actin D-value'),
          use.t = TRUE)
