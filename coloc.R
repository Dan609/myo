# Statistical analysis of colocalozation coefficients, Dan Bobkov, 2019
# Firstly, calculate coef in ImajeJ :
# Kendall's Tau-b rank correlation value = bTau
# Spearman's rank correlation value = Rs
# Manders' coefficients = tM1 and tM2
# Pearson's R value (above threshold) = Rval

# If you use this script in your publications, please cite:
# Bobkov D, Polyanskaya A, Musorina A, Lomert E, Shabelnikov S, Poljanskaya G. 
# Replicative senescence in MSCWJ-1 human umbilical cord mesenchymal stem cells is marked by characteristic changes in motility, 
# cytoskeletal organization, and RhoA localization. Mol Biol Rep. 2020;:1–17. DOI: 10.1007/s11033-020-05476-6

# Import libraries
# Load libraries-------------
library(GGally)
library(trajr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(xlsx)
library(ggplot2)
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
library(xtable)
# Load data

data1 <- read.csv('colocM9.csv')

head(data1)
fit1 <- glm(probe ~ Rval + tM1 + tM2 + bTau + Rs, data1,
           family = "quasibinomial")

fit <- glm(probe ~ Rval + tM1 + tM2 + bTau + Rs, data1,
           family = "binomial")

anova(fit, test='Chisq')
anova(fit1, test='Chisq')

newobject1 <- anova(fit, test='Chisq')
newobject2 <- anova(fit, test='Chisq')

print(xtable(newobject1, type = "latex"), file = "filename1.tex")
print(xtable(newobject2, type = "latex"), file = "filename2.tex")

data <- data1

data$probe

# Name dependent variables

data1$probe <- as.factor(data1$probe)

# order levels
data1$probe <- ordered(data1$probe,
                       levels = c("p07", "p09", "p12",
                                  "p15", "p18",
                                  "p21", "p25", "p27",
                                  "p28", "p35", "p36"))

#
# plot(data)
png()
ggpairs(data)
ggcorr(data, palette = "RdBu", label = TRUE)
dev.off()

# Kendall's Tau-b rank correlation value = bTau
# Spearman's rank correlation value = Rs
# Manders' coefficients = tM1 and tM2
# Pearson's R value (above threshold) = Rval


# Kruskal Wallis Test One Way Anova by Ranks
kruskal.test(data$bTau ~ data$probe)
kruskal.test(data$Rs ~ data$probe)
kruskal.test(data$Rval ~ data$probe)
kruskal.test(data$tM1 ~ data$probe)
kruskal.test(data$tM2 ~ data$probe)
#------bTau----
# One-way ANOVA
# Compute the analysis of variance
res.aov <- aov(bTau ~ probe, data = data1)
summary(res.aov) # Summary of the analysis
TukeyHSD(res.aov)
par(mar = c(4.5, 8, 4.5, 4.5))
plot(TukeyHSD(res.aov), las = 1)

compare_means(bTau ~ probe,  data = data, method = "kruskal.test")
data$probe

shapiro.test(data[data$probe=='p07',]$bTau)
shapiro.test(data[data$probe=='p09',]$bTau)
shapiro.test(data[data$probe=='p12',]$bTau)# data is not normally distributed
ggqqplot(data[data$probe=='p12',]$bTau)
shapiro.test(data[data$probe=='p15',]$bTau)
shapiro.test(data[data$probe=='p18',]$bTau)
shapiro.test(data[data$probe=='p21',]$bTau)
shapiro.test(data[data$probe=='p25',]$bTau)
shapiro.test(data[data$probe=='p27',]$bTau)
shapiro.test(data[data$probe=='p28',]$bTau)
shapiro.test(data[data$probe=='p35',]$bTau)
shapiro.test(data[data$probe=='p36',]$bTau)

ggqqplot(data$bTau, main = 'bTau')

compare_means(bTau ~ probe,  data = data1, method = "wilcox.test") # pairwise comparisons
compare_means(bTau ~ probe,  data = data1, method = "t.test") # pairwise comparisons

qplot(probe, bTau, data = data1,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Kendall's Tau-b rank correlation value", ) +
  labs(y = 'bTau',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1),
              xmin = c(1),
              xmax = c(9),
              annotation = "p = 0.028",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.87),
              xmin = c(1),
              xmax = c(4),
              annotation = "p = 0.096",
              tip_length = 0.04)  + theme_classic(base_size=14)


#----Rs---
ggqqplot(data$Rs, main = 'bTau')
compare_means(Rs ~ probe,  data = data1, method = "wilcox.test") # pairwise comparisons
compare_means(Rs ~ probe,  data = data, method = "t.test") # pairwise comparisons

qplot(probe, Rs, data = data1,
        geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
        main = "Spearman's rank correlation value", ) +
    labs(y = 'Rs',
         x = "Cell passage") +
    # xmin / xmax positions should match the x-axis labels' positions
    geom_signif(y_position = c(1.1),
                xmin = c(1),
                xmax = c(9),
                annotation = "p = 0.028",
                tip_length = 0.04)  + theme_classic(base_size=14)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.99),
              xmin = c(1),
              xmax = c(4),
              annotation = "p = 0.084",
              tip_length = 0.04)
#--


#-----Pearson's R value----
  ggqqplot(data$Rval, main = 'Rval')
  compare_means(Rval ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
  compare_means(Rval ~ probe,  data = data, method = "t.test") # pairwise comparisons

qplot(probe, Rval, data = data,
        geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
        main = "Pearson's R value") +
    labs(y = 'Rval',
         x = "Cell passage") +
    # xmin / xmax positions should match the x-axis labels' positions
    geom_signif(y_position = c(1),
                xmin = c(1),
                xmax = c(4),
                annotation = "p = 0.037",
                tip_length = 0.04) + theme_classic(base_size=14)
#-----

  #----
compare_means(Rval ~ probe,  data = data1, method = "wilcox.test") # pairwise comparisons
compare_means(Rval ~ probe,  data = data1, method = "t.test") # pairwise comparisons

  qplot(probe, Rval, data = data,
        geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
        main = "F-actin and Myosin-9 Colocalization in MSCWJ-1", ) +
    labs(y = 'Rval',
         x = "Cell passage") +
    # xmin / xmax positions should match the x-axis labels' positions
    geom_signif(y_position = c(1.1),
                xmin = c(1),
                xmax = c(4),
                annotation = "p = 0.026",
                tip_length = 0.04)   + theme_classic(base_size=14)  + # +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.85),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04)

#-----Manders' coefficient----



  shapiro.test(data[data$probe=='p07',]$tM1)
  shapiro.test(data[data$probe=='p09',]$tM1)
  shapiro.test(data[data$probe=='p12',]$tM1)# data is not normally distributed
  ggqqplot(data[data$probe=='p28',]$tM1)
  shapiro.test(data[data$probe=='p15',]$tM1)
  shapiro.test(data[data$probe=='p18',]$tM1)
  shapiro.test(data[data$probe=='p21',]$tM1)
  shapiro.test(data[data$probe=='p25',]$tM1)
  shapiro.test(data[data$probe=='p27',]$tM1)
  shapiro.test(data[data$probe=='p28',]$tM1)
  shapiro.test(data[data$probe=='p35',]$tM1)
  shapiro.test(data[data$probe=='p36',]$tM1)

  ggqqplot(data$tM2, main = 'tM')
  shapiro.test(data$tM1)
  shapiro.test(data$tM1)# data is not normally distributed
  ggqqplot(data[data$probe=='p15',]$tM1)
  compare_means(tM1 ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
  compare_means(tM2 ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
  compare_means(tM1 ~ probe,  data = data, method = "t.test") # pairwise comparisons

  qplot(probe, tM1, data = data1,
        geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe, #log = "y",
        main = "Manders' coefficient", ) +
    labs(y = 'tM',
         x = "Cell passage") +
    # xmin / xmax positions should match the x-axis labels' positions
    geom_signif(y_position = c(1.2),
                xmin = c(1),
                xmax = c(7),
                annotation = "p = 0.004",
                tip_length = 0.04)  +
    # xmin / xmax positions should match the x-axis labels' positions
    geom_signif(y_position = c(1.1),
                xmin = c(1),
                xmax = c(11),
                annotation = "p = 0.02",
                tip_length = 0.04)  + theme_classic(base_size=14)

  compare_means(tM1 ~ probe,  data = data1, method = "wilcox.test") # pairwise comparisons
  compare_means(tM1 ~ probe,  data = data1, method = "t.test") # pairwise comparisons
#----------

  head(data)

# Bar plot with signifiers

df.summary <- group_by(data1, probe) %>%
  summarise(
    sd = sd(bTau, na.rm = TRUE),
    bTau = mean(bTau)
  )

df.summary

##

ggplot(df.summary, aes(probe, bTau)) +
  geom_bar(stat = "identity", fill = 'gray',
           color = "black", size= 1, show.legend=TRUE) +
  geom_errorbar(aes(ymin = bTau-sd, ymax = bTau+sd), width = 0.2, size=1) +
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
  ylim(0, 1) +
  ggtitle("Colocalization of myosin-9 and F-actin in MSCWJ-1 cells") +
  labs(y="Kendall's Tau-b rank correlation value", x = "Passage number") +
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


##
# (1) Compute summary statistics for the variable probe
# (2) Bar plots of means + individual jitter points + errors
# Kendall's Tau-b rank correlation value

df.summary.bTau <- group_by(data1, probe) %>%
  summarise(
    sd = sd(bTau, na.rm = TRUE),
    bTau = mean(bTau)
  )

df.summary.bTau

df.bTau <- data1

ggplot(df.bTau, aes(probe, bTau)) +
  geom_bar(stat = "identity", data = df.summary.bTau,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2),
              color = "black") +
  geom_errorbar(
    aes(ymin = bTau-sd, ymax = bTau+sd),
    data = df.summary.bTau, width = 0.2)

# plotmeans

plotmeans(bTau ~ probe, data = data1, frame = FALSE, ylim = c(0, 1),
          mean.labels=FALSE, connect=TRUE, ccol = 'red',
          n.label=TRUE, text.n.label="n = ",
          xlab = "Passages", ylab = "Kendall's Tau-b rank correlation value",
          main="Colocalization of Myosin-9 and F-actin in WJMSC-1 cells,
          \nMean Plot with 95% CI") + scale_x_discrete(name ="Passages",
                                                       limits=c("p07", "p09", "p12",
                                                                "p15", "p18",
                                                                "p21", "p25", "p27",
                                                                "p28", "p35", "p36")) +
  scale_y_continuous(name="Kendall's Tau-b rank correlation value", limits=c(0, 1))






# where plots saved:

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)

# Now, you can copy these files to your desired directory, as follows:
# file.copy(from=plots.png.paths, to="path_to_folder")
savehistory(file='myscript.R')
#####################################

#

#
#
# boxplot(bTau ~ probe, data1)
# boxplot(Rval ~ probe, data1)
# boxplot(Rs ~ probe, data1)
# boxplot(tM1 ~ probe, data1)
# boxplot(tM2 ~ probe, data1)
#

#------------

#####################################
# Kendall's Tau-b rank correlation value
#####################################

# normality test
# perform the Shapiro-Wilk test of normality for one variable (univariate):
shapiro.test(data1$bTau)

# Kruskal Wallis Test One Way Anova by Ranks
kruskal.test(data1$bTau ~ data1$probe)

summary(data[data$probe=='p09',])
summary(data[data$probe=='p15',])
summary(data[data$probe=='p36',])
# Density plot: the density plot provides a visual judgment about whether the distribution is bell shaped.
ggdensity(data1$bTau,
          main = "Density plot of bTau",
          xlab = "bTau")


# Q-Q plot: Q-Q plot (or quantile-quantile plot) draws the correlation between a given sample and the normal distribution.
# qqPlot(data1$bTau)
ggqqplot(data1$bTau)


# Perform pairwise comparisons
compare_means(bTau ~ probe,  data = data1, method = "anova")
compare_means(bTau ~ probe,  data = data1, method = "kruskal.test")


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
# Kendall's Tau-b rank correlation value = bTau
# Spearman's rank correlation value = Rs
# Manders' coefficients = tM1 and tM2
# Pearson's R value (above threshold) = Rval
# Kruskal Wallis Test One Way Anova by Ranks
kruskal.test(data3$bTau ~ data3$probe)

kruskal.test(data2$Rs ~ data2$probe)

kruskal.test(data2$Rval ~ data2$probe)

#------bTau----
# One-way ANOVA
# Compute the analysis of variance
res.aov <- aov(bTau ~ probe, data = data3)
summary(res.aov) # Summary of the analysis
TukeyHSD(res.aov)
par(mar = c(4.5, 8, 4.5, 4.5))
plot(TukeyHSD(res.aov), las = 1)

compare_means(bTau ~ probe,  data = data3, method = "kruskal.test")
data3$probe

shapiro.test(data3$bTau)

shapiro.test(data[data$probe=='p07',]$bTau)
shapiro.test(data2[data2$probe=='p09',]$bTau)
shapiro.test(data[data$probe=='p12',]$bTau)# data is not normally distributed
ggqqplot(data3[data3$probe=='p15',]$bTau)
shapiro.test(data[data$probe=='p15',]$bTau)
shapiro.test(data[data$probe=='p18',]$bTau)
shapiro.test(data[data$probe=='p21',]$bTau)
shapiro.test(data[data$probe=='p25',]$bTau)
shapiro.test(data[data$probe=='p27',]$bTau)
shapiro.test(data[data$probe=='p28',]$bTau)
shapiro.test(data[data$probe=='p35',]$bTau)
shapiro.test(data[data$probe=='p36',]$bTau)

ggqqplot(data3$bTau, main = 'bTau')

compare_means(bTau ~ probe,  data = data3, method = "wilcox.test") # pairwise comparisons
compare_means(bTau ~ probe,  data = data3, method = "t.test") # pairwise comparisons

kruskal.test(data2$Rs ~ data2$probe)
kruskal.test(data3$Rs ~ data3$probe)

kruskal.test(data2$Rs ~ data2$probe)
kruskal.test(data3$Rs ~ data3$probe)

kruskal.test(data2$Rs ~ data2$probe)
kruskal.test(data3$Rs ~ data3$probe)

kruskal.test(data2$bTau ~ data2$probe)

compare_means(bTau ~ probe,  data = data2, method = "wilcox.test") # pairwise comparisons

kruskal.test(data3$bTau ~ data3$probe)

compare_means(bTau ~ probe,  data = data3, method = "wilcox.test") # pairwise comparisons

qplot(probe, bTau, data = data3,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "RhoA and Hoechst colocalization", ) +
  labs(y = 'bTau',
       x = "Cell passage") + theme_classic(base_size=14) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(.9),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(.8),
              xmin = c(3),
              xmax = c(4),
              annotation = "*",
              tip_length = 0.04)

#--Pearson RHoa

head(data2)


kruskal.test(data3$Rval ~ data3$probe)
kruskal.test(data2$Rval ~ data2$probe)
kruskal.test(data2$bTau ~ data2$probe)

compare_means(Rval ~ probe,  data = data3, method = "wilcox.test") # pairwise comparisons
compare_means(Rval ~ probe,  data = data2, method = "wilcox.test") # pairwise comparisons
compare_means(bTau ~ probe,  data = data2, method = "wilcox.test") # pairwise comparisons

qplot(probe, Rval, data = data2,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "RhoA and F-actin colocalization", ) +
  labs(y = 'Rval',
       x = "Cell passage") + theme_classic(base_size=14) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1),
              xmin = c(1),
              xmax = c(2),
              annotation = "***",
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(.9),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(.7),
              xmin = c(2),
              xmax = c(4),
              annotation = "**",
              tip_length = 0.04)




#--------alpha-actinin-4-----------

# Load data

data4 <- read.csv('colocA4.csv', sep = ';')
data3 <- read.csv('colocRhoAH.csv')
data4 <- read.csv('alpha.csv', sep = ';')
head(data4)

data4$Probe

# Name dependent variables
data4$Probe <- as.factor(data4$Probe)

# order levels
data4$Probe <- ordered(data4$Probe,
                       levels = c( "p07","p09","p12", "p15",
                                   "p21"))
data4[data4$Probe=='p15',]
#
# plot(data)
# png()
ggpairs(data4)
ggcorr(data4, palette = "RdBu", label = TRUE)
# dev.off()
# Kendall's Tau-b rank correlation value = bTau
# Spearman's rank correlation value = Rs
# Manders' coefficients = tM1 and tM2
# Pearson's R value (above threshold) = Rval
# Kruskal Wallis Test One Way Anova by Ranks
kruskal.test(data4$bTau ~ data4$Probe)

kruskal.test(data2$Rs ~ data2$probe)

kruskal.test(data2$Rval ~ data2$probe)

#------bTau----
# One-way ANOVA
# Compute the analysis of variance
res.aov <- aov(bTau ~ probe, data = data3)
summary(res.aov) # Summary of the analysis
TukeyHSD(res.aov)
par(mar = c(4.5, 8, 4.5, 4.5))
plot(TukeyHSD(res.aov), las = 1)

compare_means(bTau ~ probe,  data = data3, method = "kruskal.test")
data3$probe

shapiro.test(data3$bTau)

shapiro.test(data[data$probe=='p07',]$bTau)
shapiro.test(data2[data2$probe=='p09',]$bTau)
shapiro.test(data[data$probe=='p12',]$bTau)# data is not normally distributed
ggqqplot(data3[data3$probe=='p15',]$bTau)
shapiro.test(data[data$probe=='p15',]$bTau)
shapiro.test(data[data$probe=='p18',]$bTau)
shapiro.test(data[data$probe=='p21',]$bTau)
shapiro.test(data[data$probe=='p25',]$bTau)
shapiro.test(data[data$probe=='p27',]$bTau)
shapiro.test(data[data$probe=='p28',]$bTau)
shapiro.test(data[data$probe=='p35',]$bTau)
shapiro.test(data[data$probe=='p36',]$bTau)

ggqqplot(data3$bTau, main = 'bTau')

compare_means(bTau ~ Probe,  data = data4, method = "wilcox.test") # pairwise comparisons
compare_means(bTau ~ probe,  data = data3, method = "t.test") # pairwise comparisons
compare_means(Rs ~ Probe,  data = data4, method = "wilcox.test") # pairwise comparisons
compare_means(bTau ~ Probe,  data = data4, method = "wilcox.test") # pairwise comparisons
compare_means(bTau ~ Probe,  data = data4, method = "wilcox.test") # pairwise comparisons

qplot(Probe, bTau, data = data4,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = Probe,
      main = "Kendall's Tau-b rank correlation value" ) +
  labs(y = 'bTau',
       x = "Cell passage") +   theme_classic(base_size=14)+
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1),
              xmin = c(1),
              xmax = c(4),
              annotation = "0.0011",
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.9),
              xmin = c(2),
              xmax = c(4),
              annotation = "0.0424",
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(.8),
              xmin = c(3),
              xmax = c(4),
              annotation = "0.0107",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(.7),
              xmin = c(4),
              xmax = c(5),
              annotation = "0.0159",
              tip_length = 0.04)

#--aA4 Rs---

qplot(Probe, Rs, data = data4,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = Probe,
      main = "Spearman's rank correlation value", ) +
  labs(y = 'Rs',
       x = "Cell passage") +   theme_classic(base_size=14)+
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.2),
              xmin = c(1),
              xmax = c(4),
              annotation = "0.0022",
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.1),
              xmin = c(2),
              xmax = c(4),
              annotation = "0.0424",
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1),
              xmin = c(3),
              xmax = c(4),
              annotation = "0.0081",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(.9),
              xmin = c(4),
              xmax = c(5),
              annotation = "0.0159",
              tip_length = 0.04)

compare_means(Rs ~ Probe,  data = data4, method = "wilcox.test") # pairwise comparisons
kruskal.test(data4$Rs ~ data4$Probe)

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

compare_means(Rval ~ Probe,  data = data4, method = "wilcox.test") # pairwise comparisons
kruskal.test(data4$Rval ~ data4$Probe)

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

compare_means(Rval ~ Probe,  data = data4, method = "wilcox.test") # pairwise comparisons
kruskal.test(data4$Rval ~ data4$Probe)

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

compare_means(Rval ~ Probe,  data = data4, method = "wilcox.test") # pairwise comparisons
kruskal.test(data4$Rval ~ data4$Probe)



#--aA4 tM1---

qplot(Probe, tM1, data = data4,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = Probe,
      main = "Manders' coefficient" ) +
  labs(y = 'tM1',
       x = "Cell passage") +   theme_classic(base_size=14)+
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.3),
              xmin = c(1),
              xmax = c(4),
              annotation = "0.0091",
              tip_length = 0.04)  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.2),
              xmin = c(3),
              xmax = c(4),
              annotation = "0.0081",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.1),
              xmin = c(4),
              xmax = c(5),
              annotation = "0.0159",
              tip_length = 0.04)

compare_means(tM1 ~ Probe,  data = data4, method = "wilcox.test") # pairwise comparisons
kruskal.test(data4$tM1 ~ data4$Probe)
