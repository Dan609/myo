# The proportion of cells with pronounced activity of beta-galactosidase 
# during MSCWJ-1 cell line cultivation.

install.packages("DescTools")
library(DescTools)

# p9_ci90 <- BinomCI(62, 1025, conf.level = 0.90)
# Set confidence level to 95%
p9_ci95 <- BinomCI(224, 3724, conf.level = 0.95)
p15_ci95 <- BinomCI(494, 2404, conf.level = 0.95)
p20_ci95 <- BinomCI(302, 1149, conf.level = 0.95)
p28_ci95 <- BinomCI(554, 1260, conf.level = 0.95)

# Ñombining observations into matrix
b_gal <- rbind(p9_ci95, p15_ci95, p20_ci95, p28_ci95)

# Convert matrix to data frame
df <- data.frame(b_gal)
# Set names
probe <- c("p9", "p15", "p20", "p28")
df <- cbind(probe, df)
df$probe <- as.factor(df$probe)
# Order levels
df$probe <- ordered(df$probe,
                       levels = c("p9", "p15", "p20", "p28"))
# Write results
write.table(df, file = 'beta-gal.csv', sep = ',')

# Prot results
library(ggplot2)

ggplot(df, aes(probe, est)) +
  geom_bar(stat = "identity", fill = 'gray', width = 0.7,
           color = "black", size= 1, show.legend=TRUE) +
  geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci), width = 0.2, size=1) +
  geom_point() +
  ylim(0, 0.5) + 
  labs(title="Beta-galactosidase activity in MSCWJ-1 cell line",
       y="The proportion of stained cells", x = "Passage number") + theme_classic(base_size=14)# +theme(
# Change axis lines
axis.line = element_line(size = 1),
# Change axis ticks text labels: font color, size and face
axis.text.x = element_text(face = "bold",
                           size = 13, angle = 0),     # Change x axis tick labels only
axis.text.y = element_text(face = "bold", 
                           size = 13, angle = 0),     # Change y axis tick labels only
# Change axis ticks line: font color, size, linetype and length
axis.ticks = element_line(size = 1),      # Change ticks line fo all axes
axis.ticks.x = element_line(),    # Change x axis ticks only
axis.ticks.y = element_line(),    # Change y axis ticks only
axis.ticks.length = unit(5, "pt"), # Change the length of tick marks
axis.title.x = element_text(size=14, face="bold"),
axis.title.y = element_text(size=14, face="bold"),
plot.title = element_text(size=18, face="bold")
) +
  # xmin / xmax positions should match the x-axis labels' positions
  # geom_signif(y_position = c(0.85), size = 1, textsize = 12,
  #            xmin = c(1),
  #            xmax = c(2),
  #            annotation = "*", 
  #            tip_length = 0.04)

  head(df)
  
  
  qplot(probe, est, data = df, alpha = I(0.3), fill = probe,geom =  "point" ,
        main = "RhoA and nucleus colocalization", ) + 
  labs(y = 'bTau',
       x = "Cell passage") + theme_bw() +
    geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci), width = 0.2, size=1)








