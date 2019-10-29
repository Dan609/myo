# The proportion of cells with pronounced activity of beta-galactosidase 
# during MSCWJ-1 cell line cultivation.

library(ggplot2)
library(DescTools)

# p9_ci90 <- BinomCI(62, 1025, conf.level = 0.90)
# Set confidence level to 95%
p9_ci95 <- BinomCI(224, 3724, conf.level = 0.95)
p15_ci95 <- BinomCI(494, 2404, conf.level = 0.95)
# p20_ci95 <- BinomCI(302, 1149, conf.level = 0.95)
p28_ci95 <- BinomCI(554, 1260, conf.level = 0.95)
p36_ci95 <- BinomCI(556, 1009, conf.level = 0.95)
# ?ombining observations into matrix
b_gal <- rbind(p9_ci95, p15_ci95, p28_ci95, p36_ci95)

# Convert matrix to data frame
df <- data.frame(b_gal)
# Set names
probe <- c("p9", "p15", "p28", "p36")
df <- cbind(probe, df)
df$probe <- as.factor(df$probe)
# Order levels
df$probe <- ordered(df$probe,
                       levels = c("p9", "p15", "p28", "p36"))
# Write results
write.table(df, file = 'beta-gal.csv', sep = ',')

# Prot results
ggplot(df, aes(probe, est)) +
  geom_bar(stat = "identity", fill = 'gray', width = 0.7,
           color = "black", size= 1, show.legend=TRUE) +
  geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci), width = 0.2, size=1) +
  geom_point() +
  ylim(0, .7) + 
  labs(title="Beta-galactosidase activity in MSCWJ-1 cell line",
       y="The proportion of stained cells", x = "Passage number") + theme_classic(base_size=14)

  head(df)
  
  
  qplot(probe, est, data = df, alpha = I(0.3), fill = probe,geom =  "point" ,
        main = "RhoA and nucleus colocalization", ) + 
  labs(y = 'bTau',
       x = "Cell passage") + theme_bw() +
    geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci), width = 0.2, size=1)








