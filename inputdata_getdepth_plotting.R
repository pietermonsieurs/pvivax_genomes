library(ggplot2)
library(hrbrthemes)


depth_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/depth_analysis/depth_stats.summary.csv'

depth_data = read.csv(depth_file, header=FALSE)
colnames(depth_data) = c('sample', 'median_cov', 'pos_above_cutoff', 'pos_total', 'percent_callable')
head(depth_data)


hist(depth_data$median_cov, breaks=seq(0,2000,5), xlim=c(0,100))
sum(depth_data$median_cov < 5)

p_cov = ggplot(depth_data, aes(x=median_cov)) +
  geom_histogram( binwidth=5.1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() + 
  coord_cartesian(xlim=c(0,100)) + 
  labs(y= "Number of genomes", x = "Median coverage", size=20)
  
  # theme_ipsum() +
  # theme(
  #  plot.title = element_text(size=15)
  #)



hist(depth_data$percent_callable, breaks=seq(0,100,2.5), xlim=c(0,100))
sum(depth_data$percent_callable < 50)

p_call = ggplot(depth_data, aes(x=percent_callable)) +
  geom_histogram( binwidth=5.1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() + 
  coord_cartesian(xlim=c(0,100)) + 
  labs(y= "Number of genomes", x = "Percentage callable", size=20)


p_cov + p_call
