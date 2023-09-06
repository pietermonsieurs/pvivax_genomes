library(ggplot2)

k_cv_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/admixture.Kvalues.CrossValidation.csv'
k_cv_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/admixture.Kvalues.CrossValidation.SAM.csv'
k_data = read.csv(k_cv_file, 
                  sep=c(" "),
                  header=FALSE)
k_data = k_data[,3:4]
colnames(k_data) = c('K', 'CV')
k_data$K = gsub("\\(K=|\\)|:", "", k_data$K)
k_data$K = as.integer(k_data$K)
k_data

ggplot(k_data, aes(x=K, y=CV)) + 
  geom_point() + 
  geom_line() + 
  theme_bw()
