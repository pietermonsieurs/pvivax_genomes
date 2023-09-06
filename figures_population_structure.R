library(ggplot2)
library(cowplot)
library(ggpubr)
library(gridExtra)

picture_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/pictures_pub/'

#### 1. SAM ####

## import the picture from the script plink_admixture_plotting.R: this can 
## be done via running directly in the script, or loading the object
admixture_file = paste0(picture_dir, 'admixture_SAM.rds')
p_admix = readRDS(admixture_file)
p_admix

pca_file = paste0(picture_dir, 'pca_1_2_SAM.rds')
p_pca = readRDS(pca_file)
p_pca


# ggtree_file 
tree_file = paste0(picture_dir, 'tree_sam.rds')
p_tree = readRDS(tree_file)
p_tree


# Combine plots in two rows
grid.arrange(p_pca, p_tree, p_admix, ncol=1, nrow=3, top="A, B",
             widths=c(3), heights=c(3,3,3))
p_admix


p = plot_grid(plot_grid(p_pca, p_tree, ncol = 2, labels=c("A", "B")), plot_grid(p_admix, ncol = 1, labels="C"), nrow=2, label_size = 16)
output_file = paste0(output_dir, 'pop_structure_sam.png')
ggsave(output_file, plot=p, width=12, height = 8, dpi=300)





#### 2. world ####


## import the picture from the script plink_admixture_plotting.R: this can 
## be done via running directly in the script, or loading the object
admixture_file = paste0(picture_dir, 'admixture_world.rds')
p_admix = readRDS(admixture_file)
p_admix

pca_file = paste0(picture_dir, 'pca_1_2_world.rds')
p_pca = readRDS(pca_file)
p_pca


# ggtree_file 
tree_file = paste0(picture_dir, 'tree_world.rds')
p_tree = readRDS(tree_file)
p_tree


# Combine plots in two rows
# grid.arrange(p_pca, p_tree, p_admix, ncol=1, nrow=3, top="A, B",
#              widths=c(3), heights=c(3,3,3))


p = plot_grid(plot_grid(p_pca, p_tree, ncol = 2, labels=c("A", "B")), plot_grid(p_admix, ncol = 1, labels="C"), nrow=2, label_size = 16)
p
output_file = paste0(output_dir, 'pop_structure_world.png')
ggsave(output_file, plot=p, width=12, height = 8, dpi=300)






library(gridExtra)

plot1 <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
plot2 <- ggplot(mtcars, aes(x = wt, y = hp)) + geom_point()
plot3 <- ggplot(mtcars, aes(x = wt, y = drat)) + geom_point()

grid.arrange(plot1, plot2, plot3, nrow = 3, ncol = 1, heights = c(2, 2, 1))
plot_grid(plot_grid(plot1, plot2, ncol = 2), plot_grid(plot3, ncol = 2), labels = c("A", "B", "C"), label_size = 16)
