# Load necessary libraries
library(ggplot2)

#### Elements to Sr Relationships ####
element_to_Sr = read.csv('Element_to_Sr.csv')

element_to_Sr_BH = element_to_Sr[element_to_Sr$Drainage %in% "Big Hole", ]

element_to_Sr_BH_few = element_to_Sr_BH[-c(11:13, 16:37),]

element_to_Sr_BV = element_to_Sr[element_to_Sr$Drainage %in% "Beaverhead", ]

element_to_Sr_RB = element_to_Sr[element_to_Sr$Drainage %in% "Ruby", ]

# Calcium to Sr
ggplot(data = element_to_Sr_BH_few, aes(x = X87Sr.86Sr, y = Mg.Ca, color = Tributary)) +
  geom_point(size = 8) +
  facet_wrap(~Drainage) 

# overall strontium differentiation
ggplot(data = element_to_Sr, aes(x = Tributary, y = X87Sr.86Sr)) +
  geom_point() +
  facet_wrap(~Drainage)


#### PCA Analysis ####
elemental_data <- read.csv('Element_to_Sr.csv')
elemental_data = elemental_data[-6, ] # removing canyon creek
elemental_ratios = elemental_data[, -c(1:8, 14:15)] # extracting just the ratios we want

# Standardize the data
standardized_data <- scale(elemental_ratios)  

# Perform PCA
pca_result <- prcomp(standardized_data, scale = TRUE)

# Create a data frame with the principal components
pc_df <- data.frame(pca_result$x[, 1:2], Tributary = elemental_data$Tributary, Drainage = elemental_data$Drainage)

# Biplot
ggplot(pc_df, aes(x = PC1, y = PC2, color = Drainage)) +
  geom_point()
  