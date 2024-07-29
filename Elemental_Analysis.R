# Load necessary libraries
library(ggplot2)

#### Elements to Sr Relationships ####
element_to_Sr = read.csv('Element_to_Sr.csv')

# Calcium to Sr
ggplot(data = element_to_Sr, aes(x = X87Sr.86Sr, y = Ca_mgL, color = Tributary)) +
  geom_point() +
  facet_wrap(~Drainage) 

# Barium to Sr
ggplot(data = element_to_Sr, aes(x = X87Sr.86Sr, y = Ba_mgL, color = Tributary)) +
  geom_point() +
  facet_wrap(~Drainage)

# Magnesium to Sr
ggplot(data = element_to_Sr, aes(x = X87Sr.86Sr, y = Mg_mgL, color = Tributary)) +
  geom_point() +
  facet_wrap(~Drainage)

# Manganses to Sr
ggplot(data = element_to_Sr, aes(x = X87Sr.86Sr, y = Mn_mgL, color = Tributary)) +
  geom_point() +
  facet_wrap(~Drainage)

# Strontium element to ratio
ggplot(data = element_to_Sr, aes(x = X87Sr.86Sr, y = Sr_mgL, color = Tributary)) +
  geom_point() +
  facet_wrap(~Drainage)


#### PCA Analysis ####
elemental_data <- read.csv('Elemental_Ratios.csv')
elemental_ratios = elemental_data[, -c(1:11)] # extracting just the ratios we want

# Standardize the data
standardized_data <- scale(elemental_ratios)  

# Perform PCA
pca_result <- prcomp(standardized_data, scale = TRUE)

# Create a data frame with the principal components
pc_df <- data.frame(pca_result$x[, 1:2], Tributary = elemental_data$Tributary, Drainage = elemental_data$Drainage)

# Biplot
ggplot(pc_df, aes(x = PC1, y = PC2, color = Drainage)) +
  geom_point()
  