#######################################################################################################
# This script is created to implement unsupervised clustering algorithms on COVID-19 dataset
# Author: Tanya Jain/Kavya Anantha Rao
# Date Created: April 23th 2020
#######################################################################################################

#######################################################################################
# Install the required packages and load the libraries
#######################################################################################
rm(list=ls())
library(kohonen)
library(readxl)
library(janitor)
library(imputeTS)
library(ggplot2)
library(ggfortify)
library(cluster)
library(tidyverse)
library(fpc)

##################################################
#Load the data
##################################################
## choose full data
comp_data <- read_excel("original_dataset.xlsx") %>% clean_names() 
comp_data<-as.data.frame(comp_data)
colnames(comp_data)
view(comp_data)
## remove records having 90% of na values
na_check <-sapply(comp_data, function(y) sum(length(which(is.na(y)))))
na_check <- data.frame(na_check)
na_pct <- na_check$na_check < 5644 * 0.9 
new_data <- comp_data[,na_pct]

#keep numerics
new_data1 <- new_data[,1:20]
sum(rowSums(is.na(new_data1)) > 13)
table(new_data1$sars_cov_2_exam_result)
head(new_data1)

#retain only records that have all data
my_data <- new_data1[(rowSums(is.na(new_data1)) < 1),]
table(new_data2$sars_cov_2_exam_result)
table(rowSums(is.na(new_data2)))
my_data<-as.data.frame(my_data)
head(my_data)
my_datasum<-my_data

## Data Normalization
#Clean the data for columns giving complete blood count values
scaled_data <- scale(my_datasum[,c(7:20)])#only numeric columns
head(scaled_data)
dim(scaled_data)
colnames(scaled_data)

##################################################
# Data Statistics and Visualizations
##################################################
treat_names<-names(comp_data) 
rep_comp_data<-comp_data

covid_stats <- rep_comp_data %>% 
  mutate(
    sars_cov_2_exam_result = factor(rep_comp_data$sars_cov_2_exam_result, levels = c('positive','negative')),
    adm_type = ifelse(rep_comp_data$patient_addmited_to_regular_ward_1_yes_0_no == 1, 'regular_ward',
                      ifelse(rep_comp_data$patient_addmited_to_semi_intensive_unit_1_yes_0_no == 1, 'semi_intensive',
                             ifelse(rep_comp_data$patient_addmited_to_intensive_care_unit_1_yes_0_no == 1, 'intensive_care_unit','discharged'))),
    adm_type = factor(adm_type,levels = c('discharged','regular_ward','semi_intensive','intensive_care_unit'))
  )

covid_stats
test_result<- rep_comp_data %>% count(sars_cov_2_exam_result)

#display total number of cases
ggplot(test_result, aes(sars_cov_2_exam_result, n, fill = sars_cov_2_exam_result))+
  geom_col()+
  geom_label(aes(label = n))+
  theme_classic()+
  scale_fill_manual(values = c('green','red'))+
  labs(x = '', title = 'Cases tested to COVID-19 in Einstein',
       fill = 'COVID')

ggplot(rep_comp_data, aes(patient_age_quantile, fill = sars_cov_2_exam_result,alpha=0.5))+
  geom_density()+
  scale_fill_manual(values = c('green','red'))+
  theme_classic()+
  labs(title = 'Distribution of Age Quantles by COVID results')

###############################################################################################
### Self Organized Maps
###############################################################################################
set.seed(111)
getgrid_som <- somgrid(xdim = 5, ydim = 5, topo = "hexagonal")
cbc_som <- som(scaled_data, grid = getgrid_som, rlen = 1000)
summary(cbc_som)
code_val <- cbc_som$codes[[1]]

#Training Progress Plot
plot(cbc_som, type="changes", main = "Training Progress")
#Fan Diagram
plot(cbc_som, type="codes", main = "Fan Diagram")
#Distance Neighbor Plot
plot(cbc_som, type="dist.neighbours", main = "Average distance to neighbor neuron") 
#Node Count Plot
plot(cbc_som, type="counts", main="Node count")
#Mapping Plot
plot(cbc_som, type="mapping", main = "Mapping Plot")

## Determine if the suspect is positive or negative for COVID-19 
#negative/res=0
idx0<- as.numeric(as.factor(my_datasum$sars_cov_2_exam_result))-1
for(i in 1:598){
  if(idx0[i]==0){
    idx0[i]=0}
  else{idx0[i]=" "}
}

par(mfrow = c(1,2))
plot(cbc_som, type="mapping",labels =  idx0, col = 3,main = "COVID-19 Negative")

#positive/res=1
idx1<- as.numeric(as.factor(my_datasum$sars_cov_2_exam_result))-1
for(i in 1:598){
  if(idx1[i]==1){
    idx1[i]=1}
  else{idx1[i]=" "}
}
plot(cbc_som, type="mapping",labels =  idx1, col = 2,main = "COVID-19 Positive")

graphics.off()

## Determine if the patient was admitted to regular, semi-intensive or intensive ward
#regular ward
id_reg<- as.numeric(my_datasum$patient_addmited_to_regular_ward_1_yes_0_no)

for (x in 1:598){
  if (id_reg[x] ==1){
    id_reg[x]   
  }else{
    id_reg[x]=" "
  }
}
plot(cbc_som, type="mapping",labels =  id_reg, col = 6,main = "Patients admitted to Regular Ward")

#semi-intensive ward
id_semi<- as.numeric(my_datasum$patient_addmited_to_semi_intensive_unit_1_yes_0_no)
for (x in 1:598){
  if (id_semi[x] ==1){
    id_semi[x]   
  }else{
    id_semi[x]=" "
  }
}

plot(cbc_som, type="mapping",labels =  id_semi, col = 4,main = "Patients admitted to Semi-Intensive Ward")

#intensive ward
id_icu<- as.numeric(my_datasum$patient_addmited_to_intensive_care_unit_1_yes_0_no)
for (x in 1:598){
  if (id_icu[x] ==1){
    id_icu[x]   
  }else{
    id_icu[x]=" "
  }
}
plot(cbc_som, type="mapping",labels =  id_icu, col =2,main = "Patients admitted to Intensive Care Ward")

graphics.off()

####################################################################################################
# Principal Component Analysis
####################################################################################################
num_data <- as.data.frame(my_datasum[,c(2,4:17)]) # Select only numeric variables
row.names(num_data) <- my_datasum$patient_id
num_data$res_val<- as.numeric(as.factor(my_datasum$sars_cov_2_exam_result))-1
head(num_data)

#check data to scale
colMeans(num_data) # View column means
apply(num_data, 2, sd) # View column standard deviations

# Fit PCA
pca_out <- prcomp(scale(num_data[-16]), center = TRUE, scale = FALSE)
summary(pca_out) 
loadings<-pca_out$rotation
autoplot(pca_out, data = num_data,loadings = TRUE, loadings.colour = 'deeppink4',
         loadings.label = TRUE, loadings.label.size = 3,frame = TRUE, frame.type = 'norm')
biplot(pca_out) 

num_data$res_val <- as.factor(num_data$res_val)

# principal components
plot(pca_out,col="royalblue4")

# Check the Eigen values
eig_vals <- (pca_out$sdev)^2 

#show variation in Eigen values
pve <- eig_vals / sum(eig_vals)

## Plot variance explained for each principal component
plot(pve, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     ylim = c(0, 1), type = "b",col="royalblue4")

## cumulative proportion of variance
plot(cumsum(pve), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     ylim = c(0, 1), type = "b",col="deeppink4")

## Plot for Eigen Values
plot(eig_vals, xlab = 'Principal Components', ylab = 'Eigenvalue', main = 'Scree Plot')
lines(eig_vals, col="green",type="o")

autoplot(pca_out, data = num_data, col = 'res_val',
              loadings = TRUE, loadings.colour = 'blue',
              loadings.label = TRUE, loadings.label.size = 3,frame = TRUE, frame.type = 'norm')



####################################################################################################
# Hierarchical Clustering
####################################################################################################
num_scaleddata <- scale(num_data[,-16])
index <- sample(c(1:length(num_scaleddata[,1])))
scale_data_new <- num_scaleddata[index, ]
dist_val_scale <- dist(scale_data_new)
dim(as.matrix(dist_val_scale))
hclust_comp <- hclust(dist(dist_val_scale), method = "complete") 
summary(hclust_comp) 

# Plot dendogram
plot(hclust_comp) 
rect.hclust(hclust_comp,k=6, border="red")

dendogram_cut <- cutree(hclust_comp, k = 6) # Cut tree so that it has 4 clusters 
x11()
plot(dendogram_cut,col="maroon")
table(dendogram_cut,num_data$res_val) # Compare cluster membership to actual diagnoses
table(dendogram_cut)

#######################################################################################
# Use code_val from SOM for distance in hierarchical
#####################################################################################
dist_val_new<-dist(code_val1)

clust_new <- hclust(dist_val_new)
dendogram_cut_new<-cutree(clust_new,h=7) ## Clusters = 6
dendogram_cut_new

my_pal <- c("red", "royalblue4", "yellow","cyan","forestgreen","deeppink4")
my_bhcol <- my_pal[dendogram_cut_new]

graphics.off()

plot(cbc_som1, type = "mapping", col = "black", bgcol = my_bhcol)
add.cluster.boundaries(cbc_som1, dendogram_cut_new)

####################################################################################################
# Determine the number of clusters
####################################################################################################
plot_val <- c()
for (k in 1:15) {
  numout_val <- kmeans(scaled_data, k, nstart = 20)
  plot_val[k] <- rbind(numout_val$tot.withinss)
}
plot_val
plot(plot_val, type = "o", col = 'blue' ,xlab = "k",ylab='Total WSS',main="Elbow plot")#6

# gap statistics for kmeans
gapstat_kmeans <- clusGap(num_data[,-16], kmeans, nstart = 20, K.max = 15, B = 100)

####################################################################################################
# Implement K-means
####################################################################################################
kmeans_fit <- kmeans(scaled_data, centers = 6)

table(dendogram_cut,num_data$res_val) # Hierarchical to actual resultss
table(kmeans_fit$cluster, num_data$res_val) #k-means to actual results
table(kmeans_fit$cluster, dendogram_cut) #k-means wrt hierarchical clusters


## Plotting Clusters of k-means using data after PCA
PC1 <- pca_out$x[,1]
PC2 <- pca_out$x[,2] 
PCA_data <- cbind(PC1, PC2)
km_pca <- kmeans(PCA_data, centers = 6, nstart = 10)
plot(PCA_data, col = km_pca$cluster, main = "k-means with PCA")
points(km_pca$centers, col = 1:6, pch = 8, cex= 2)

## Plotting Clusters of k-means using processed data
km <- kmeans(scaled_data, centers =6, nstart = 10)
plot(num_data[ ,c("hemoglobin", "leukocytes")], col = km$cluster, main = "k-means without PCA")
points(km$centers, col = 1:6, pch = 8, cex= 2)

####################################################################################################
# Implement K-medoids
####################################################################################################
mediods.out<-pamk(scaled_data,5)
table(mediods.out$pamobject$clustering,num_data$res_val)
#layout(matrix(c(1,2),1,2)) 
plot(mediods.out$pamobject)
 