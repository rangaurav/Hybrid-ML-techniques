###############################################################
# Credit: Gaurav Sharma                            #
# Email : ergauravrandev@gmail.com                                    #
###############################################################


#Include Library
install.packages("cluster")
install.packages("devtools")
install.packages("NbClust")
devtools::install_github("kassambara/factoextra")
library(cluster)
library(factoextra)
library(NbClust)
df = read.csv("dataset.csv",TRUE,",")    #Load data from file
ncol(df)
sdf <- scale(df)
head(sdf, n = 3)

kmeans(sdf, 10, iter.max = 10, nstart = 1)    #Run K-mean for cluster size ranging 1-10 iteratively 


fviz_nbclust(df, kmeans, method = "wss")  + 
  geom_vline(xintercept = 4, linetype = 2)     #Plot WSSE garph for 10 iterations


set.seed(123)

km.res <- kmeans(df, 4, nstart = 25)   #Run K-means for cluster size 4         

print(km.res)

dd <- cbind(df, cluster = km.res$cluster)

modelName <- "Kmeans"

write.csv(dd, file=paste(modelName,"4.csv",sep=''), row.names=FALSE)
