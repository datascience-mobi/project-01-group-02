#Versuch Hierarchical clustering zu machen mit targetexpression
#Using the most variabl
topVar = apply(targetexpression, 1, var)
summary(topVar)
targetexpression.topVar = targetexpression[topVar > quantile(topVar, probs = 0.75), ]
dim(targetexpression.topVar)
rm(topVar)
#Creating a correlation based distance matrix
cor.mat = cor(targetexpression.topVar, method = "spearman")
cor.dist = as.dist(1 - cor.mat)
cor.hc = hclust(cor.dist, method = "ward.D2")
cor.hc = as.dendrogram(cor.hc)
#install packages dendextend für dendogramm
install.packages('dendextend')
library(dendextend)
#Plotting the dendrogram
plot(cor.hc, las = 2, cex.lab = 0.7)
