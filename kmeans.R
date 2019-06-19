

# k means 
# mit welchen daten ist es sinnvoll einen kmeans zu machen? 
# wir vegleichen ja immer die Zelllinien untereinander also ein Punkt ist eine Zelllinen
# also villt die Zelllinien so umbenen dass man sieht welche Driver Mutations dort Mutiert sind
# wie machen wir das weil die Drivermutations tauchen ja nicht in der Mutationsmatrix auf also 
# villt die Zelllinien nach den Drivermutations benennnen die die stärkste Abweichung von 
# den anderen Zelllinien haben in der Drivermutation
# oder wir sagen wir benennen sie nach den Mutationen der Driver Mutations obwohl diese so selten sind 
# es gibt so viele verschiedene RAS deswegen ist es villt egal das sie so selten auftauchen
# Expressionsmatrix von taget genes mit den driver mutations
# von der variation der Gene ? also die Expressionsmatrix normalisieren und dann davon K-means 
# von copynumber daten 
# die Hoffnung bei diesen beiden ist das wir ein k-means rausbekommen bei denen die 
# mit den gleichen DriverMuatations in die gleichen Cluster kommen


# Zelllinien welche eine mutation in den Drivermutations haben 
namesdrivers <- c(rownames(BRAFexpression),rownames(NF1expression), rownames(RASexpression), rownames(WTexpression))
zuordnenZelllinienDM <- Mutation_1dataframe[which(Mutation_1dataframe$Hugo_Symbol %in% namesdrivers),]


km = kmeans(x =t(targetexpression), centers = 6, nstart = 30)
wss = sapply(1:11, function(k) {
  kmeans(x = t(targetexpression), centers = k)$tot.withinss
})
plot(1:11, wss, type = "b", pch = 19, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares")


D = dist(t(targetexpression))
library(cluster)
km = kmeans(x = t(targetexpression), centers = 8, nstart = 100)

fviz_cluster(km, data = t(targetexpression))

D = dist(t(targetexpression))
plot(silhouette(km$cluster, D))
mean(silhouette(km$cluster, D))
a = silhouette(km$cluster, D)
mean(a[,3])


# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(t(targetexpression), centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(t(targetexpression)))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")


library(tidyverse)  
library(cluster)    
library(factoextra)

fviz_nbclust(t(targetexpression), kmeans, method = "silhouette")

km2 <- kmeans(t(targetexpression), centers = 2, nstart = 25)
km3 <- kmeans(t(targetexpression), centers = 6, nstart = 25)
km4 <- kmeans(t(targetexpression), centers = 4, nstart = 25)
km5 <- kmeans(t(targetexpression), centers = 9, nstart = 25)

# plots to compare
p1 <- fviz_cluster(km2, geom = "text", labelsize = 8, data = t(targetexpression)) + ggtitle("k = 2")
p2 <- fviz_cluster(km3, geom = "text", labelsize = 8, data = t(targetexpression)) + ggtitle("k = 6")
p3 <- fviz_cluster(km4, geom = "text", labelsize = 8, data = t(targetexpression)) + ggtitle("k = 4")
p4 <- fviz_cluster(km5, geom = "text", labelsize = 8, data = t(targetexpression)) + ggtitle("k = 9")

library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)

