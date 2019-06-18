#Skript für PCA Expressionsdaten alle Zelllinien für Matrix targetexpression und drivermutations.expression
#bind targetexpression und drivermutations.expression
rbind(drivermutations.expression, targetexpression)
PCA <- rbind(drivermutations.expression, targetexpression)
#Ausführung PCA 
pca = prcomp(t(PCA), center = F, scale. = F)
show(pca)
#Summary pf pca
summary(pca)
#man sieht das 90% im PC1 ist, so gut wie alle Information können durch PC1 festgehalten werden 
#str() to have a look at our PCA objects 
str(pca)
#Plotting PCA
plot(pca, type = "l")
#nun zum plotten neues package runterladen 
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

ggbiplot(pca)
#pca Punkte als Zelllinien darstellen 
ggbiplot(pca, labels=rownames(pca$x))
#PCA für PC2vsPC3
ggbiplot(pca, choices = 2:3)