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
#PCA für drivermutations.knockdown
pca = prcomp(drivermutations.knockdown, center = F, scale. = F)
summary(pca)
str(pca)
plot(pca, type = "l")
library(ggbiplot)
ggbiplot(pca)
ggbiplot(pca, choices = 2:3)
ggbiplot(pca, labels=rownames(pca$x))
#PCA für drivermutations.knockdown.prob
pca = prcomp(drivermutations.knockdown.prob, center = F, scale. = F)
summary(pca)
str(pca)
plot(pca, type = "l")
library(ggbiplot)
ggbiplot(pca)
ggbiplot(pca, choices = 2:3)
#PCA für drivermutations.copynumber
pca = prcomp(drivermutations.copynumber, center = F, scale. = F)
summary(pca)
str(pca)
plot(pca, type = "l")
#fällt auf, dass der Knick etwas anders ist als die pca´s davor. Bedeutung ?
library(ggbiplot)
ggbiplot(pca)
#kann man gut erkennen
ggbiplot(pca, choices = 2:3)
#PCA von driver knockdown und prob.
rbind(drivermutations.knockdown.prob, drivermutations.knockdown)
PCA2 <- rbind(drivermutations.knockdown.prob, drivermutations.knockdown)
pca = prcomp(PCA2, center = F, scale. = F)
summary(pca)
plot(pca, type = "l")
library(ggbiplot)
ggbiplot(pca)
ggbiplot(pca, choices = 2:3)
ggbiplot(pca, labels=rownames(pca$x))