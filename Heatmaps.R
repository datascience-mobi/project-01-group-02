#Heat map in der knochdown-Matrix
heatmap(allDepMap_kd.ceres_SkinCancer, scale = "rows")

#does not work -> matrix not numeric -> convert matric into numeric matrix
as.data.frame(lapply(allDepMap_kd.ceres_SkinCancer,as.numeric))

#converting the mtrix did not work
heatmap(allDepMap_kd.ceres_SkinCancer, scale = "none")

#using the matrix for the heatmap without converting the matrix
#into a numeric matrix
heatmap(as.matrix(allDepMap_kd.ceres_SkinCancer[, -1]))
#klammer zeigt an bis wohin die daten einbezogen werden -> [, -1] sagt, dass dei erste column nicht mit einbezogen wird aber alle rows

#create smaller matrix from the Knockout matrix
matrix.klein = allDepMap_kd.ceres_SkinCancer [1:20, 1:34]

#creating heatmap with this smaller matrix
heatmap(as.matrix(matrix.klein))

#creating smaller heatmap with probability-matrix
matrix.kleinprob = allDepMap_kd.prob_SkinCancer [1:20, 1:34]

#heatmap with smaller probability matrix
heatmap(as.matrix(matrix.kleinprob))
