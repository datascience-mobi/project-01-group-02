
# mit dem T-test kann geschaut werden ob sich zwei Stichproben seginifikant 
# voneinander unterscheinden
# Daten: wir brauchen die Namen der Driver Genes welche wir uns anschauen

load("C:/Users/LeoTh/Documents/GitHub/project-01-group-02/driverdata.RDS")
ttestgenes <- rownames(driverkd.ceres)
length(ttestgenes)
rm(countingMutations,drivercopynumber, driverexpression, driverkd.ceres, driverkd.prob, namesdrivermutations, Mutation_1dataframe)
# wir brauen das komplette ceres dataframe: 
load("C:/Users/LeoTh/Documents/GitHub/project-01-group-02/CellCancerLines.RDS")
rm(allDepMap_annotation_SkinCancer, allDepMap_copynumber_SkinCancer, allDepMap_expression_SkinCancer, allDepMap_kd.prob_SkinCancer, allDepMap_mutation_SkinCancer)

# Frage: weichen die Mittelwerte der Überlebenswerte von den Genen seknifikant voneinander ab
# zwischen den verschiedenen Genen, also werden immer die Zeilen verglichen, je gr
# Hypothese: es gibt keine Abweichen des Überlebens wenn in den gleichen Zellen einmal ein Driver und 
# einmal ein Target gene ausgenockt wird.

is.numeric(allDepMap_kd.ceres_SkinCancer["RASA3",])
is.numeric(as.numeric(allDepMap_kd.ceres_SkinCancer["RASA3",]))

t.test(as.numeric(allDepMap_kd.ceres_SkinCancer["RASA3",]), as.numeric(allDepMap_kd.ceres_SkinCancer[3,]), paired = TRUE,alternative = c("two.sided"))$statistic



a <- as.list(c(ttestgenes))


for (i in c(1:length(ttestgenes))){
  k <- ttestgenes[i]
  a[[i]] <- sapply(1:length(rownames(allDepMap_kd.ceres_SkinCancer)), function(y){t.test(as.numeric(allDepMap_kd.ceres_SkinCancer[ttestgenes[i],]), as.numeric(allDepMap_kd.ceres_SkinCancer[y,]), paired = TRUE,alternative = c("two.sided"))$statistic
  })}

for (i in c(1:length(ttestgenes))){
  names(a)[[i]] <- ttestgenes[i]
}
ttestdata <- a

rm(a)

save(file= "C:/Users/LeoTh/Documents/GitHub/project-01-group-02/ttestdataframe.RDS", list="ttestdata", "allDepMap_kd.ceres_SkinCancer", "ttestgenes")
load("C:/Users/LeoTh/Documents/GitHub/project-01-group-02/ttestdataframe.RDS")



