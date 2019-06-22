
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
# zwischen den verschiedenen Genen, also werden immer die Zeilen verglichen in einem paired ttest
# Hypothese: es gibt keine Abweichen des Überlebens wenn in den gleichen Zellen einmal ein Driver und 
# einmal ein Target gene ausgenockt wird.

#das Problem ist das die Werte in unsere Dataframe nicht numerisch sind 
is.numeric(allDepMap_kd.ceres_SkinCancer["RASA3",])
is.numeric(as.numeric(allDepMap_kd.ceres_SkinCancer["RASA3",]))

#ttest mit einem einzelnen Drivergene über alle "potentiellen target genes"

t.test(as.numeric(allDepMap_kd.ceres_SkinCancer["RASA3",]), as.numeric(allDepMap_kd.ceres_SkinCancer[3,]), paired = TRUE,alternative = c("two.sided"))$statistic


# übertragen auf alle Drivergenes aus dem ttestgene vektor 
# erstellen eines Dataframes das alle Drivermutation Überlebenswerte mit denen der targetgenes vergleicht 
a <- as.list(c(ttestgenes))


for (i in c(1:length(ttestgenes))){
  k <- ttestgenes[i]
  a[[i]] <- sapply(1:length(rownames(allDepMap_kd.ceres_SkinCancer)), function(y){t.test(as.numeric(allDepMap_kd.ceres_SkinCancer[ttestgenes[i],]), as.numeric(allDepMap_kd.ceres_SkinCancer[y,]), paired = TRUE,alternative = c("two.sided"))$statistic
  })}

# unser Dataframe ist nun eine List aus ganz vielen dataframes
# jedes dieser Dataframes wiederum gibt die Kombination eines Drivergenes mit allen targetgenen wieder 
# im folgenden werden die Elemente des dataframes benannt

for (i in c(1:length(ttestgenes))){
  names(a)[[i]] <- ttestgenes[i]
}
ttestdata <- a

rm(a)

for (i in c(1:length(ttestdata))) {
  names(ttestdata[[i]]) <- rownames(allDepMap_kd.ceres_SkinCancer)
}

# der nächste Schritt wäre nun die t Werte so zu ordnen das wir die wo am meisten von 0 
# abweichen rauszeihen können
# damit wir später sagen können welche Werte signifikant von unserer H0 Hypthese abweichen 





save(file= "C:/Users/LeoTh/Documents/GitHub/project-01-group-02/ttestdataframe.RDS", list="ttestdata", "allDepMap_kd.ceres_SkinCancer", "ttestgenes")
load("C:/Users/LeoTh/Documents/GitHub/project-01-group-02/ttestdataframe.RDS")

# an dieser Stelle bin ich dann ein wenig verloren...
# wir haben ja hier die Mittelwerte über alle Zelllinie bei einem Knockout eines 
# spezifischen Genes betrachtet und diese Mittelwerte verglichen 
# das heist für stark von 0 abweichende Werte könenn wir sagen , dass wir die H0 Hypothese 
# verwerfen können
# aber wie sagt mir das etwas über die Interaktion ? wenn die Überlebenschancen signifikant 
# voneinander abweichen? oder müssen wir die Genen betrachten für welche die H0 Hypothese 
# zutrifft ? weil das Zellwachstum bei den beiden beiden Genen in gleicher weise voneinander abweicht
# oder müssen wir mit der Korrelation arbeiten ? 
# ist der Ansatz so wie wir Ihn gewählt haben so überhaupt richtig ? 





