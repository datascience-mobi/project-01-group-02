# welche Gene mutieren viel (oben im Boxplot) und welche wenig (unten im Boxplot) wie sind sie generell verteilt 

# zu sehen wie sind Expressionsvalues verteilt in einer Zelllinie

boxplot(allDepMap_expression_SkinCancer$`ACH-000014`, main = "ACH-000014", xlab = "Expressionvalue", horizontal = T)
boxplot(allDepMap_expression_SkinCancer , main = "overview", xlab= "cell  lines", ylab = "values")

     


         



# Besprechung mit David: 
# alle Hugo name in einen Vektor 
# mit Funktion count, unique und sapply neue Matrix erstellen mit den Counts von jedem einzelnen Gen
# Beispiel f?r apply family ist in der Gruppe 
# clean up ist eig nicht n?tig aber k?nnen es auch so lassen
# nachdem mit den Counts die gene ausgesucht wurden die wichtig sind neue Expression matrix
# extra matrix mit den Gene die viel mutieren und den genen die nie mutiert sind(die die in der unges?urberten Mutationsmatrix nicht auftauchen)
# mit melt funktion daten bearbeiten die mit ggplot dargestellt werden sollen 
# f?r SDL k-means clustering
# die packeges f?r regressionsanalyse etc erfragen 
     
     #erstellen der Matrizen mit den Driver Mutations (BRAF, NFI, WT, RAS) aus der Expression Daten
     #Davor erstmal alphabetisch Ordnen der rows 
     allDepMap_expression_SkinCancer = allDepMap_expression_SkinCancer[order(rownames(allDepMap_expression_SkinCancer)),]
     #für BRAF
     which(rownames(allDepMap_expression_SkinCancer) == "BRAF")
     #um herauszufinden in welchen Zelllinien BRAF exprimiert wird, [1] 10195 das kommt dabei raus 
     BRAFexpression <- allDepMap_expression_SkinCancer[which(rownames(allDepMap_expression_SkinCancer) == "BRAF"),]
     #alle Zeilen mit BRAF habe ich dann BRAFexpression genannt
     boxplot(BRAFexpression)
     #Boxplot zu BRAF
     #Als nächstes RAS, habe erstmal in der Liste geschaut welche Kürzel zu RAS passen würden, denn wenn man nur nach RAS kommt nichts raus 
     which(rownames(allDepMap_expression_SkinCancer) == "RARS")
     which(rownames(allDepMap_expression_SkinCancer) == "RASSF9")
     #um die Zeilen, in denen die RAS Gene sind zu identifizieren, war dan von 27136 bis 27181
     RASexpression <- allDepMap_expression_SkinCancer[27136:27181,]
     
     #habe dann alles als RASexpression zusammengefasst und Boxplot zeigen lasssen
     boxplot(RASexpression)
     #für NF1 der gleiche Ablauf
     which(rownames(allDepMap_expression_SkinCancer) == "NF1")
     which(rownames(allDepMap_expression_SkinCancer) == "NF1P8")
     NF1expression <- allDepMap_expression_SkinCancer[22721:22729,]
     #NF1expression mit allen NF1 vorhanden
     boxplot(NF1expression)
     #Last but not least: WT 
     which(rownames(allDepMap_expression_SkinCancer) == "WT1")
     which(rownames(allDepMap_expression_SkinCancer) == "WTIP")
     WTexpression <- allDepMap_expression_SkinCancer[47802:47806,]
     boxplot(WTexpression)
     #Nun haben wir alle unsere Infos einzelt extrahiert, zum zusammenfügen hab ich dann diesen Code benutzt:
     rbind(BRAFexpression, RASexpression, NF1expression, WTexpression)
     #Das Neue habe ich dann als Drivermutation.expr umbenannt
     Drivermutation.expr <- rbind(BRAFexpression, RASexpression, NF1expression, WTexpression)
     # dann hab ich noch einen wunderschönen Boxplot aus den Daten erstellt 
     boxplot(Drivermutation.expr)
     #das wars mic-drop :)
     
     #RASexpression varianz
     RASexpression1 = apply(RASexpression, 1, function(x){var(x)})
     RASexpression1
     hist(RASexpression1)
     #NF1expression varianz
     NF1expression1 = apply(NF1expression, 1, function(x){var(x)})
     hist(NF1expression1)
     #WTexpression varianz
     WTexpression1 = apply(WTexpression, 1,function(x){var(x)})
     hist(WTexpression1)
     
     
     #creating a Barplot with the data of the Drivermutationmatrix
     #so we can see the expression of the different genes over all celllines
     
     # first barplot has the cellines in the x-axis and the expression of all mutations in the y-axis 
     #so we see the expression of all genes in the celllines
     barplot(as.matrix(Drivermutation.expr[-1,-1]))
    
     #We want to have the genes in the x-axis so i created a new matrix with the rows and columns of the drivermutation.expr swapped
     df2 <- data.frame(t(Drivermutation.expr[-2, -2]))
     
     #creating a barplot with this new matrix df2
     #the genes are the x-axis and the expression is the y-axis over all celllines
     #the lines in the bars are the celllines -> so we can see in which celline how much of this gene is expressed
     barplot(as.matrix(df2[-1,-1]))
     
     #boxplot over the same data 
     boxplot(as.matrix(df2[-1, -1]))
     
     