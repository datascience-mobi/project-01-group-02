# welche Gene mutieren viel (oben im Boxplot) und welche wenig (unten im Boxplot) wie sind sie generell verteilt 

# zu sehen wie sind Expressionsvalues verteilt in einer Zelllinie

boxplot(allDepMap_expression_SkinCancer$`ACH-000014`, main = "ACH-000014", xlab = "Expressionvalue", horizontal = T)
boxplot(allDepMap_expression_SkinCancer , main = "overview", xlab= "cell  lines", ylab = "values")

ggplot(data = allDepMap_expression_SkinCancer[1:4, ])  
     geom_smooth(mapping = aes(x = colnames, y = rownames(allDepMap_expression_SkinCancer[1:4, ]),  group = rownames(allDepMap_expression_SkinCancer[1:4, ]) )
                   )

     


         



# Besprechung mit David: 
# alle Hugo name in einen Vektor 
# mit Funktion count, unique und sapply neue Matrix erstellen mit den Counts von jedem einzelnen Gen
# Beispiel für apply family ist in der Gruppe 
# clean up ist eig nicht nötig aber können es auch so lassen
# nachdem mit den Counts die gene ausgesucht wurden die wichtig sind neue Expression matrix
# extra matrix mit den Gene die viel mutieren und den genen die nie mutiert sind(die die in der ungesäurberten Mutationsmatrix nicht auftauchen)
# mit melt funktion daten bearbeiten die mit ggplot dargestellt werden sollen 
# für SDL k-means clustering
# die packeges für regressionsanalyse etc erfragen 
