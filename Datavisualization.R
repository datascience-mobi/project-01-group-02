# welche Gene mutieren viel (oben im Boxplot) und welche wenig (unten im Boxplot) wie sind sie generell verteilt 

# zu sehen wie sind Expressionsvalues verteilt in einer Zelllinie

boxplot(allDepMap_expression_SkinCancer$`ACH-000014`, main = "ACH-000014", xlab = "Expressionvalue", horizontal = T)
boxplot(allDepMap_expression_SkinCancer , main = "overview", xlab= "cell  lines", ylab = "values")

ggplot(data = allDepMap_expression_SkinCancer[1:4, ])  
     geom_smooth(mapping = aes(x = colnames, y = rownames(allDepMap_expression_SkinCancer[1:4, ]),  group = rownames(allDepMap_expression_SkinCancer[1:4, ]) )
                   )
#wie mache ich das hier es gibt nich so schöne kategorien wie bei dem Beispiel also wie benenne ich mien x und y ???
# waren ein wenig ratlos wie genau man das ggplot paket benutzt
# klar mit den beispielen aus dem R for datascience funktioniert das alles wunderbar aber wir haben hier ja nicht so 
# kategorien die ich definieren kann , also was setze ich für die aesthetics ein ? 
     

# der nächste Schritt wäre jtzt das erstellen eines Histogramms, welches zeigt welche Gene am öftesten mutieren 
# wie machen wir das ? 
# haben bisher im Clean up eine dataframe erstellt der uns alle Mutationen über alle Zelllinien wiedergibt 
     # also Gene die gleich heisen stehen direkt untereinander 
# das sind aber echt nicht viele die genau gleich sind ? 
# sollen wir also wirklich nur alle die genau gleich heißen als eine Mutation im Gen zusammenfassen oder heißen 
# villt Mutaitonen anderster obwohl sie im selben Gen liegen? wie sehen wir das dann weil wir haben zwar 
# die Startposition und endposition von den Mutationen aber ja nicht von den Genen
# oder gibt es einfach so viel Passenger mutations das wir uns dann einfach nur auf die nicht uniquen Genemutationen konzentrieren sollen ??
     
         
length(unique(Mutation_1dataframe$Hugo_Symbol))


# im Grunde sind unsere Hauptfragen:
# wie nutze ich ggplot 
# wie können wir die mutationmatrix so auslesen das wir nur die anzahl der Mutationen über alle Zelllinien bekommen 
# wie ziehen wir eine verbindung zwischen den Bezeichungen der Mutationen und den Namen der Gene in der Expressionsmatrix für die Arbeit nach der Visualisierung ? 
# ist die Säuberung so okee es wurden halt für alle Zelllinien die Annotations beispielsweise gelöscht 
# zum erstellen der großen Matrix da einige Zelllinien diese nicht hatten

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
