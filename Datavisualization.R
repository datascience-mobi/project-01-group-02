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
     










which(allDepMap_mutation_SkinCancer$`ACH-000014`[,2]=="AL627309.1")

elif = names(allDepMap_mutation_SkinCancer)





#gibt mir aus in welchen Zelllines es an welchen Positionen genau diese Mutation gibt
i = 0
while(i<35){i = i +1; print(which(allDepMap_mutation_SkinCancer[[elif[i]]][ ,2]== "SHROOM2"));}

# gibt mir an in welcher Anzahl diese mutation auftritt in allen Zelllinien
while(i<length(elif)+1){i = i +1; print(length(which(allDepMap_mutation_SkinCancer[[elif[i]]][ ,2]== "SHROOM2")));}

# gibt mir an in welcher Anzahl bestimmte Mutationen in einer Zelllinie auftreten
a = rownames(allDepMap_expression_SkinCancer)
while(i<length(a)+1){i = i +1; print(length(which(allDepMap_mutation_SkinCancer$`ACH-000274`$Hugo_Symbol== a[i])));}

#ich brauche jtzt also als erstes einen Vektor a der wirklich alle möglichen Mutationen enthält 
#als a werde ich die aufgelisteten Mutationen aus der Expressionsmatrix nehmen


#wenn ich das dan habe werde ich die 2. Schleif als äußere Schleife benutzen und die 1. als innere
#dadurch kann ich mir Reihe pro Reihe eine neue Dataframe erstellen welche mir wiedergibt wie 
#viele Mutationen eine bestimmten Gens in jeder Celllinie vorliegt 
# die Zelllinien werden die Spalten sein 
# und die Reihen die Mutationen 
#aber 

# die namen der Celllinien hat eig die gleiche auswirkung wie Vektor c(1:34)



# Namen aller Gene welche auch in unsere Expressionsmatrix eine Rolle spielen 



#while(b< length(a)+1){b = b+1; i = 0;  while(i<length(elif)+1){i = i +1; print(length(which(allDepMap_mutation_SkinCancer[[elif[i]]][ ,2]== a[b])));}}

# mir wird hier allerdings deutlich zu wenig ausgegeben
# war ein fehler von mir da ich das i nicht jedes mal wieder auf 0 gesetzt habe sondern nur einmal auserhalb 
# while schleife und dann nichtmehr daher ist die innere Schleife nur genau einmal abgelaufen weil i schon zu anfang nicht mehr <34 war

a = rownames(allDepMap_expression_SkinCancer)
elif = names(allDepMap_mutation_SkinCancer)
i = 1
b = 1
Anzahlmutation <- c()

while(b< length(a)){
    i = 1
    Spalte <- c()
    while(i<length(elif)){
        Spalte[i] <- length(which(allDepMap_mutation_SkinCancer[[elif[i]]][ ,2]== a[b]))
        i = i +1
        }
    Anzahlmutation[b] <- Spalte
    b = b+1
    }

# villt stimmen die Reihennamen der Expressionsmatrix einfach nicht mit denen der mutationsmatrix überein 
#wie soll ich sie sonst raussuchen, wir wollen ja nämlich ein Histogramm erstellen die uns die Anzahl de Mutationen
#von ienem bestimmten gen über alles Zelllinien anzeigen 
# villt mit der Annotations spalte ??

# im obigen Text wurde von uns versucht durch die einzelnen dataframes jeder einzelnen Zelllinie zu gehen und von jedem einzelnen Gen 
# rauszusuchen wie oft es mutiert ist 
# also schleife in einer schleife innere schleife geht mit einem bestimmten gen durch alle Zelllinien und gibt die Anzahl der Mutataionen wieder 
# dies soll in einem vektor gespeichert werden  
# in der äußeren schleife soll dann dieser Vektor zu einem dataframe ergänzt werden und durch alle Gennamen gegangen werden das bei jedem 
#durchlauf der inneren schleife nach einem anderen Gen gesucht wird 
# Problem ist hierbei das wir ja keine Liste aller mutierter Gene haben wie erstellen wir eine solche List ? 
# denken wir villt viel zu kompliziert? gibt es einen einfacherern WEg das Histogramm darzustellen?

# im Grunde sind unsere Hauptfragen:
# wie nutze ich ggplot 
# wie können wir die mutationmatrix so auslesen das wir nur die anzahl der Mutationen über alle Zelllinien bekommen 
# wie ziehen wir eine verbindung zwischen den Bezeichungen der Mutationen und den Namen der Gene in der Expressionsmatrix für die Arbeit nach der Visualisierung ? 


