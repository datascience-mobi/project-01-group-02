# Unsere Originalen Ausgangsdaten laden: 

# allDepMapData = readRDS("C:/Users/LeoTh/Documents/GitHub/project-01-group-02- gesamte Daten/DepMap19Q1_allData.RDS")

# wollen wir hier aber nicht machen sondern wir befassen uns mit den Daten welche nur unsere Zelllinien betreffen und säubern diese

# das ist unser Datensatz mit unseren Zelllinien, er ist im Github Ordner gespeicher 
# Achung : ihr habt einen anderen Pfad als ich 

load("C:/Users/LeoTh/Documents/GitHub/project-01-group-02/CellCancerLines.RDS")

## Data Cleaning 

## Mutation Dataframe 

# alle Reihen entfernen die nicht Deleterious sind 

a = 1
while (a < length(names(allDepMap_mutation_SkinCancer))) { 
  allDepMap_mutation_SkinCancer[[a]] = allDepMap_mutation_SkinCancer[[a]][allDepMap_mutation_SkinCancer[[a]][,"isDeleterious"]==TRUE, ]
  a = a+1
}




# als erstes haben wir uns die NAs pro Reihe also die Missing values pro Mutation angeschauen 
# aber wie bei der annotation matrix sollten wir nicht schauen welche Mutationen viele Nas haben 
# sondern welche observations in verschiedenen Mutationen missing value lifern 
# also nicht alle NAs pro Reihe(was dann die Mutationen wären) sondern pro Observables (hier also die Spalten )
#rmv.rows = apply (allDepMap_mutation_SkinCancer$`ACH-000274`, 1, function (x){sum(is.na(x))})
#rmv.rows = apply (allDepMap_mutation_SkinCancer$`ACH-000274`, 2, function (x){sum(is.na(x))})

# rmv.rows speichert nun alle spalten für die Nas auftretten 
#allDepMap_mutation_SkinCancer$`ACH-000274` = allDepMap_mutation_SkinCancer$`ACH-000274`[, -which(rmv.rows>0)]
# diese Spalten werden aus der Datenmatrix entfernt 
# das ist aber nur eine zellline jtzt müssen wir es noch für alle 34 Zelllinien machen 




i = 1
while (i < length(names(allDepMap_mutation_SkinCancer))) {

  rmv.rows = apply (allDepMap_mutation_SkinCancer[[i]], 2, function (x){sum(is.na(x))})
  allDepMap_mutation_SkinCancer[[i]] = allDepMap_mutation_SkinCancer[[i]][, -which(rmv.rows>0)]
  i = i+1
}
# R script hatte ein Problem damit das wir immer mit i = 0 angefangen haben
# daher haben wir jtzt mit  i= 0 gemacht und das +1 in der schleife erst die Aktion ausgeführt ist positioniert

# müssen noch Spalten entfernt werden??

# die Mutationen in alphabetischer Reihenfolge
a = 1
while(a< 35){
  allDepMap_mutation_SkinCancer[[a]] = allDepMap_mutation_SkinCancer[[a]][order(allDepMap_mutation_SkinCancer[[a]][,2]), ]
   a = a+1
}

# um die einzelnen Listen in eine Dataframe zu machen müssen wir erst checken ob sie die gleichen Spaten haben 
a = 1
 while (a<35) {print(length(colnames(allDepMap_mutation_SkinCancer[[a]]))); a = a+1;}
# um zu sehen was die niedrigste Spaltenanzahl ist , die niedrigste Spalten anzahl hatt Zellinie ACH-001550
# wie beziehen also die Spalten auf alle anderne Zelllinien und löschen den Rest raus 
# diese Zelllinie hat den Index 33 in der Zelllinienliste
a = 1
while (a<35) {
  allDepMap_mutation_SkinCancer[[a]] = allDepMap_mutation_SkinCancer[[a]][,colnames(allDepMap_mutation_SkinCancer[[33]])]
  a = a+1
}


#alles in eine Dataframe packen und dann in alphabetischer Reihenfolge angschauen

mymatrix <- rbind(allDepMap_mutation_SkinCancer[[1]],allDepMap_mutation_SkinCancer[[2]])

a = 3
 while (a<35) {mymatrix <- rbind(mymatrix,allDepMap_mutation_SkinCancer[[a]])
    a = a+1   
 }

Mutation_1dataframe <- mymatrix
rm(mymatrix)


# kontrolle ob die Dimmension der neuen Dataframe mit den Anzahl der Mutationnen 
# in allen Zelllinien übereinstimmt 
a = 1
Reihen = 0 
while (a<35) { 
    Reihen = Reihen + length(rownames(allDepMap_mutation_SkinCancer[[a]]))
    a = a+1
}

# jtzt nach alphabet ordnen damit man sieht welche Mutation die gleich sind untereinander stehen

Mutation_1dataframe = Mutation_1dataframe[order(Mutation_1dataframe$Chromosome),]

Mutation_1dataframe = Mutation_1dataframe[order(Mutation_1dataframe$Start_position),]

Mutation_1dataframe = Mutation_1dataframe[order(Mutation_1dataframe$Hugo_Symbol),]


## Annotation Data Matrix

#haben uns den Dataframe angeschaut
View(allDepMap_annotation_SkinCancer)

# Reihennamen nach den Zelllinien benennnen
 rownames(allDepMap_annotation_SkinCancer)= allDepMap_annotation_SkinCancer$DepMap_ID
 
# damit wird die Spalte DepMap_ID überflüssig 
# wir haben außerdem beschlossen auch die Spalten Aliases ,Gender ,Source aus dem Dataframer herauszunehmen
# weil wir die Info erstmal nicht brauchen
 
allDepMap_annotation_SkinCancer = allDepMap_annotation_SkinCancer[ , -which(colnames(allDepMap_annotation_SkinCancer) %in% c("DepMap_ID", "Aliases", "Gender", "Source")) ]

# dann wollen wir die Reihen umsortieren damit die gleichen Typen von Cancer beieinande stehen 
which(allDepMap_annotation_SkinCancer$Subtype.Disease != "Melanoma")
allDepMap_annotation_SkinCancer = allDepMap_annotation_SkinCancer[c(1:13, 15:19, 21:26, 28:30, which(allDepMap_annotation_SkinCancer$Subtype.Disease != "Melanoma")), ]

# die daten nominal machen, d.h wir schmeißen einfach die zusätzlichen Kategorien raus die wir nichtmehr brauchen
summary(allDepMap_annotation_SkinCancer)

allDepMap_annotation_SkinCancer$Subtype.Disease = factor(allDepMap_annotation_SkinCancer$Subtype.Disease)
allDepMap_annotation_SkinCancer$Primary.Disease = factor(allDepMap_annotation_SkinCancer$Primary.Disease)
summary(allDepMap_annotation_SkinCancer)

## kd.ceres Data Matrix 

dim(allDepMap_kd.ceres_SkinCancer)
# wir speichern auf der variable rmv.rows die Summe aller nicht NA Werte pro Reihe
# Es werden uns für alle Reihen also für jedes Gen angegeben in wie vielen Zelllinen keine Werte vorhanden sind 
# da unser Vektor uns nur die Anzahl angibt müssen wir nach jenen suchen  für die der Wert der fehlenden Werte gleich 0 ist und diese dann aussortieren 
rmv.rows = apply(allDepMap_kd.ceres_SkinCancer, 1, function(x) {sum(is.na(x))})
length(allDepMap_kd.ceres_SkinCancer$`ACH-000014`) == length(rmv.rows)
length(which(rmv.rows ==0 )) == length(rownames(allDepMap_kd.ceres_SkinCancer))
# Vektor aus allen Reihen welche keine Missing Values aufweisen ist genauso lange wie der Vektor aller Reihennamen 
# also haben wir keine NAs in unserer Dataframe 

allDepMap_kd.ceres_SkinCancer = allDepMap_kd.ceres_SkinCancer[order(rownames(allDepMap_kd.ceres_SkinCancer)), ]

## kd.prob Data Matrix

# dasselbe wie bei der kd.ceres Dataframe 
# also keine Missing values 
rmv.rows = apply(allDepMap_kd.prob_SkinCancer, 1, function(x) {sum(is.na(x))})
  length(allDepMap_kd.prob_SkinCancer$`ACH-000014`) == length(rmv.rows)
  length(which(rmv.rows ==0 )) == length(rownames(allDepMap_kd.prob_SkinCancer))


  allDepMap_kd.prob_SkinCancer = allDepMap_kd.prob_SkinCancer[order(rownames(allDepMap_kd.prob_SkinCancer)), ]
  



    
## Copynumber Dataframe 

# wir schmeisen die Gene mit NAs raus da sie dann nicht über die verschiedene Zelllinien vergleichbar sind 
# Gene = Zeilen 
    
elif = apply(allDepMap_copynumber_SkinCancer, 1, function(x) {sum(is.na(x))})
allDepMap_copynumber_SkinCancer = allDepMap_copynumber_SkinCancer[-which(elif>0),]

# ist so fertig denke ich 
# generell mal gedanken darüber machen ob wir niedrige Varianzen herausfiltern sollten

#alphabetische Sortierung für bessere Übersicht 
allDepMap_copynumber_SkinCancer = allDepMap_copynumber_SkinCancer[order(rownames(allDepMap_copynumber_SkinCancer)), ]


## Expressions Dataframe

elif = apply(allDepMap_expression_SkinCancer, 1, function(x) {sum(is.na(x))})
sum(isTRUE(which(elif>0)))

# als Wert wird hier 0 ausgegeben also gibt es keine NAs in diesem Dataframe
#alphabetische Sortierung für bessere Übersicht 

allDepMap_expression_SkinCancer = allDepMap_expression_SkinCancer[order(rownames(allDepMap_expression_SkinCancer)), ]



## Speichern der gesäuberten Daten

save(file= "C:/Users/LeoTh/Documents/GitHub/project-01-group-02/CellCancerLinesafterCleanup.RDS", list="allDepMap_annotation_SkinCancer", "allDepMap_copynumber_SkinCancer", "allDepMap_expression_SkinCancer", "allDepMap_kd.ceres_SkinCancer", "allDepMap_kd.prob_SkinCancer", "allDepMap_mutation_SkinCancer")



#load("C:/Users/LeoTh/Documents/GitHub/project-01-group-02/CellCancerLinesafterCleanup.RDS")


# Fragen: 

# low variance Ideen ? also müssen wir uns noch unterhalten wie wir das mit der Varianz machen 
# villt nur mit dem Histogramm aus den oft mutierten Genen arbeiten und da die Gene raussuchen die Wichtig sind 
#und die Wichtigkeit nicht von der Varianzen abhängig machen




