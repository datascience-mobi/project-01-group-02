

# k means 
# mit welchen daten ist es sinnvoll einen kmeans zu machen? 
# wir vegleichen ja immer die Zelllinien untereinander also ein Punkt ist eine Zelllinen
# also villt die Zelllinien so umbenen dass man sieht welche Driver Mutations dort Mutiert sind
# wie machen wir das weil die Drivermutations tauchen ja nicht in der Mutationsmatrix auf also 
# villt die Zelllinien nach den Drivermutations benennnen die die stärkste Abweichung von 
# den anderen Zelllinien haben in der Drivermutation
# oder wir sagen wir benennen sie nach den Mutationen der Driver Mutations obwohl diese so selten sind 
# es gibt so viele verschiedene RAS deswegen ist es villt egal das sie so selten auftauchen
# Expressionsmatrix von taget genes mit den driver mutations
# von der variation der Gene ? also die Expressionsmatrix normalisieren und dann davon K-means 
# von copynumber daten 
# die Hoffnung bei diesen beiden ist das wir ein k-means rausbekommen bei denen die 
# mit den gleichen DriverMuatations in die gleichen Cluster kommen


# Zelllinien welche eine mutation in den Drivermutations haben 
namesdrivers <- c(rownames(BRAFexpression),rownames(NF1expression), rownames(RASexpression), rownames(WTexpression))
zuordnenZelllinienDM <- Mutation_1dataframe[which(Mutation_1dataframe$Hugo_Symbol %in% namesdrivers),]
