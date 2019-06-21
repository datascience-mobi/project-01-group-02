


allDepMapData = readRDS("C:/Users/LeoTh/Documents/GitHub/project-01-group-02- gesamte Daten/DepMap19Q1_allData.RDS")


namescelllines = allDepMapData$annotation$DepMap_ID[which(allDepMapData$annotation$Primary.Disease == 'Skin Cancer')]

allDepMap_expression_SkinCancer =  allDepMapData$expression[ ,which(colnames(allDepMapData$expression) %in% namescelllines)]
allDepMap_mutation_SkinCancer =  allDepMapData$mutation[match(x = namescelllines, table = names(allDepMapData$mutation))]
allDepMap_copynumber_SkinCancer =  allDepMapData$copynumber[ ,which(colnames(allDepMapData$copynumber) %in% namescelllines)]
allDepMap_kd.ceres_SkinCancer =  allDepMapData$kd.ceres[ ,which(colnames(allDepMapData$kd.ceres) %in% namescelllines)]
allDepMap_kd.prob_SkinCancer =  allDepMapData$kd.prob[ ,which(colnames(allDepMapData$kd.prob) %in% namescelllines)]
allDepMap_annotation_SkinCancer =  allDepMapData$annotation[which(allDepMapData$annotation$DepMap_ID %in% namescelllines), ]



save(file= "C:/Users/LeoTh/Documents/GitHub/project-01-group-02/CellCancerLines.RDS", list="allDepMap_annotation_SkinCancer", "allDepMap_copynumber_SkinCancer", "allDepMap_expression_SkinCancer", "allDepMap_kd.ceres_SkinCancer", "allDepMap_kd.prob_SkinCancer", "allDepMap_mutation_SkinCancer")
load("C:/Users/LeoTh/Documents/GitHub/project-01-group-02/CellCancerLines.RDS")
rm(allDepMapData)


a = 1
while (a < length(names(allDepMap_mutation_SkinCancer))) { 
  allDepMap_mutation_SkinCancer[[a]] = allDepMap_mutation_SkinCancer[[a]][allDepMap_mutation_SkinCancer[[a]][,"isDeleterious"]==TRUE, ]
  a = a+1
}

i = 1
while (i < length(names(allDepMap_mutation_SkinCancer))) {
  
  rmv.rows = apply (allDepMap_mutation_SkinCancer[[i]], 2, function (x){sum(is.na(x))})
  allDepMap_mutation_SkinCancer[[i]] = allDepMap_mutation_SkinCancer[[i]][, -which(rmv.rows>0)]
  i = i+1
}

a = 1
while(a< 35){
  allDepMap_mutation_SkinCancer[[a]] = allDepMap_mutation_SkinCancer[[a]][order(allDepMap_mutation_SkinCancer[[a]][,2]), ]
  a = a+1
}

a = 1
while (a<35) {print(length(colnames(allDepMap_mutation_SkinCancer[[a]]))); a = a+1;}

a = 1
while (a<35) {
  allDepMap_mutation_SkinCancer[[a]] = allDepMap_mutation_SkinCancer[[a]][,colnames(allDepMap_mutation_SkinCancer[[33]])]
  a = a+1
}

mymatrix <- rbind(allDepMap_mutation_SkinCancer[[1]],allDepMap_mutation_SkinCancer[[2]])

a = 3
while (a<35) {mymatrix <- rbind(mymatrix,allDepMap_mutation_SkinCancer[[a]])
a = a+1   
}

Mutation_1dataframe <- mymatrix
rm(mymatrix)

length(unique(Mutation_1dataframe$Hugo_Symbol))

Mutation_1dataframe = Mutation_1dataframe[order(Mutation_1dataframe$Chromosome),]

Mutation_1dataframe = Mutation_1dataframe[order(Mutation_1dataframe$Start_position),]

Mutation_1dataframe = Mutation_1dataframe[order(Mutation_1dataframe$Hugo_Symbol),]

Mutation_1dataframe = Mutation_1dataframe[-which(Mutation_1dataframe$isDeleterious == FALSE),]

# wir wollen jtzt die gleichen Gene in eine Bezeichung zusammenfassen damit 
# sich damit besser arbeiten lässt

Mutation_1dataframe <- cbind(Mutation_1dataframe,Mutation_1dataframe$Hugo_Symbol)
Mutation_1dataframe <- Mutation_1dataframe[, c(1,2,17,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]
colnames(Mutation_1dataframe)[3] <- "Hugo_Symbols2"

#Mutation_1dataframe[grep('^ABCA.*', as.vector(Mutation_1dataframe$Hugo_Symbol), value = F),3] = "ABCA"

grep('^WT.*', as.vector(Mutation_1dataframe$Hugo_Symbol), value = F)
Mutation_1dataframe[grep("^WT.*", as.vector(Mutation_1dataframe$Hugo_Symbol)),"Hugo_Symbol"] <- "WT"

grep('^RAS.*', as.vector(Mutation_1dataframe$Hugo_Symbol), value = F)
Mutation_1dataframe[grep("^RAS.*", as.vector(Mutation_1dataframe$Hugo_Symbol)),"Hugo_Symbol"] <- "RAS"

grep('^RB.*', as.vector(Mutation_1dataframe$Hugo_Symbol), value = F)
Mutation_1dataframe[grep("^RB.*", as.vector(Mutation_1dataframe$Hugo_Symbol)),"Hugo_Symbol"] <- "RB"

grep('^NF.*', as.vector(Mutation_1dataframe$Hugo_Symbol), value = F)
Mutation_1dataframe[grep("^NF.*", as.vector(Mutation_1dataframe$Hugo_Symbol)),"Hugo_Symbol"] <- "NF"

grep('^BRAF.*', as.vector(Mutation_1dataframe$Hugo_Symbol), value = F)
Mutation_1dataframe[grep("^BRAF.*", as.vector(Mutation_1dataframe$Hugo_Symbol)),"Hugo_Symbol"] <- "BRAF"


#Mutation_1dataframe[grep("^[A-Za-z-0-9_]+(\@)?$", as.vector(Mutation_1dataframe$Hugo_Symbol)),"Hugo_Symbol"] <- "ABCA"

allDepMap_copynumber_SkinCancer = allDepMap_copynumber_SkinCancer[order(rownames(allDepMap_copynumber_SkinCancer)),]
allDepMap_expression_SkinCancer = allDepMap_expression_SkinCancer[order(rownames(allDepMap_expression_SkinCancer)),]

rownames(allDepMap_annotation_SkinCancer)= allDepMap_annotation_SkinCancer$DepMap_ID
allDepMap_annotation_SkinCancer = allDepMap_annotation_SkinCancer[ , -which(colnames(allDepMap_annotation_SkinCancer) %in% c("DepMap_ID", "Aliases", "Gender", "Source")) ]
which(allDepMap_annotation_SkinCancer$Subtype.Disease != "Melanoma")
allDepMap_annotation_SkinCancer = allDepMap_annotation_SkinCancer[c(1:13, 15:19, 21:26, 28:30, which(allDepMap_annotation_SkinCancer$Subtype.Disease != "Melanoma")), ]
allDepMap_annotation_SkinCancer$Subtype.Disease = factor(allDepMap_annotation_SkinCancer$Subtype.Disease)
allDepMap_annotation_SkinCancer$Primary.Disease = factor(allDepMap_annotation_SkinCancer$Primary.Disease)


rmv.rows = apply(allDepMap_kd.ceres_SkinCancer, 1, function(x) {sum(is.na(x))})
allDepMap_kd.ceres_SkinCancer = allDepMap_kd.ceres_SkinCancer[order(rownames(allDepMap_kd.ceres_SkinCancer)), ]


rmv.rows = apply(allDepMap_kd.prob_SkinCancer, 1, function(x) {sum(is.na(x))})
allDepMap_kd.prob_SkinCancer = allDepMap_kd.prob_SkinCancer[order(rownames(allDepMap_kd.prob_SkinCancer)), ]



save(file= "C:/Users/LeoTh/Documents/GitHub/project-01-group-02/CellCancerLinesafterCleanup.RDS", list="allDepMap_annotation_SkinCancer", "allDepMap_copynumber_SkinCancer", "allDepMap_expression_SkinCancer", "allDepMap_kd.ceres_SkinCancer", "allDepMap_kd.prob_SkinCancer", "allDepMap_mutation_SkinCancer", "Mutation_1dataframe")
load("C:/Users/LeoTh/Documents/GitHub/project-01-group-02/CellCancerLinesafterCleanup.RDS")

