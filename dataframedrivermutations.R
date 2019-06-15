
load("C:/Users/LeoTh/Documents/GitHub/project-01-group-02/CellCancerLinesafterCleanup.RDS")


# the next task is to extract the occupancy of each gene mutation 

#mutationoccupancy <- as.data.frame(table(Mutation_1dataframe$Hugo_Symbol))


#mutationoccupancy <- mutationoccupancy[order(mutationoccupancy$Freq, decreasing = TRUE),]


library(tidyverse)
library(reshape2)


# building a dataframe with the counts of mutated cells per mutation 



countingMutations <- table(Mutation_1dataframe$Hugo_Symbol)
countingMutations <- melt(countingMutations)
countingMutations = countingMutations[order(countingMutations$value),]

dim(countingMutations)

rownames(countingMutations)
rownames(countingMutations) <- c(1:3526)

colnames(countingMutations) <- c("Mutations", "Occupancy")


countingMutations$Mutations <- factor(countingMutations$Mutations, levels = countingMutations$Mutations)

ggplot(data = countingMutations)+
  (geom_bar(mapping = aes(x = Mutations, y = Occupancy), stat = "identity"))

targetmutations <- countingMutations[which(countingMutations$Occupancy >6),]
targetmutations$Var1 <- factor(targetmutations$Mutations, levels = targetmutations$Mutations)

ggplot(data = targetmutations)+
  (geom_bar(mapping = aes(x = Mutations, y = Occupancy), stat = "identity"))



namestargetmutations <- targetmutations$Mutations 

targetexpression <- allDepMap_expression_SkinCancer[which(rownames(allDepMap_expression_SkinCancer)%in% namestargetmutations),]
targetcopynumber <- allDepMap_copynumber_SkinCancer[which(rownames(allDepMap_copynumber_SkinCancer)%in% namestargetmutations),]
targetkd.ceres <- allDepMap_kd.ceres_SkinCancer[which(rownames(allDepMap_kd.ceres_SkinCancer)%in% namestargetmutations),]
targetkd.prob <- allDepMap_kd.prob_SkinCancer[which(rownames(allDepMap_kd.prob_SkinCancer)%in% namestargetmutations),]


# notargets <- countingMutations$Mutations[which(countingMutations$Occupancy < 7)]

# targetexpression2 <- allDepMap_expression_SkinCancer[-which(rownames(allDepMap_expression_SkinCancer) %in% notargets),]


save(file= "C:/Users/LeoTh/Documents/GitHub/project-01-group-02/targetdata.RDS", list="countingMutations", "targetexpression", "targetcopynumber", "namestargetmutations", "targetkd.ceres", "targetkd.prob","Mutation_1dataframe")
load("C:/Users/LeoTh/Documents/GitHub/project-01-group-02/targetdata.RDS")

