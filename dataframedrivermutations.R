# the next task is to extract the occupancy of each gene mutation 

mutationoccupancy <- as.data.frame(table(Mutation_1dataframe$Hugo_Symbol))


mutationoccupancy <- mutationoccupancy[order(mutationoccupancy$Freq, decreasing = TRUE),]

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

b <- countingMutations[3500:3526,]
b$Var1 <- factor(b$Mutations, levels = b$Mutations)

ggplot(data = b)+
  (geom_bar(mapping = aes(x = Mutations, y = Occupancy), stat = "identity"))


Targetmutations <- countingMutations



