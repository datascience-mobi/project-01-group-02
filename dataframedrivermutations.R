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

b<-countingMutations[3520:3526,]

ggplot(data = b) +
  geom_bar(mapping = aes(x = Var1, y = value), stat = "identity")


