


#show rownames of matrix to compare if they all show the genes in rownames
rownames(Drivermutation.expr)
rownames(allDepMap_kd.ceres_SkinCancer)
rownames(allDepMap_kd.prob_SkinCancer)
rownames(allDepMap_copynumber_SkinCancer)
rownames(allDepMap_expression_SkinCancer)


# we generated a new matrix (drivermutations.knockdown) which contains the knockdown data of the mutated genes from Drivermutation.expr
drivermutations.knockdown = allDepMap_kd.ceres_SkinCancer[which(rownames(allDepMap_kd.ceres_SkinCancer) %in% rownames(Drivermutation.expr)),]

dim(drivermutations.knockdown)
# thats what we get out: [1] 42 34

# we generated a new matrix (drivermutations.knockdown.prob) which contains the knockdown.prob data of the mutated genes from Drivermutation.expr#
drivermutations.knockdown.prob = allDepMap_kd.prob_SkinCancer[which(rownames(allDepMap_kd.prob_SkinCancer) %in% rownames(Drivermutation.expr)),]

dim(drivermutations.knockdown.prob)
# thats what we get out [1] 42 34

# we generated a new matrix (drivermutations.copynumber) which contains the copynumber data of the mutated genes from Drivermutation.expr
drivermutations.copynumber = allDepMap_copynumber_SkinCancer[which(rownames(allDepMap_copynumber_SkinCancer) %in% rownames(Drivermutation.expr)),]

dim(drivermutations.copynumber)
# thats what we get out: [1] 48 34

# we generated a new matrix (drivermutations.expression) which contains the expression data of the mutated genes from Drivermutation.expr
drivermutations.expression = allDepMap_expression_SkinCancer[which(rownames(allDepMap_expression_SkinCancer) %in% rownames(Drivermutation.expr)),]
 
dim(drivermutations.expression)
# thats what we get out: [1] 61 34

#generating a heatmap of the knockdown data and the knockdown.prob data
heatmap(as.matrix(drivermutations.knockdown))
heatmap(as.matrix(drivermutations.knockdown.prob))

