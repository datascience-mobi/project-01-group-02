

#Importing the libs
library(ggplot2)
library(gridExtra)
library(reshape2)
library(data.table)
library(cluster)
library(rstudioapi)
library(pheatmap)
library(caret)
library(tidyverse)
library(dendextend)
library(factoextra)
library(devtools)
library(ggfortify)
library(rstudioapi)  
library(data.table) 
library(ggplot2) 
library(scales) 
library(stats)

#Setting the sys-path
root.dir = dirname(rstudioapi::getSourceEditorContext()$path)

###################################################################################################
#Part 1: Importing the data and cleaning it (data preprocessing)
###################################################################################################

#Extract and split the data
data = readRDS(paste0(root.dir, "/DepMap19Q1_allData.RDS")) #read in the data (with paste0 you can concatinate strings; here some paths this makes your code very flexible)
mut <- data$mutation #pick the mutation data (this has a specific format and neets to be treated seperatelly)

'%!in%' <- function(x,y)!('%in%'(x,y)) #define an operator that will only pick the data that is NOT defined in the list; so the data that needs to be excluded so in our case: DNase, Methylation and RNA data should be excluded; so it will extract notin
dt_new <- lapply(which(names(data) %!in% "mutation"), function(a) data[[a]]) #extract only non-mutation data
names(dt_new) <- names(data)[which(names(data) %!in% "mutation")] #rename the data with the original names

sample_case = c("Skin Cancer") #define the samples you want to compare

#Get the data for the specific celltype
samples = data$annotation$DepMap_ID[which(data$annotation$Primary.Disease == sample_case)]

processed_data <- lapply(1:length(dt_new), function(a) { #pick the data for your sample 
  dat_picker <- dt_new[[a]] #pick one file at each iteration 
  output <- dat_picker[,which(colnames(dat_picker) %in% samples)]
  output <- output[complete.cases(output),]
  output <- output[order(rownames(output)),]#make a datacleanup
  return(output)
})
names(processed_data) <- names(dt_new)

processed_data$annotation <- data$annotation[which(data$annotation$DepMap_ID %in% samples),]

ids = which(names(mut) %in% samples) #extract the mutation data (take care; this is a list of lists; therefore special treatment needed)
allDepMap_mutation_SkinCancer = lapply(ids, function(a) {
  mut[[a]]})
# putting the names of the matrixes in tne allDep_mutation 
names(allDepMap_mutation_SkinCancer) <- samples

# losing the mutations which are not deleterious
allDepMap_mutation_SkinCancer = lapply(1:34, function(a) {allDepMap_mutation_SkinCancer[[a]][which(allDepMap_mutation_SkinCancer[[a]][,"isDeleterious"]== TRUE), ]})

# losing all the genes which are not in every data frame ######################

#first we need to pick all Gene names we have out of our data
 a <- unique(c(rownames(processed_data[[1]]),rownames(processed_data[[2]]),rownames(processed_data[[3]]),rownames(processed_data[[4]])))

# now we pick these genes which are in all 4 dataframes we need for further analysis
i <- 1
out <- vector("character", length(seq_along(1:16970)))
for (x in seq_along(a)) {
  if(a[x] %in% rownames(processed_data$expression) & a[x] %in% rownames(processed_data$copynumber) & a[x] %in% rownames(processed_data$kd.ceres) & a[x] %in% rownames(processed_data$kd.prob))
  {out[i] <- a[x]
  i <- i+1
  } 
}

processed_data2 <- processed_data

processed_data2 <- lapply(processed_data, function(a) {
  a <- a[which(rownames(a) %in% out),]
  return(a)
})
processed_data2$annotation <- processed_data$annotation
processed_data <- processed_data2

allDepMap_mutation_SkinCancer <- lapply(allDepMap_mutation_SkinCancer, function(a) {
  a <- a[which(as.character(a[,2]) %in% out),]
  return(a)
  })



###################################################################################################
#Part 2: Visualize the data
###################################################################################################
#Prepare the data for plotting

finalPlottingData <- lapply(1:(length(processed_data)-1), function(a) { #take care to take not all the processed data becasue we will not need annotation
  dtPicker <- processed_data[[a]]
  out <- melt(dtPicker) #bind the data togehter that we have samples and values as columns
  out$Gene <- rep(rownames(dtPicker), ncol(dtPicker)) #add the genes; probably this might be useful in a later stage
  out$Case <- names(processed_data)[1:(length(processed_data)-1)][a] #add a labelling column
  colnames(out) <- c("Sample", "Value", "Gene", "Case") #rename the columns
  return(out)
})
names(finalPlottingData) <- names(processed_data)[1:(length(processed_data)-1)] #rename the data 
lapply(finalPlottingData, head)

#Make some nice boxplots
ggplot(data = finalPlottingData$expression, aes(x=Sample, y=Value)) +
  geom_boxplot(aes(fill = Sample), outlier.size = 0.1, outlier.alpha = 0.2) + #reconstruct the outliers a bit (so reduce them in size; because we are interested in the boxplots and not the outliers)
  theme_bw(base_size = 7) + #format the size of the theme nicely
  theme(legend.position= "none", #define the legend position (here no leghend will be needed)
        legend.direction="horizontal", #define the legend direction if one is there
        plot.title = element_text(hjust = 0.5), #make the title of the plot into the middle
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #define the orientation of the text on the x-axis
        legend.title= element_blank(), #no title of the legend should be plotted
        axis.title.x = element_blank(), #no title of the x-axis is relevant; because that would be samples and that is cleare due to the naming
        strip.text.y = element_text(angle = 90)) #define the orientation of the text of the y-axis

#Other things can be done in the same format

#Make a heatmap of the kdCERES values
pheatmap(as.matrix(processed_data$kd.ceres[1:10, 1:34]), clustering_method = "ward.D2",border_color = "white", fontsize = 10, 
         main = paste0("kdCERES for potential 2nd site targets"),
         show_rownames = F, show_colnames = T,
         cutree_rows = 4,
         cutree_cols = 2, 
         fontsize_row=10) #take care blue means really, really significant p-value!

#Make a hyrachial clustering

cor.mat = cor(processed_data$kd.ceres[1:50,], method = "spearman")
cor.dist = as.dist(1 - cor.mat)
cor.hc = hclust(cor.dist, method = "ward.D2")
cor.hc = as.dendrogram(cor.hc)
plot(cor.hc, las = 2, cex.lab = 0.7)
 


#Make a histogram
singleGenes <- as.vector(unique(as.data.frame(rbindlist(lapply(seq_along(allDepMap_mutation_SkinCancer), function(a) {out <- as.data.frame(as.vector(unique(allDepMap_mutation_SkinCancer[[a]]$Hugo_Symbol)))}))))[,1])


geneCounts <- sapply(seq_along(singleGenes), function(a) {
  genePicker <- singleGenes[a] #pick one gene
  print(paste0("I am doing: ", a))
  sumGene <- lapply(seq_along(allDepMap_mutation_SkinCancer), function(b) {
    mutPicker <- allDepMap_mutation_SkinCancer[[b]] #pick one of the 34 mutation lists
    out <- as.data.frame(length(which(mutPicker$Hugo_Symbol == genePicker))) #look how often an entry is is in the mutation list
    return(out)
  })
  geneCount <- colSums(as.data.frame(rbindlist(sumGene))) #sum it up to get the total count for each gene
  return(geneCount)
})
names(geneCounts) <- singleGenes #rename it nicely
geneCounts <- as.data.frame(geneCounts) #make a nice dataframe
colnames(geneCounts) <- c("Value")
sortedGenCounts <- geneCounts[order(-geneCounts$Value), , drop = FALSE] #sort the data frame
head(sortedGenCounts) #be amazed

#Plot only the top 50 genes
plotData <- sortedGenCounts[1:40, ,drop = FALSE]

plotData$Gene <- rownames(plotData)

plotData$Gene <- factor(plotData$Gene, levels = plotData$Value)



ggplot(data = plotData) +
  (geom_bar(mapping = aes(x = Gene, y = Value), stat = "identity")) +
  theme_bw(base_size = 7) + #format the size of the theme nicely
  theme(legend.position= "none", #define the legend position (here no leghend will be needed)
        legend.direction="horizontal", #define the legend direction if one is there
        plot.title = element_text(hjust = 0.5), #make the title of the plot into the middle
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #define the orientation of the text on the x-axis
        legend.title= element_blank(), #no title of the legend should be plotted
        axis.title.x = element_blank(), #no title of the x-axis is relevant; because that would be samples and that is cleare due to the naming
        strip.text.y = element_text(angle = 90)) #define the orientation of the text of the y-axis


###################################################################################################
#Part 3: Making some kind of dimensionality reduction (PCA and kMeans)
###################################################################################################
#Extact the data for the top 40
topDriverGenes <- rownames(plotData)

dataTopDriverGenes <- lapply(1:length(processed_data), function(a) { #pick the data for your sample 
  dat_picker <- processed_data[[a]] #pick one file at each iteration 
  output <- dat_picker[which(rownames(dat_picker) %in% topDriverGenes),]
  return(output)
})
names(dataTopDriverGenes) <- names(dt_new)
lapply(dataTopDriverGenes, head)

#Extrac the data for epression and ceres
driverexpression <- dataTopDriverGenes$expression
driverkd.ceres <- dataTopDriverGenes$kd.ceres

#targetexpression ist driverexpression --> könnt ihr selbst weiter machen
#ich weiß nicht so recht, worauf ihr da alles selektiert habt bei euch; also ich selektiere hier einfach auf die TOP 40 eurer Driver mutations....


##########################################kmeans#####################################################
lapply(1:34, function(v){ "ABAT" %in% allDepMap_mutation_SkinCancer[[5]][,"Hugo_Symbol"]})



km = kmeans(x =t(driverexpression), centers = 6, nstart = 30)
wss = sapply(1:11, function(k) {
  kmeans(x = t(driverexpression), centers = k)$tot.withinss
})
plot(1:11, wss, type = "b", pch = 19, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares")


D = dist(t(driverexpression))
km = kmeans(x = t(driverexpression), centers = 8, nstart = 100)

fviz_cluster(km, data = t(driverexpression))

D = dist(t(driverexpression))
plot(silhouette(km$cluster, D))
mean(silhouette(km$cluster, D))
a = silhouette(km$cluster, D)


# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(t(driverexpression), centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(t(driverexpression)))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")




fviz_nbclust(t(driverexpression), kmeans, method = "silhouette")

km2 <- kmeans(t(driverexpression), centers = 2, nstart = 25)
km3 <- kmeans(t(driverexpression), centers = 6, nstart = 25)
km4 <- kmeans(t(driverexpression), centers = 4, nstart = 25)
km5 <- kmeans(t(driverexpression), centers = 9, nstart = 25)

# plots to compare
p1 <- fviz_cluster(km2, geom = "text", labelsize = 8, data = t(driverexpression)) + ggtitle("k = 2")
p2 <- fviz_cluster(km3, geom = "text", labelsize = 8, data = t(driverexpression)) + ggtitle("k = 6")
p3 <- fviz_cluster(km4, geom = "text", labelsize = 8, data = t(driverexpression)) + ggtitle("k = 4")
p4 <- fviz_cluster(km5, geom = "text", labelsize = 8, data = t(driverexpression)) + ggtitle("k = 9")


grid.arrange(p1, p2, p3, p4, nrow = 2)

rm(km,km2,km3,km4,km5,p1,p2,p3,p4)

######PCA #####################################################################################

pca = prcomp(t(driverexpression), center = F, scale. = F)
summary(pca)
autoplot(prcomp(t(driverexpression)), data = t(driverexpression), colour = 'blue')
#zum Anzeigen von labels (Zelllinien)
autoplot(prcomp(t(driverexpression)), data = t(driverexpression), colour = 'blue', label = TRUE, label.size = 2)
str(pca)
plot(pca, type = "l")

###################################################################################################
#Part 4: Making a statistical test
###################################################################################################
###driverdata <- load("/Users/davidschwarzenbacher/Downloads/driverdata.RDS") #load your data (this should be redundant)

driverGenes <- topDriverGenes[1:10] #only use the TOP 10 driver genes
ttestgenes <- rownames(driverkd.ceres)

potSecondSites <- lapply(seq_along(driverGenes), function(a) {
  genePicker <- driverGenes[a] #pick one driver gene
  print(paste0("I am doing driver mut: ", a))
  output <- sapply(seq_along(rownames(processed_data$kd.ceres)), function(b) { #the kdCERES matrix is of interest take its' rownames as refrence
    secondSitePicker <- rownames(processed_data$kd.ceres)[b] #pick a potetnial 2nd site target
    if (secondSitePicker != genePicker) {
      drMUT <- processed_data$kd.ceres[which(rownames(processed_data$kd.ceres) == genePicker),] #pick the driver mut data
      sndMUT <- as.vector(processed_data$kd.ceres[which(rownames(processed_data$kd.ceres) == secondSitePicker),]) #pick the 2nd site data
      cor.val <- cor.test(unlist(drMUT, use.names=FALSE) , unlist(sndMUT, use.names=FALSE), method = "spearman") #make a spearman correlation
      return(cor.val$p.value) #return the p-values
    } else {
      return(1)
    }
  })
  names(output) <- rownames(processed_data$kd.ceres) #rename all
  output <- as.data.frame(output) #get a nice dataframe
  return(output)
})
names(potSecondSites) <- driverGenes #rename the list of lists
lapply(potSecondSites, head) #look at the nice data

# order the data according to their p-values
potSecondSites <- lapply(potSecondSites, function(a){
  a <- as.data.frame(cbind(a$output, rownames(a)))
  a <- a[order(a[1]), ]
})


# selecting the 20 Genes out of every DriverGene List with the lowest p score
potSecondSitestop20 <- lapply(seq_along(potSecondSites), function (a){
  output <- potSecondSites[[a]][1:20,]
  return(output)
})
names(potSecondSitestop20) <- driverGenes


potSecondSitestop20plot <- lapply(seq_along(potSecondSitestop20), function(a){
  output <- ggplot(data = potSecondSitestop20[[a]]) +
    (geom_bar(mapping = aes(x = V2, y = V1), stat = "identity")) +
    theme_bw(base_size = 7) + #format the size of the theme nicely
    theme(legend.position= "none", #define the legend position (here no leghend will be needed)
          legend.direction="horizontal", #define the legend direction if one is there
          plot.title = element_text(hjust = 0.5), #make the title of the plot into the middle
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #define the orientation of the text on the x-axis
          legend.title= element_blank(), #no title of the legend should be plotted
          axis.title.x = element_blank(), #no title of the x-axis is relevant; because that would be samples and that is cleare due to the naming
          strip.text.y = element_text(angle = 90)) #define the orientation of the text of the y-axis
  return(output)
})





###################################################################################################
#Part 5: Making a multiple linear regression (MLR)
###################################################################################################



sessionInfo() #finally done:)