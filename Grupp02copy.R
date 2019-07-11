

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
library(caTools)

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
rm(processed_data2)



######################Plotting Data ###################################################


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



###################plotting Data - mutation ###########################################

singleGenes <- as.vector(unique(as.data.frame(rbindlist(lapply(seq_along(allDepMap_mutation_SkinCancer), function(a) {out <- as.data.frame(as.vector(unique(allDepMap_mutation_SkinCancer[[a]]$Hugo_Symbol)))}))))[,1])


geneCounts <- sapply(seq_along(singleGenes), function(a) {
  genePicker <- singleGenes[a] #pick one gene
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

#Extact the data for the top 40

plotData <- sortedGenCounts[1:10, ,drop = FALSE]
topDriverGenes <- rownames(plotData)

dataTopDriverGenes <- lapply(1:length(processed_data), function(a) { #pick the data for your sample 
  dat_picker <- processed_data[[a]] #pick one file at each iteration 
  output <- dat_picker[which(rownames(dat_picker) %in% topDriverGenes),]
  return(output)
})
names(dataTopDriverGenes) <- names(dt_new)


#Extrac the data for epression and ceres
driverexpression <- dataTopDriverGenes$expression
driverkd.ceres <- dataTopDriverGenes$kd.ceres




# putting the names of the matrixes in one allDep_mutation 
names(allDepMap_mutation_SkinCancer) <- samples

OneMatrix <- data.frame()
for (i in c(1:34)) {
  OneMatrix <- rbind(OneMatrix, allDepMap_mutation_SkinCancer[[i]][,Hugo_Symbol:DepMap_ID])
}

ZelllinesMutations <- OneMatrix[which(OneMatrix$Hugo_Symbol %in% topDriverGenes ),]
ZelllinesMutations <- cbind(ZelllinesMutations$Hugo_Symbol, ZelllinesMutations$DepMap_ID)


###################################################################################################
#Part 2: Visualize the data
###################################################################################################


#Make a heatmap of the kdCERES values
pheatmap(as.matrix(processed_data$kd.ceres[25:45, 1:34]), clustering_method = "ward.D2",border_color = "white", fontsize = 10, 
         main = paste0("kdCERES for potential 2nd site targets"),
         show_rownames = F, show_colnames = T,
         cutree_rows = 4,
         cutree_cols = 2, 
         fontsize_row=10) #take care blue means really, really significant p-value!


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



#Plotting the spread in the number of Mutations per Gene in a Boxplot
geneCounts <- cbind(geneCounts, "Mutations per Gene")

ggplot(data = geneCounts, aes(x="Mutations per Gene", y=Value)) +
  geom_boxplot(aes(fill = "Mutations per Gene"), outlier.size = 2, outlier.alpha = 0.2) + #reconstruct the outliers a bit (so reduce them in size; because we are interested in the boxplots and not the outliers)
  theme_bw(base_size = 7) + #format the size of the theme nicely
  theme(legend.position= "none", #define the legend position (here no leghend will be needed)
        legend.direction="horizontal", #define the legend direction if one is there
        plot.title = element_text(hjust = 0.5), #make the title of the plot into the middle
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust= 0.5, size = 10), #define the orientation of the text on the x-axis
        legend.title= element_blank(), #no title of the legend should be plotted
        axis.title.x = element_blank(), #no title of the x-axis is relevant; because that would be samples and that is cleare due to the naming
        strip.text.y = element_text(angle = 90)) #define the orientation of the text of the y-axis

#Plot only the top 50 genes


plotData <- sortedGenCounts[1:10, ,drop = FALSE]

plotData$Gene <- rownames(plotData)

ggplot(data = plotData) +
  (geom_bar(mapping = aes(x = Gene, y = Value), stat = "identity")) +
  theme_bw(base_size = 7) + #format the size of the theme nicely
  theme(legend.position= "none", #define the legend position (here no leghend will be needed)
        legend.direction="horizontal", #define the legend direction if one is there
        plot.title = element_text(hjust = 0.5), #make the title of the plot into the middle
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), #define the orientation of the text on the x-axis
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), #define the orientation of the text on the x-axis
        legend.title= element_blank(), #no title of the legend should be plotted
        axis.title.x = element_blank(), #no title of the x-axis is relevant; because that would be samples and that is cleare due to the naming
        strip.text.y = element_text(angle = 90)) #define the orientation of the text of the y-axis


###################################################################################################
#Part 3: Making some kind of dimensionality reduction (PCA and kMeans)
###################################################################################################


###########################  hirachial clustering #################################################



drivergene <- 4 # determines which of the Drivermutations will be seen in the cluster at the x axis
dataset <- processed_data$kd.ceres # determines which dataset we use
realcelllinenames <- processed_data$copynumber # is just defined so we can switch the colnames back so that they are showing the celllines after the clustering

colnames(dataset)[which(colnames(dataset) %in% unique(ZelllinesMutations[which(ZelllinesMutations[,1] == topDriverGenes[drivergene]),2]))] <- topDriverGenes[drivergene]

cor.mat = cor(dataset[1:50,], method = "spearman")
cor.dist = as.dist(1 - cor.mat)
cor.hc = hclust(cor.dist, method = "ward.D2")
cor.hc = as.dendrogram(cor.hc)
plot(cor.hc, las = 2, cex.lab = 0.7)

colnames(dataset) <- realcelllinenames
rm(drivergene, realcelllinenames, dataset)


View(ZelllinesMutations)



##########################################kmeans#####################################################



km = kmeans(x =t(driverexpression), centers = 6, nstart = 30)
wss = sapply(1:11, function(k) {
  kmeans(x = t(driverexpression), centers = k)$tot.withinss
})
plot(1:11, wss, type = "b", pch = 19, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares")

D = dist(t(driverexpression))
plot(silhouette(km$cluster, D))

# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(t(driverexpression), centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(t(driverexpression)))
  mean(ss[, 3])
}

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

drivergene <- 4 # determines which of the Drivermutations will be seen in the cluster at the x axis
dataset <- processed_data$expression # determines which dataset we use
realcelllinenames <- processed_data$copynumber # is just defined so we can switch the colnames back so that they are showing the celllines after the clustering

colnames(dataset)[which(colnames(dataset) %in% unique(ZelllinesMutations[which(ZelllinesMutations[,1] == topDriverGenes[drivergene]),2]))] <- topDriverGenes[drivergene]


pca = prcomp(t(dataset), center = F, scale. = F)
summary(pca)
#zum Anzeigen von labels (Zelllinien)
autoplot(prcomp(t(dataset)), data = t(dataset), colour = 'blue', label = TRUE, label.size = 3)
str(pca)
plot(pca, type = "l")

colnames(dataset) <- realcelllinenames
rm(drivergene, realcelllinenames, dataset)

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

ggplot(data = melt(potSecondSitestop20)) +
  (geom_bar(mapping = aes(x = V2, y = V1), stat = "identity")) +
  theme_bw(base_size = 7) + #format the size of the theme nicely
  theme(legend.position= "none", #define the legend position (here no leghend will be needed)
        legend.direction="horizontal", #define the legend direction if one is there
        plot.title = element_text(hjust = 0.5), #make the title of the plot into the middle
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), #define the orientation of the text on the x-axis
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), #define the orientation of the text on the x-axis
        legend.title= element_blank(), #no title of the legend should be plotted
        axis.title.x = element_blank(), #no title of the x-axis is relevant; because that would be samples and that is cleare due to the naming
        strip.text.y = element_text(angle = 90)) #define the orientation of the text of the y-axis


###################################################################################################
#Part 5: Making a multiple linear regression (MLR)
###################################################################################################
rm(a, cor.dist, D, driverGene, i , plotData, mut, finalPlottingData, ids, out, R_GC_MEM_GROW, root.dir, sample_case, samples, singleGenes, split, topDriverGenes, wss, x, allDepMap_mutation_SkinCancer, cor.hc, cor.mat, data, driverexpression, driverkd.ceres, dt_new, geneCounts, OneMatrix, pca, potSecondSites, potSecondSitestop20, sortedGenCounts)

MRegressionanalysis <- t(processed_data$copynumber)[,1:16000]

MRegressionanalysis <-as.data.frame(MRegressionanalysis)
colnames(MRegressionanalysis) <- as.vector(colnames(MRegressionanalysis))





##################################################################################

# Splitting the dataset into the Training set and Test set
set.seed(123) #initialize the random numbers
split = sample.split(MRegressionanalysis[,"A1BG"], SplitRatio = 0.8) #split the dataset into 4/5 Training and 1/5 Testing dataset
training_set = subset(MRegressionanalysis, split == TRUE) #use the labels to get the training data
test_set = subset(MRegressionanalysis, split == FALSE) #dim(test_set) will give you know 10 --> 50/5*1 = 10; wuhu train/test split worked



# Fitting Multiple Linear Regression to the Training set
regressor = lm(A1BG ~. , data = training_set) #predict profit based on ALL (=.) the input variables for one company 




# Predicting the Test set results
y_pred = predict(regressor, newdata = test_set, se.fit = TRUE) #predict the profit based on your testing data (this data the model did NEVER see and highly usefull to evaluat the performance)
test_set$Prediction = y_pred #add your predictions to the dataset
test_set #Now you can compare your Predictions (last column) with the real values of the startups (2nd last column)

cor.test(test_set$A1BG, test_set$Prediction, method = "spearman")

y_pred$fit

sessionInfo() #finally done:) 
