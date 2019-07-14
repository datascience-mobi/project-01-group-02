

#Importing the libs
library(ggplot2)
library(relaimpo)
library(factoextra)
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
data = readRDS(paste0(root.dir, "/DepMap19Q1_allData.RDS")) #read in the data (with paste0 you can concatinate strings; here some paths this makes the code very flexible)
mut <- data$mutation #pick the mutation data (this has a specific format and neets to be treated seperatelly)

'%!in%' <- function(x,y)!('%in%'(x,y)) #define an operator that will only pick the data that is NOT defined in the list; so the data that needs to be excluded 
dt_new <- lapply(which(names(data) %!in% "mutation"), function(a) data[[a]]) #extract only non-mutation data
names(dt_new) <- names(data)[which(names(data) %!in% "mutation")] #rename the data with the original names

sample_case = c("Skin Cancer") #defining at which samples we will take out of the  original dataset

#Get the data for Skin Cancer Cell lines
samples = data$annotation$DepMap_ID[which(data$annotation$Primary.Disease == sample_case)]
# we are looking at the column which tells us which ist telling us to which cancer type this cell line belongs and only taking the 
# Cell lines for which the sample case ist true 
# we are getting a vector of all the cell line names we want to look at 

processed_data <- lapply(1:length(dt_new), function(a) { #pick the data for your sample 
  dat_picker <- dt_new[[a]] #pick one file at each iteration 
  if(names(dt_new[a])== "annotation"){ # treating the annotations differnetly because the cell line names are in a colum and are not the columnames like in the other matrices
    output <- dat_picker[which(dat_picker[,1] %in% samples),]
  } else {
  output <- dat_picker[,which(colnames(dat_picker) %in% samples)]# only taking the skin cancer cell lines 
  output <- output[complete.cases(output),] # only taking rows without NAs 
  output <- output[order(rownames(output)),] # reordering the Genes acording to their name
  }
  return(output)
})
names(processed_data) <- names(dt_new) # rename the objects according to the original data
rm(dt_new,sample_case) # remove objects we don`t need anymore

#extract the mutation data 
ids = which(names(mut) %in% samples) 
allDepMap_mutation_SkinCancer = lapply(ids, function(a) {
  mut[[a]]})

rm(mut, ids, data) #tidying

# losing the mutations which are not deleterious
allDepMap_mutation_SkinCancer = lapply(1:34, function(a) {
  allDepMap_mutation_SkinCancer[[a]][which(allDepMap_mutation_SkinCancer[[a]][,"isDeleterious"]== TRUE), ]
    })
names(allDepMap_mutation_SkinCancer) <- samples


# losing all the genes which are not in every data frame ######################

#first we need to pick all Gene names we have out of our data
Genenames <- unique(c(rownames(processed_data[[1]]),rownames(processed_data[[2]]),rownames(processed_data[[3]]),rownames(processed_data[[4]])))

# now we pick these genes which are in all 4 dataframes we need for further analysis
i <- 1
out <- vector("character", length(seq_along(1:16970)))
for (x in seq_along(Genenames)) {
  if(Genenames[x] %in% rownames(processed_data$expression) & Genenames[x] %in% rownames(processed_data$copynumber) & Genenames[x] %in% rownames(processed_data$kd.ceres) & Genenames[x] %in% rownames(processed_data$kd.prob))
  {out[i] <- Genenames[x]
  i <- i+1
  } 
}

allDepMap_annotation_SkinCancer <- processed_data$annotation # saving the annotation object in a seperate dataframe
# because it doesnt contain any information about the genes 

processed_data <- lapply(processed_data[1:4], function(a) {
  a <- a[which(rownames(a) %in% out),]
  return(a)
})

processed_data$mutation <- allDepMap_mutation_SkinCancer
processed_data$annotation <- allDepMap_annotation_SkinCancer
rm(i,out, Genenames,x, allDepMap_annotation_SkinCancer, samples, allDepMap_mutation_SkinCancer)



###################### general Plotting Data with all the data ###################################################


generalPlottingData <- lapply(1:(length(processed_data)-2), function(a) { #take care to take not all the processed data becasue we will not need annotation
  dtPicker <- processed_data[[a]]
  out <- melt(dtPicker) #bind the data togehter that we have samples and values as columns
  out$Gene <- rep(rownames(dtPicker), ncol(dtPicker)) #add the genes; probably this might be useful in a later stage
  out$Case <- names(processed_data)[1:(length(processed_data)-1)][a] #add a labelling column
  colnames(out) <- c("Sample", "Value", "Gene", "Case") #rename the columns
  return(out)
})
names(generalPlottingData) <- names(processed_data)[1:(length(processed_data)-2)] #rename the data 




############################ Plotting Data - Driver Mutations ###########################################

#producing a vector which encompases every gene which at least mutated once 
singleGenes <- as.vector(unique(as.data.frame(rbindlist(lapply(seq_along(processed_data$mutation), function(a) {
  out <- as.data.frame(as.vector(unique(processed_data$mutation[[a]]$Hugo_Symbol)))}))))[,1])

# creating a dataframe which contains how often every gene is mutated 
geneCounts <- sapply(seq_along(singleGenes), function(a) {
  genePicker <- singleGenes[a] #pick one gene
  sumGene <- lapply(seq_along(processed_data$mutation), function(b) {
    mutPicker <- processed_data$mutation[[b]] #pick one of the 34 mutation lists
    out <- as.data.frame(length(which(mutPicker$Hugo_Symbol == genePicker))) #look how often an entry is is in the mutation list
    return(out)
  })
  geneCount <- colSums(as.data.frame(rbindlist(sumGene))) #sum it up to get the total count for each gene
  return(geneCount)
})
names(geneCounts) <- singleGenes #rename 
geneCounts <- as.data.frame(geneCounts) #make a nice dataframe
colnames(geneCounts) <- c("Value")
geneCounts <- geneCounts[order(-geneCounts$Value), , drop = FALSE] #sort the data frame
head(geneCounts) 

#Extact the data for the top 10 which will be our Drivermutations in our further investigation

dataTopDriverGenes <- lapply(1:(length(processed_data)-2), function(a) { #pick the data for your sample 
  dat_picker <- processed_data[[a]] #pick one file at each iteration 
  output <- dat_picker[which(rownames(dat_picker) %in% rownames(geneCounts)[1:10]),] # compare the rownames of the picked data with the names of the 10 most mutated genes
  return(output)
})
names(dataTopDriverGenes) <- names(processed_data)[1:4]

rm(singleGenes)


#################### extracting the drivermutations for every Celline ##################################

#put all Mutationdata in one Matrix 
OneMatrix <- data.frame()
for (i in c(1:34)) {
  OneMatrix <- rbind(OneMatrix,processed_data$mutation[[i]][,Hugo_Symbol:DepMap_ID])
}

#extract just the column of the Gene name and the Zellline
ZelllinesMutations <- OneMatrix[which(OneMatrix$Hugo_Symbol %in% rownames(geneCounts)[1:10] ),]
ZelllinesMutations <- cbind(ZelllinesMutations$Hugo_Symbol, ZelllinesMutations$DepMap_ID)



# extracting the drivermuations for every Zellline out of the dataframe and putting it into another dataframe so it can be used for the plotting

Genes <- c("COL11A1,TMTC2,TTN", " HMCN1", "COL11A1,HMCN1,SLC510", "HMCN1,TMTC2", "COL11A1,TP53,TTN","none","ZNF292","RYR2","HMCN" ,"none2","none3", "TP53, TTN","HMCN1", "TTN,ZNF292","TMTC2,TP53,NEB","TP53", "TMTC2,NEB","none4","TMTC2,TTN,ZNF292", "none5","CACNA1I","HMCN1,TP53,ZNF292","none6","none7","HMCN1,TMTC2,ZNF292","RYR2,TMTC2,NEB","RYR2,NEB,TTN,CACNA1I","HMCM1,TP53","TTN","COL11A1,SLC5A10","COL11A1,CACNA1I","TTN,CACNA1I","RYR2,CACNA1I,ZNF292","TP53,TTN,CACNA1I" )
Zelllines <- c(colnames(processed_data$expression))
zellinesMutations <- as.data.frame(cbind(Zelllines, Genes))

rm(OneMatrix, Genes, ZelllinesMutations, Zelllines,i)

# why we do these kind of data extraction we will outline in the following data Visualization part 

###################################################################################################
#Part 2: Visualize the data
###################################################################################################


############################ HEATMAP with the knock down data ################################################


pheatmap(as.matrix(processed_data$kd.ceres[1:50,]), clustering_method = "ward.D2",border_color = "white", fontsize = 10, 
         main = paste0("kdCERES for potential 2nd site targets"),
         show_rownames = F, show_colnames = T,
         cutree_rows = 4,
         cutree_cols = 2, 
         fontsize_row=10) 

# we can see that there are clear differences between the knockdown data depending on what gene you is kncked down in a specific cell 
# so if there are different behaviours between the cells when knocking down the same Genes then there may be a pattern in 
# This matrix consist of gene knockdown scores. The score is a measure of how essential/important is a particular gene for the cell survival. This score reflects 
# whether upon knocking down that genes does the cell reduce its proliferation or increases it or has no change.
# Smaller values refers to higher essentiality.

############################ Distribution of the Expression values between the different cell lines ##################################################


# Here we see how the the expression values of the genes are distributed in every Cell line 
# So as you can see there are many genes between the 25 and 75 quantile but there are also many outliers 
# These outliers will be form special interest in the following data analysis 
# For now we can just say that the data is differnetly distributed between the celllines 
# based on different mutations in the different cell lines 

ggplot(data = generalPlottingData$expression, aes(x=Sample, y=Value)) +
  geom_boxplot(aes(fill = Sample), outlier.size = 0.1, outlier.alpha = 0.2) + #reconstruct the outliers a bit (so reduce them in size; because we are interested in the boxplots and not the outliers)
  theme_bw(base_size = 7) + #format the size of the theme nicely
  theme(legend.position= "none", #define the legend position (here no leghend will be needed)
        legend.direction="horizontal", #define the legend direction if one is there
        plot.title = element_text(hjust = 0.5), #make the title of the plot into the middle
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #define the orientation of the text on the x-axis
        legend.title= element_blank(), #no title of the legend should be plotted
        axis.title.x = element_blank(), #no title of the x-axis is relevant; because that would be samples and that is cleare due to the naming
        strip.text.y = element_text(angle = 90)) #define the orientation of the text of the y-axis

######################## Plotting how often a gene is mutated over all Cell lines #############################


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

# So we have seen that the Expression values shows a different distribution for the different cell lines 
# different Expressionvalues can arise from gene mutaitons of specific genes
# so the question is if there are Mutations which occur more often than others 
# we suspect that these may be one of the reasons for the differing Expressionvalues
# So as we can see: yes there are Mutations which occure significantly more often



# Now we want to see which Mutations are the top 10 mutated Genes 
# These 10 Genes will be our Driver Genes for which we want to characterize interactions with other genes 

plotData <- geneCounts[1:10, ,drop = FALSE]

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
rm(plotData)

###################################################################################################
#Part 3: Making some kind of dimensionality reduction (PCA and kMeans)
###################################################################################################

# general Questions in this part: 
#   - can we group the different Driver mutations together so that we can see in which  other genes the 
#     Cell lines with a specific driver mutation differentiate 
#   - through that we could gain insight which other genes are our secound targets 

###########################  hirachial clustering #################################################




#drivergene <- 3 # determines which of the Drivermutations will be seen in the cluster at the x axis
dataset <- processed_data$expression # determines which dataset we use

#colnames(dataset)[which(colnames(dataset) %in% unique(zellinesMutations[which(zellinesMutations[,1] == rownames(geneCounts)[drivergene]),2]))] <- rownames(geneCounts)[drivergene]
colnames(dataset) <- zellinesMutations$Genes

cor.mat = cor(dataset[1:50,], method = "spearman")
cor.dist = as.dist(1 - cor.mat)
cor.hc = hclust(cor.dist, method = "ward.D2")
cor.hc = as.dendrogram(cor.hc)
plot(cor.hc, las = 2, cex.lab = 2, main = "Clustering of the expression values of all Zelllines")
#wie bekomme ich es hin das man die ganzen Namen sieht ?


rm(drivergene, realcelllinenames, dataset, cor.hc, cor.mat, cor.dist)


##########################################kmeans#####################################################

dataset <- t(processed_data$expression[-which(rownames(processed_data$expression) %in% rownames(geneCounts)[1:10]),]) 
# determines which dataset we use
# we are trying to cluster the cell lines with the same drivermutations in the same cluster according to the 
# expression data without the expression of the Drivermutaitons
# because we want to see what is driving the differences betweent the cell lines except for the Drivermutation expression values

rownames(dataset) <- zellinesMutations$Genes

dataset <- dataset[,-which(apply(dataset, 2, function(x) {
  var(x)
}) == 0)]


# kick method for choosing the "right" count of centers for the clustering
wss = sapply(1:9, function(k) {
  kmeans(x = dataset, centers = k)$tot.withinss
})
plot(1:9, wss, type = "b", pch = 19, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares")
# theres no knick in this curve so we need to use other methods to tell us how much centers to choose 

# so we are trying to choose by using the silhoutte method 

fviz_nbclust(dataset, kmeans, method = "silhouette")
# the clustering with two centers seems to be the best by far according to the 

# were taking a look at the silhouette values for the clustering with two centers and 
km = kmeans(x =dataset, centers = 2, nstart = 100)
plot(silhouette(km$cluster,dist(dataset)), main = "Silhouette Values for 2 clusters")

km2 <- kmeans(dataset, centers = 2, nstart = 100)
km3 <- kmeans(dataset, centers = 5, nstart = 100)
km4 <- kmeans(dataset, centers = 4, nstart = 100)
km5 <- kmeans(dataset, centers = 10, nstart = 100)

p1 <- fviz_cluster(km2,geom = "text", labelsize = 8, data = dataset) + ggtitle("k = 2")
p2 <- fviz_cluster(km3, geom = "text", labelsize = 8, data = dataset) + ggtitle("k = 5")
p3 <- fviz_cluster(km4, geom = "text", labelsize = 8, data = dataset) + ggtitle("k = 4")
p4 <- fviz_cluster(km5, geom = "text", labelsize = 8, data = dataset) + ggtitle("k = 10")

grid.arrange(p1, p2, p3, p4, nrow = 2)

plot(p4) # clearly the clustering with 10 centers does not conclude in clusters with the same Drivermutations
# the reason for that may be that most of our cell lines have more than one Driver mutation

plot(p1) # the clustering with 2 centers seems to be the best one 
# our next step in the pca will be to see which of the genes drive the differentation of the celllines
# in this plot because they will be the most variable and thus interessting ones 

rm(km,km2,km3,km4,km5,p1,p2,p3,p4, dataset,wss)

############################ PCA #####################################################################################

# now after we saw how the data is clustering together we want to see what is driving the differences 
# for that we are looking at the first two Principal Components 

#drivergene <- 4 # determines which of the Drivermutations will be seen in the cluster at the x axis
dataset <- processed_data$expression # determines which dataset we use


#colnames(dataset)[which(colnames(dataset) %in% unique(ZelllinesMutations[which(ZelllinesMutations[,1] == topDriverGenes[drivergene]),2]))] <- topDriverGenes[drivergene]
colnames(dataset)<- zellinesMutations$Genes

pca = prcomp(t(dataset), center = F, scale. = F)
summary(pca)
#zum Anzeigen von labels (Zelllinien)


fviz_eig(pca)
str(pca)
autoplot(pca, colour = 'blue')
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
# we see once again two clusters
# the first principal component contains the most information about the data 

var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}

loadings <- pca$rotation
sdev <- pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 

var.cos2 <- var.coord^2
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}

var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[, 1:4])
top100var.contrib <- var.contrib[,1]
top100var.contrib <- as.data.frame(top100var.contrib[order(-top100var.contrib)])
top100var.contrib$Genes <- rownames(top100var.contrib)
top100var.contrib <- top100var.contrib[1:100,]
colnames(top100var.contrib)[1] <- "Contribution"


ggplot(data = top100var.contrib) +
  (geom_bar(mapping = aes(x = Genes, y = Contribution), stat = "identity")) +
  theme_bw(base_size = 7) + #format the size of the theme nicely
  theme(legend.position= "none", #define the legend position (here no leghend will be needed)
        legend.direction="horizontal", #define the legend direction if one is there
        plot.title = element_text(hjust = 0.5), #make the title of the plot into the middle
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), #define the orientation of the text on the x-axis
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), #define the orientation of the text on the x-axis
        legend.title= element_blank(), #no title of the legend should be plotted
        axis.title.x = element_blank(), #no title of the x-axis is relevant; because that would be samples and that is cleare due to the naming
        strip.text.y = element_text(angle = 90)) #define the orientation of the text of the y-axis


# these are the Components which are contributing the most to our variation in the data 
# may be we will find some of these in our result of the p-test 

rm(top100var.contrib, drivergene, realcelllinenames, dataset, loadings, pca, realcelllinenames, var.contrib, var.coord, var.cos2, comp.cos2, sdev)
   



###################################################################################################
#Part 4: Making a statistical test
###################################################################################################
###driverdata <- load("/Users/davidschwarzenbacher/Downloads/driverdata.RDS") #load your data (this should be redundant)

driverGenes <- rownames(geneCounts)[1:10] #only use the TOP 10 driver genes
ttestgenes <- rownames(processed_data$kd.ceres)

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

ggplot(data = potSecondSitestop20$TTN) +
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

rm(potSecondSites, ttestgenes)


###################################################################################################
#Part 5: Making a multiple linear regression (MLR)
###################################################################################################

### Predicting the expression of our Drivergenes with all the data #############

#First creating the dataframe for the multiple linear Regression with all the dataframes as Columns and the rows are every gene in every Cell line
a <- generalPlottingData$expression[,1:3]
a <-a[,c(1,3,2)]
copynumber <- generalPlottingData$copynumber[,2]
kd.ceres <- generalPlottingData$kd.ceres[,2]
kd.prob <- generalPlottingData$kd.prob[,2]

RegData <- cbind(a,copynumber,kd.ceres,kd.prob)

# doing the multiple linear Regression 
# then comparint the predicted values of our model with the real values 
#of the test_data by spearman correlaiton
# doing this for every Driver Gene

Regressionanalysis <-lapply(1:10, function(x){
  RegData <- cbind(a,copynumber,kd.ceres,kd.prob)
  Driverexpression <- c()
  for (i in 1:34) {
    a <- 16970*i
    c <- (16970* (i-1))+1
    b <- colnames(processed_data$expression)[i]
    Driverexpression[c:a] <- processed_data$expression[rownames(geneCounts)[x],b]
  }
  print(paste0("I am doing driver mut: ", rownames(geneCounts)[x]))
  RegData <- cbind(RegData,Driverexpression)
  RegData <-as.data.frame(RegData)
  colnames(RegData) <- as.vector(colnames(RegData))
  set.seed(123) #initialize the random numbers
  split = sample.split(RegData, SplitRatio = 0.8) #split the dataset into 4/5 Training and 1/5 Testing dataset
  training_set = subset(RegData, split == TRUE) #use the labels to get the training data
  test_set = subset(RegData, split == FALSE) 
  rm(RegData)
  # Fitting Multiple Linear Regression to the Training set
  regressor = lm(Driverexpression ~ Value + copynumber + kd.ceres + kd.prob , data = training_set) #predict profit based on ALL (=.) the input variables for one company 
  return(regressor)
  # Predicting the Test set results
  y_pred = predict(regressor, newdata = test_set, se.fit = TRUE) #predict the expression based on your testing data 
  test_set$Prediction = y_pred$fit #add your predictions to the dataset
  #Now compare the Predictions (last column) with the real values of the startups (2nd last column)
  Results <- cor.test(test_set$Driverexpression, test_set$Prediction, method = "spearman", exact=FALSE)
  return(Results)
})
names(Regressionanalysis) <- rownames(geneCounts)[1:10]
Regressionanalysis <- as.vector(Regressionanalysis)
rm(RegData,kd.ceres,kd.prob,copynumber,a)




######## Predicting the Expression of the Drivermutations with the expression of the other genes ########

dataset <- t(processed_data$expression)

#dataset <- dataset[-which(apply(dataset, 1, function(x) {
#  var(x)
#}) == 0),]
dataset <- as.data.frame(dataset)



set.seed(123) #initialize the random numbers
split = sample.split(dataset, SplitRatio = 0.5) #split the dataset into 4/5 Training and 1/5 Testing dataset
training_set = subset(dataset, split == TRUE) #use the labels to get the training data
test_set = subset(dataset, split == FALSE) 
rm(dataset)


# Fitting Multiple Linear Regression to the Training set
regressor = lm(TP53 ~ ., data = training_set) #predict profit based on ALL (=.) the input variables for one company 

# Predicting the Test set results
y_pred = predict(regressor, newdata = test_set, se.fit = TRUE) #predict the expression based on your testing data 
test_set$Prediction = y_pred$fit #add your predictions to the dataset
#Now compare the Predictions (last column) with the real values of the startups (2nd last column)
Results <- cor.test(test_set$Driverexpression, test_set$Prediction, method = "spearman", exact=FALSE)







sessionInfo() #finally done:) 
