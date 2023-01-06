# Set the script's title and author
title <- "Gene Expression Analysis"
author <- "JosephAgiM"

# Get the current date
date <- Sys.time()

# Loading necessary packages
library(BiocManager)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(pathview)

# Read in data table and store it in a variable 
Roels_GeneName <- read.table("AnnotedgeR_normcounts.tabular", header = T) 
# Read in data table and store it in a variable 
Roels_GeneID <- read.table("edgeR_normcounts.tabular", header = T) 
# Read in data table and store it in a variable 
Verboom_GeneName <- read.table("AnnotedgeR_Verboom.tabular", header = T)
# Read in data table and store it in a variable 
Verboom_GeneID <- read.table("edgeR_Verboom.tabular", header = T)

# Replace the GeneID column with the SYMBOL column in Roels_GeneID data table
Roels_GeneID$GeneID <- Roels_GeneName$SYMBOL
# Replace the GeneID column with the SYMBOL column in Verboom_GeneID data table
Verboom_GeneID$GeneID <- Verboom_GeneName$SYMBOL


# Remove NA values from Roels_GeneID data table using the na.omit() function
Roels_GeneID <- na.omit(Roels_GeneID)
# Remove NA values from Verboom_GeneID data table using the na.omit() function
Verboom_GeneID <- na.omit(Verboom_GeneID)
# Merge the first columns of both data tables and add the other columns next to each other 
merged_GeneID <- merge(Roels_GeneID, Verboom_GeneID, c(1))

# Remove the logFC, logCPM, F, P-value, and FDR columns from the merged data table
merged_GeneID <- merged_GeneID[-c(10:14)]
# View the first few rows of the Roels_GeneID data table
head(Roels_GeneID)
# View the first few rows of the Verboom_GeneID data table
head(Verboom_GeneID)

# Move the GeneID column to the rownames of the merged data table
rownames(merged_GeneID) <- merged_GeneID$GeneID
# Remove the GeneID column from the merged data table
merged_GeneID <- subset(merged_GeneID, select = -GeneID)
# View the first few rows of the merged data table
head(merged_GeneID)

# Rename the first four columns of the merged data table
colnames(merged_GeneID)[c(1:4)] <- c("Verboom1","Verboom2","Verboom3","Verboom4")
# Rename the last four columns of the merged data table
colnames(merged_GeneID)[c(5:8)] <- c("Roels1","Roels2","Roels3","Roels4")

# Creating new column with the means of Roels and Verboom data frames 
merged_GeneID$RoelsMean <- rowMeans(merged_GeneID[,5:8])
merged_GeneID$VerboomMean <- rowMeans(merged_GeneID[,1:4])


# Getting p-value column for the GeneID row in a new data frame P_Value  
Pvalue <-  data.frame(Pvalue=rep(0,nrow(merged_GeneID)))

# Assigning the row names of Pvalue to be the same as those of merged_GeneID
rownames(Pvalue) <- rownames(merged_GeneID)

# Looping through the t.test on each row (i) for each sample type (Roels and Verboom)
# and storing the resulting p-values in Pvalue$Pvalue
for (i in 1:nrow(merged_GeneID)) {
  Pvalue$Pvalue[i] <- t.test(merged_GeneID[i,1:4], merged_GeneID[i,5:8], alternative = "two.sided")$p.value
}

# Removing rows with P-values over 0.05 and storing the remaining rows in signPvalue
signPvalue <- Pvalue[Pvalue$Pvalue < 0.05, , drop= F]

# Adjusting p-values using the Benjamini-Hochberg method and storing the results in a new data frame
AdjustedPvalue <- data.frame(AdjustedPvalue = p.adjust(Pvalue[,"Pvalue"], method = "BH"))

# Assigning the row names of AdjustedPvalue to be the same as those of Pvalue
rownames(AdjustedPvalue) <- rownames(Pvalue)

# Removing rows with adjusted P-values over 0.05 and storing the remaining rows in AdjustedPvalue
AdjustedPvalue <- AdjustedPvalue[AdjustedPvalue$AdjustedPvalue < 0.05,  ,drop = F] 

# Creating a vector of the significant genes based on the significant p-values
signGenes <- rownames(AdjustedPvalue)

# Subsetting merged_GeneID based on the significant genes and storing the result in signGenesTable
signGenesTable <- subset(merged_GeneID, rownames(merged_GeneID) %in% signGenes)

# Converting signGenesTable to a matrix
signGenesMatrix <- data.matrix(signGenesTable)

# Visualizing signGenesMatrix using a heatmap
heatmap(signGenesMatrix)


# Visualizing signGenesMatrix using a prettier heatmap
pheatmap(signGenesMatrix)

# Adjusting p-values using the Benjamini-Hochberg method and storing the results in a new data frame
AdjustedPvalueOFClog <- data.frame(AdjustedPvalue = p.adjust(Pvalue[,"Pvalue"], method = "BH"))

# Assigning the row names of AdjustedPvalueOFClog to be the same as those of Pvalue
rownames(AdjustedPvalueOFClog) <- rownames(Pvalue)

# Calculating the fold change (FClog) and storing it in AdjustedPvalueOFClog$FClog
AdjustedPvalueOFClog$FClog <- rowMeans(merged_GeneID[,1:4]) - rowMeans(merged_GeneID[,5:8])

#viewing too dataframe
head(AdjustedPvalueOFClog)

# Adding the first eight columns of merged_GeneID to AdjustedPvalueOFClog
AdjustedPvalueOFClog[,3:10] <- merged_GeneID[,1:8]

# Filtering AdjustedPvalueOFClog to only include rows with adjusted P-values less than 0.05
signPvalue1 <- AdjustedPvalueOFClog[AdjustedPvalueOFClog$AdjustedPvalue < 0.05, , drop = F]

# Re-ordering signPvalue1 based on FClog
signPvalue1 <- signPvalue1[order(signPvalue1$FClog), , drop = F]

# Selecting the top and bottom 10 rows of signPvalue1
topBottomGenes <- signPvalue1[c(1:10, 4229:4238),1:10]

# Removing the first two columns of topBottomGenes
topBottomGenes_tabl <- topBottomGenes[-c(1:2)]

# Converting topBottomGenes_tabl to a matrix
topBottomGenes_Matr <- data.matrix(topBottomGenes_tabl)

# Visualizing topBottomGenes_Matr using a heatmap
pheatmap(topBottomGenes_Matr)

# Preparing the row names (gene names) of topBottomGenes_tabl as a column called GeneNames
topBottomGenes_tabl$GeneNames <- rownames(topBottomGenes_tabl)

# Calculate the mean of the first four columns (belonging to the Verboom sample set) 
# and store the result in topBottomGenes_tabl$VerboomMean
topBottomGenes_tabl$VerboomMean <- rowMeans(topBottomGenes_tabl[,1:4]) 

# Calculate the mean of the last four columns (belonging to the Roels sample set) 
# and store the result in topBottomGenes_tabl$RoelsMean
topBottomGenes_tabl$RoelsMean <- rowMeans(topBottomGenes_tabl[,5:8])

# Plot the means of the Verboom and Roels sample sets using a scatterplot
ggplot(topBottomGenes_tabl, aes(x=VerboomMean, y=RoelsMean, label=GeneNames))+
  geom_point()+  # add points to the plot
  geom_abline(intercept = 0)+  # add a line with intercept 0
  theme_bw()+  # use a black and white theme
  ggtitle(label="Top 10 up and downregulated genes")+  # add a title to the plot
  aes(color = VerboomMean)+  # color the points based on VerboomMean
  geom_text_repel()+  # add labels to the points, avoiding overlap
  theme(plot.title = element_text(hjust = 0.5, size=20))  # adjust the title size and alignment

# Remove eight columns of AdjustedPvalueOFClog and store the result of 2 first columns in AdjustedPvalueOFClog_1
AdjustedPvalueOFClog_1 <- AdjustedPvalueOFClog[-c(3:10)]

# Standard volcano plot
Vulkan <- AdjustedPvalueOFClog_1 %>%
  ggplot(mapping = aes(x=FClog, y=-log10(AdjustedPvalue)))+
  geom_point()+
  theme_minimal()


# Add a line to highlight the padjusted value of 0.05 and a threshold for the logFC. 

Vulkan2 <- Vulkan +
  geom_hline(yintercept = -log10(0.05), col="red") +
  geom_vline(xintercept= c(-2, 2), col="red")

Vulkan2

# Color the dots based on their differential expression
AdjustedPvalueOFClog_1$diffexpressed <- "NO" # create a new column in the data frame and fill it with NO
AdjustedPvalueOFClog_1$diffexpressed[AdjustedPvalueOFClog_1$FClog > 2 & AdjustedPvalueOFClog_1$AdjustedPvalue < 0.05] <- "UP"
AdjustedPvalueOFClog_1$diffexpressed[AdjustedPvalueOFClog_1$FClog < -2 & AdjustedPvalueOFClog_1$AdjustedPvalue < 0.05] <- "DOWN"

Vulkan <- AdjustedPvalueOFClog_1 %>%
  ggplot(mapping=aes(x=FClog, y=-log10(AdjustedPvalue), col=diffexpressed))+
  geom_point()+
  theme_minimal()

Vulkan2 <- Vulkan +
  geom_hline(yintercept = -log10(0.05), col="red")+
  geom_vline(xintercept= c(-2,2), col="red")

mycolors <- c("red", "green", "black") # specify the colors you want to use

names(mycolors) <- c("DOWN", "UP", "NO") # add column names

Vulkan3 <- Vulkan2+
  scale_color_manual(values=mycolors)

Vulkan3


# Viewing data
dim(merged_GeneID) # number of rows and columns    
head(merged_GeneID) # console logs dataframe


# Plot two genes only
TwoGenes <- merged_GeneID[c("ATF3", "NDNF"),]
TwoGenes$genenames <- rownames(TwoGenes)

# Get my data into a long format

TwoGenes_boxPlot <-gather(TwoGenes, key="samples", value="values", -c(VerboomMean, RoelsMean, genenames))
TwoGenes_boxPlot$group <- rep(c("Roels", "Verboom"), each=8)


TwoGenes_boxPlot <- stack(TwoGenes[,1:8])
TwoGenes_boxPlot$genenames <- rep(c("ATF3", "NDNF"), times=8)
TwoGenes_boxPlot$group <- rep(c("Roels", "Verboom"), each=8)



# Use the long data for plotting a boxplot

ggplot(TwoGenes_boxPlot, aes(x=genenames, y=values, fill=group))+
  geom_boxplot()+
  scale_fill_manual(values=c("green", "red"))+
  theme_light()


# Fixing dataframe
# Data table pasted to object
keggObject <- read.table("edgeR_Verboom.tabular", header = T)
# Differential Expressed Genes
# Keep only DEGs using a threshold of adjusted P < 0.05
degs <- keggObject[keggObject$FDR < 0.05, ]

# Creating a character object for the backround genes
keggGene <- keggObject$GeneID
keggGene <- as.character(keggGene)

# KEGG pathway enrichment for top 1000 DEGs
kegg_enrich <- enrichKEGG(gene = degs$GeneID[1:1000], 
                          organism = "hsa", # Homo sapiens
                          keyType = "kegg", # KEGG database
                          pvalueCutoff = 0.05, 
                          pAdjustMethod = "BH", 
                          universe = keggGene) # Background genes
View(kegg_enrich@result)

# Plotting the top 10 significant KEGG pathways 
dotplot(kegg_enrich, x = "Count", showCategory = 10, orderBy = "Count", title = "Top significant KEGG pathways")


# Create a list of logFCs with their corresponding genes as names
FClog <- keggObject$logFC
names(FClog) <- keggObject$GeneID

# Run pathview for a particular pathway ID that was the highest represented in the Kegg dotplot
pathview(gene.data = FClog, 
         pathway.id = "hsa04151", 
         low = "red", mid = "grey", high = "green",
         species = "hsa", 
         limit = c(-10, 10))


