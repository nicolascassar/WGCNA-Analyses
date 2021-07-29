library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("tidyverse")
library("viridis")
library("reshape")
library("WGCNA")


# extract SAV table
OTUs = as(otu_table(data18S.log), "matrix")

# transpose 
OTUs <- t(OTUs)

datExpr0 = as.data.frame(OTUs)

# check data for excessive missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# cluster the samples to check for outliners
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12, 9)
par(cex = 0.6);
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Load trait data
traitData <- read.table ("allstn18S.txt", header=T, sep="\t", na.strings = "NaN")

Samples = rownames(datExpr0)
traitRows = match(Samples, traitData$SampleID)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage()

# Visualize how env traits related to clustered samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits[10:23], signed = FALSE);

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits[10:23]),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr0, datTraits, file = "LTER18S-datainput.RData")

# Set up 
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# This line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Co-expression similarity and adjacency
softPower = 4;
adjacency = adjacency(datExpr0, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab = "", sub = "", mean = "ASV clustering on TOM-based dissimilarity",
     labels = F, hang = 0.04)

# Set as 10 to generate medium to large modules
minModuleSize = 10

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram with module colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "ASV dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits[10:23], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits[10:23]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$vncp);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

module = "turquoise"
column = match(module, modNames);
moduleGenes = dynamicColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for vncp",
                   main = paste("Module membership vs. taxa significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

annot = read.table("Combined_taxonomy.txt", header=T, check.names=F, na.strings="", sep="\t")
probes = names(datExpr0)
probes2annot = match(probes, annot$ASVID)

# check for the number of probes without annotation (should return zero)
sum(is.na(probes2annot))

# Create the starting data frame
TaxaInfo0 = data.frame(Taxon = probes,
                       TaxaSymbol = annot$ASVID[probes2annot],
                       LinkID = annot$Taxonomy[probes2annot],
                       moduleColor = dynamicColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for weight (vncp)
modOrder = order(-abs(cor(MEs, weight, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(TaxaInfo0)
  TaxaInfo0 = data.frame(TaxaInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(TaxaInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the OTUs in the geneInfo variable first by module color, then by geneTraitSignificance
TaxaOrder = order(TaxaInfo0$moduleColor, -abs(TaxaInfo0$GS.weight));
TaxaInfo = TaxaInfo0[TaxaOrder, ]
