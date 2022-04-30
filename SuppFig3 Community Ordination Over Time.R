library(phyloseq)
library("ggplot2"); packageVersion("ggplot2")

#################### Supplementary Figure 3 ###########################
#Plotting Directly From Source Data File
BrayPCoAData <- read.table(file="/Users/tommyauchtung/Desktop/SuppFig3.txt", header = TRUE, sep = "\t")

#Bacteria
ggplot(BrayPCoAData, aes(BactBrayPCoA1,BactBrayPCoA2))+geom_point(size=20, aes(colour = Month)) + 
  theme_bw() +
  theme(axis.text= element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_gradientn(colors=rainbow(5))

#Fungi
ggplot(BrayPCoAData, aes(FungiBrayPCoA1,FungiBrayPCoA2))+geom_point(size=20, aes(colour = Month)) + 
  theme_bw() +
  theme(axis.text= element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_gradientn(colors=rainbow(5))

#Plants
ggplot(BrayPCoAData, aes(PlantBrayPCoA1,PlantBrayPCoA2))+geom_point(size=20, aes(colour = Month)) + 
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_gradientn(colors=rainbow(5)) 

#Export images as 2000 x 2000 tiffs

#################### HOW DATA WAS GENERATED ###############################################
#############################################################################################################################################################################################################
BACTERIA
#############################################################################################################################################################################################################
#Merge Samples by Month
biom.FilePath <- "/Users/tommyauchtung/Desktop/16S_OTU_Table.r3K.biom"
meta.FilePath <- "/Users/tommyauchtung/Desktop/agemap.txt"

phylo <- import_biom( biom.FilePath, 
                      parallel = TRUE,
                      parseFunction = function(x){
                        x <- gsub("^[\ _]+", "", x)
                        x <- gsub("[\ ]+$",  "", x)
                        names(x) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:length(x)]
                        return(x)
                      })
MappingFile <- read.table(meta.FilePath, header=TRUE, comment.char='', sep="\t",
                          check.names=TRUE, strip.white=TRUE, blank.lines.skip=TRUE)
rownames(MappingFile) <- MappingFile[['SampleID']]
phylo2 <- phyloseq(otu_table(phylo), sample_data(MappingFile))

merged_table2 <- merge_samples(phylo2, "AgeinMonths")
write.csv(otu_table(merged_table2),file="/Users/tommyauchtung/Desktop/bacteriasamplesmerged.csv")
#Use Excel and Macqiime to convert into biom file

#############################################################################################################################################################################################################
#Determine Bacterial Ordination Coordinates
biom.FilePath <- "/Users/tommyauchtung/Desktop/BacteriaR3Kmerged.biom"

# Read in the OTU table
phylo <- import_biom( biom.FilePath, 
                      parallel = TRUE,
                      parseFunction = function(x){
                        x <- gsub("^[\ _]+", "", x)
                        x <- gsub("[\ ]+$",  "", x)
                        names(x) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:length(x)]
                        return(x)
                      })

# Rarefy to month with the lowest number of bacterial reads
phylorare <- rarefy_even_depth(phylo, 51000, replace=FALSE)

#For PCoA:
DistanceMatrixBray <- distance(phylorare, method="bray", binary=FALSE)
OrdinationBrayPCoA     <- ordinate(phylorare, "PCoA", distance=DistanceMatrixBray)
CoordsBrayPCoA <- plot_ordination(phylorare, OrdinationBrayPCoA, axes=c(1:2), justDF=TRUE)
write.table(CoordsBrayPCoA, "/Users/tommyauchtung/Desktop/CoordsBrayPCoA16S.txt")
#Clean up file and save as BacteriaCommunityOrdinationCoordinates.txt

########################################################################################################################
#Plotting Bacterial Ordinations
bacteriadata <- read.table(file="/Users/tommyauchtung/Desktop/BacteriaCommunityOrdinationCoordinates.txt", header = TRUE, sep = "\t")

ggplot(bacteriadata, aes(BrayPCoA1,BrayPCoA2))+geom_point(size=20, aes(colour = Month)) + 
  theme_bw() +
  theme(axis.text= element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_gradientn(colors=rainbow(5))

#############################################################################################################################################################################################################
FUNGI
#############################################################################################################################################################################################################
#Merge Samples by Month
biom.FilePath <- "/Users/tommyauchtung/Desktop/OTU_Table.Merged.99.justFungi.NoBleed.r3K.biom"
meta.FilePath <- "/Users/tommyauchtung/Desktop/agemap.txt"

phylo <- import_biom( biom.FilePath, 
                      parallel = TRUE,
                      parseFunction = function(x){
                        x <- gsub("^[\ _]+", "", x)
                        x <- gsub("[\ ]+$",  "", x)
                        names(x) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:length(x)]
                        return(x)
                      })
MappingFile <- read.table(meta.FilePath, header=TRUE, comment.char='', sep="\t",
                          check.names=TRUE, strip.white=TRUE, blank.lines.skip=TRUE)
rownames(MappingFile) <- MappingFile[['SampleID']]
phylo2 <- phyloseq(otu_table(phylo), sample_data(MappingFile))

merged_table2 <- merge_samples(phylo2, "AgeinMonths")
write.csv(otu_table(merged_table2),file="/Users/tommyauchtung/Desktop/fungisamplesmerged.csv")
#Use Excel and Macqiime to convert into biom file

#############################################################################################################################################################################################################
#Determine Fungal Ordination Coordinates
biom.FilePath <- "/Users/tommyauchtung/Desktop/TEDDY_OTU_Table.Merged.99.justFungi.NoBleed.r3K.CollapsedByAge.biom"
# Read in the OTU table
phylo <- import_biom( biom.FilePath, 
                      parallel = TRUE,
                      parseFunction = function(x){
                        x <- gsub("^[\ _]+", "", x)
                        x <- gsub("[\ ]+$",  "", x)
                        names(x) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:length(x)]
                        return(x)
                      })

# Rarefy to month with the lowest number of fungal reads
phylorare <- rarefy_even_depth(phylo, 45000, replace=FALSE)

#For PCoA:
DistanceMatrixBray <- distance(phylorare, method="bray", binary=FALSE)
OrdinationBrayPCoA     <- ordinate(phylorare, "PCoA", distance=DistanceMatrixBray)
CoordsBrayPCoA <- plot_ordination(phylorare, OrdinationBrayPCoA, axes=c(1:2), justDF=TRUE)
write.table(CoordsBrayPCoA, "/Users/tommyauchtung/Desktop/CoordsBrayPCoAFungi.txt")
#Clean up file and save as FungiCommunityOrdinationCoordinates.txt

########################################################################################################################
#Plotting Fungal Ordinations
fungidata <- read.table(file="/Users/tommyauchtung/Desktop/FungiCommunityOrdinationCoordinates.txt", header = TRUE, sep = "\t")

ggplot(fungidata, aes(BrayPCoA1,BrayPCoA2))+geom_point(size=20, aes(colour = Month)) + 
  theme_bw() +
  theme(axis.text= element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_gradientn(colors=rainbow(5))

#############################################################################################################################################################################################################
PLANTS
#############################################################################################################################################################################################################
#Merge Samples by Month
biom.FilePath <- "/Users/tommyauchtung/Desktop/OTU_Table.Merged.99.justPlants.NoBleed.r100.biom"
meta.FilePath <- "/Users/tommyauchtung/Desktop/agemap.txt"

phylo <- import_biom( biom.FilePath, 
                      parallel = TRUE,
                      parseFunction = function(x){
                        x <- gsub("^[\ _]+", "", x)
                        x <- gsub("[\ ]+$",  "", x)
                        names(x) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:length(x)]
                        return(x)
                      })
MappingFile <- read.table(meta.FilePath, header=TRUE, comment.char='', sep="\t",
                          check.names=TRUE, strip.white=TRUE, blank.lines.skip=TRUE)
rownames(MappingFile) <- MappingFile[['SampleID']]
phylo2 <- phyloseq(otu_table(phylo), sample_data(MappingFile))

merged_table2 <- merge_samples(phylo2, "AgeinMonths")
write.csv(otu_table(merged_table2),file="/Users/tommyauchtung/Desktop/plantsamplesmerged.csv")
#Use Excel and Macqiime to convert into biom file

#############################################################################################################################################################################################################
#Determine Plant Ordination Coordinates
biom.FilePath <- "/Users/tommyauchtung/Desktop/MergedPlant.biom"

# Read in the OTU table
phylo <- import_biom( biom.FilePath, 
                      parallel = TRUE,
                      parseFunction = function(x){
                        x <- gsub("^[\ _]+", "", x)
                        x <- gsub("[\ ]+$",  "", x)
                        names(x) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:length(x)]
                        return(x)
                      })

# Rarefy to month with the lowest number of plant reads [40m had 2600, 43 had 1900, 48 had just 1100]
phylorare <- rarefy_even_depth(phylo, 2600, replace=FALSE)

#For PCoA:
DistanceMatrixBray <- distance(phylorare, method="bray", binary=FALSE)
OrdinationBrayPCoA     <- ordinate(phylorare, "PCoA", distance=DistanceMatrixBray)
CoordsBrayPCoA <- plot_ordination(phylorare, OrdinationBrayPCoA, axes=c(1:2), justDF=TRUE)
write.table(CoordsBrayPCoA, "/Users/tommyauchtung/Desktop/CoordsBrayPCoAPlants.txt")
#Clean up file and save as PlantCommunityOrdinationCoordinates.txt

########################################################################################################################
#Plotting Plant Ordinations
plantdata <- read.table(file="/Users/tommyauchtung/Desktop/PlantCommunityOrdinationCoordinates.txt", header = TRUE, sep = "\t")

ggplot(plantdata, aes(BrayPCoA1R2600,BrayPCoA2R2600))+geom_point(size=20, aes(colour = Month)) + 
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_gradientn(colors=rainbow(5)) 

