library(phyloseq)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(plyr)

#################### Final Figure 2b ###########################
BothStats <- read.table(file="/Users/tommyauchtung/Desktop/Fig2b.txt", header = TRUE, sep = "\t")
ggplot(BothStats)  +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line    = element_line(color='black'),
        plot.margin=margin(l=20,b=15,t=15,r=20,unit = "pt")) + 
  labs(x = "Age (Months)") + 
  geom_smooth(aes(x=Time, y=OneMinusMean, color=Category), method="loess", size=1, se=FALSE) +
  geom_errorbar(aes(x=Time, ymin=OneMinuslower, ymax=OneMinusupper,color=Category), width=0.1) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48)) +
  coord_cartesian(xlim = c(0, 48), ylim = c(0.0, 0.7), expand = FALSE) +
  geom_point(aes(x=Time, y=OneMinusMean, shape = BetweenWithin, color=Category), size = 3) +
  scale_shape_manual(values = c(19, 1))+
  scale_color_manual(values = c("red","red","blue","blue")) +
  guides(col = guide_legend(reverse = TRUE))

#Export image as 12 x 6 pdf 

#################### HOW DATA WAS GENERATED ###########################
#### FUNGI #########################################################################################################
biom.FilePath <- "/Users/CMMR/Desktop/OTU_Table.Merged.97.JustFungi.ForExport.biom"

# Read in the OTU table
phylo <- import_biom( biom.FilePath, 
                      parallel = TRUE,
                      parseFunction = function(x){
                        x <- gsub("^[\ _]+", "", x)
                        x <- gsub("[\ ]+$",  "", x)
                        names(x) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:length(x)]
                        return(x)
                      })
# Rarefy
phylo3000 <- rarefy_even_depth(phylo, 3000, replace=FALSE)

# Calculate the distance matrix and save
DistanceMatrixBray<- distance(phylo3000, method="bray", binary=FALSE)
write.table(as.matrix(DistanceMatrixBray), file = "/Users/CMMR/Desktop/FungiDistanceMatrixBray.txt", sep = "\t", row.names = TRUE, col.names = NA)
#Get rid of quotes by: cat FungiDistanceMatrixBray.txt | sed -e 's/"//g' > FungiDistanceMatrixBrayNoQuotes.txt 

# File names to read from / write to
MappingFile <- "/Users/CMMR/Desktop/MP139_SUPERMETADATA.txt"
GroupingCol <- "SubjectMaskId"
DistMatrix  <- "/Users/CMMR/Desktop/FungiDistanceMatrixBrayNoQuotes.txt" 

# Read in the mapping file and distance matrix
md <- read.table(MappingFile, header=TRUE, row.names=1, sep="\t")
dmfull <- as.matrix(read.table(DistMatrix, header=TRUE, row.names=1, as.is=TRUE))

# To subset by age
md2 <- subset(md, ExactAgeinMonths >= 2 & ExactAgeinMonths < 3)
md2[['SampleID']] <- rownames(md2)
md2 <- plyr::ddply(md2, 'SubjectMaskId', head, n=1)
rownames(md2) <- md2[['SampleID']]

# To subset by age
md3 <- subset(md, ExactAgeinMonths >= 3 & ExactAgeinMonths < 4)
md3[['SampleID']] <- rownames(md3)
md3 <- plyr::ddply(md3, 'SubjectMaskId', head, n=1)
rownames(md3) <- md3[['SampleID']]

md4 <- subset(md, ExactAgeinMonths >= 4 & ExactAgeinMonths < 5)
md4[['SampleID']] <- rownames(md4)
md4 <- plyr::ddply(md4, 'SubjectMaskId', head, n=1)
rownames(md4) <- md4[['SampleID']]

md5 <- subset(md, ExactAgeinMonths >= 5 & ExactAgeinMonths < 6)
md5[['SampleID']] <- rownames(md5)
md5 <- plyr::ddply(md5, 'SubjectMaskId', head, n=1)
rownames(md5) <- md5[['SampleID']]

md6 <- subset(md, ExactAgeinMonths >= 6 & ExactAgeinMonths < 7)
md6[['SampleID']] <- rownames(md6)
md6 <- plyr::ddply(md6, 'SubjectMaskId', head, n=1)
rownames(md6) <- md6[['SampleID']]

md7 <- subset(md, ExactAgeinMonths >= 7 & ExactAgeinMonths < 8)
md7[['SampleID']] <- rownames(md7)
md7 <- plyr::ddply(md7, 'SubjectMaskId', head, n=1)
rownames(md7) <- md7[['SampleID']]

md8 <- subset(md, ExactAgeinMonths >= 8 & ExactAgeinMonths < 9)
md8[['SampleID']] <- rownames(md8)
md8 <- plyr::ddply(md8, 'SubjectMaskId', head, n=1)
rownames(md8) <- md8[['SampleID']]

md9 <- subset(md, ExactAgeinMonths >= 9 & ExactAgeinMonths < 10)
md9[['SampleID']] <- rownames(md9)
md9 <- plyr::ddply(md9, 'SubjectMaskId', head, n=1)
rownames(md9) <- md9[['SampleID']]

md10 <- subset(md, ExactAgeinMonths >= 10 & ExactAgeinMonths < 11)
md10[['SampleID']] <- rownames(md10)
md10 <- plyr::ddply(md10, 'SubjectMaskId', head, n=1)
rownames(md10) <- md10[['SampleID']]

md11 <- subset(md, ExactAgeinMonths >= 11 & ExactAgeinMonths < 12)
md11[['SampleID']] <- rownames(md11)
md11 <- plyr::ddply(md11, 'SubjectMaskId', head, n=1)
rownames(md11) <- md11[['SampleID']]

md12 <- subset(md, ExactAgeinMonths >= 12 & ExactAgeinMonths < 13)
md12[['SampleID']] <- rownames(md12)
md12 <- plyr::ddply(md12, 'SubjectMaskId', head, n=1)
rownames(md12) <- md12[['SampleID']]

md13 <- subset(md, ExactAgeinMonths >= 13 & ExactAgeinMonths < 14)
md13[['SampleID']] <- rownames(md13)
md13 <- plyr::ddply(md13, 'SubjectMaskId', head, n=1)
rownames(md13) <- md13[['SampleID']]

md14 <- subset(md, ExactAgeinMonths >= 14 & ExactAgeinMonths < 15)
md14[['SampleID']] <- rownames(md14)
md14 <- plyr::ddply(md14, 'SubjectMaskId', head, n=1)
rownames(md14) <- md14[['SampleID']]

md15 <- subset(md, ExactAgeinMonths >= 15 & ExactAgeinMonths < 16)
md15[['SampleID']] <- rownames(md15)
md15 <- plyr::ddply(md15, 'SubjectMaskId', head, n=1)
rownames(md15) <- md15[['SampleID']]

md16 <- subset(md, ExactAgeinMonths >= 16 & ExactAgeinMonths < 17)
md16[['SampleID']] <- rownames(md16)
md16 <- plyr::ddply(md16, 'SubjectMaskId', head, n=1)
rownames(md16) <- md16[['SampleID']]

md17 <- subset(md, ExactAgeinMonths >= 17 & ExactAgeinMonths < 18)
md17[['SampleID']] <- rownames(md17)
md17 <- plyr::ddply(md17, 'SubjectMaskId', head, n=1)
rownames(md17) <- md17[['SampleID']]

md18 <- subset(md, ExactAgeinMonths >= 18 & ExactAgeinMonths < 19)
md18[['SampleID']] <- rownames(md18)
md18 <- plyr::ddply(md18, 'SubjectMaskId', head, n=1)
rownames(md18) <- md18[['SampleID']]

md19 <- subset(md, ExactAgeinMonths >= 19 & ExactAgeinMonths < 20)
md19[['SampleID']] <- rownames(md19)
md19 <- plyr::ddply(md19, 'SubjectMaskId', head, n=1)
rownames(md19) <- md19[['SampleID']]

md20 <- subset(md, ExactAgeinMonths >= 20 & ExactAgeinMonths < 21)
md20[['SampleID']] <- rownames(md20)
md20 <- plyr::ddply(md20, 'SubjectMaskId', head, n=1)
rownames(md20) <- md20[['SampleID']]

md21 <- subset(md, ExactAgeinMonths >= 21 & ExactAgeinMonths < 22)
md21[['SampleID']] <- rownames(md21)
md21 <- plyr::ddply(md21, 'SubjectMaskId', head, n=1)
rownames(md21) <- md21[['SampleID']]

md22 <- subset(md, ExactAgeinMonths >= 22 & ExactAgeinMonths < 23)
md22[['SampleID']] <- rownames(md22)
md22 <- plyr::ddply(md22, 'SubjectMaskId', head, n=1)
rownames(md22) <- md22[['SampleID']]

md23 <- subset(md, ExactAgeinMonths >= 23 & ExactAgeinMonths < 24)
md23[['SampleID']] <- rownames(md23)
md23 <- plyr::ddply(md23, 'SubjectMaskId', head, n=1)
rownames(md23) <- md23[['SampleID']]

md24 <- subset(md, ExactAgeinMonths >= 24 & ExactAgeinMonths < 25)
md24[['SampleID']] <- rownames(md24)
md24 <- plyr::ddply(md24, 'SubjectMaskId', head, n=1)
rownames(md24) <- md24[['SampleID']]

md25 <- subset(md, ExactAgeinMonths >= 25 & ExactAgeinMonths < 26)
md25[['SampleID']] <- rownames(md25)
md25 <- plyr::ddply(md25, 'SubjectMaskId', head, n=1)
rownames(md25) <- md25[['SampleID']]

md26 <- subset(md, ExactAgeinMonths >= 26 & ExactAgeinMonths < 27)
md26[['SampleID']] <- rownames(md26)
md26 <- plyr::ddply(md26, 'SubjectMaskId', head, n=1)
rownames(md26) <- md26[['SampleID']]

md27 <- subset(md, ExactAgeinMonths >= 27 & ExactAgeinMonths < 28)
md27[['SampleID']] <- rownames(md27)
md27 <- plyr::ddply(md27, 'SubjectMaskId', head, n=1)
rownames(md27) <- md27[['SampleID']]

md28 <- subset(md, ExactAgeinMonths >= 28 & ExactAgeinMonths < 29)
md28[['SampleID']] <- rownames(md28)
md28 <- plyr::ddply(md28, 'SubjectMaskId', head, n=1)
rownames(md28) <- md28[['SampleID']]

md29 <- subset(md, ExactAgeinMonths >= 29 & ExactAgeinMonths < 30)
md29[['SampleID']] <- rownames(md29)
md29 <- plyr::ddply(md29, 'SubjectMaskId', head, n=1)
rownames(md29) <- md29[['SampleID']]

md30 <- subset(md, ExactAgeinMonths >= 30 & ExactAgeinMonths < 31)
md30[['SampleID']] <- rownames(md30)
md30 <- plyr::ddply(md30, 'SubjectMaskId', head, n=1)
rownames(md30) <- md30[['SampleID']]

md31 <- subset(md, ExactAgeinMonths >= 31 & ExactAgeinMonths < 32)
md31[['SampleID']] <- rownames(md31)
md31 <- plyr::ddply(md31, 'SubjectMaskId', head, n=1)
rownames(md31) <- md31[['SampleID']]

md32 <- subset(md, ExactAgeinMonths >= 32 & ExactAgeinMonths < 33)
md32[['SampleID']] <- rownames(md32)
md32 <- plyr::ddply(md32, 'SubjectMaskId', head, n=1)
rownames(md32) <- md32[['SampleID']]

md33 <- subset(md, ExactAgeinMonths >= 33 & ExactAgeinMonths < 34)
md33[['SampleID']] <- rownames(md33)
md33 <- plyr::ddply(md33, 'SubjectMaskId', head, n=1)
rownames(md33) <- md33[['SampleID']]

md34 <- subset(md, ExactAgeinMonths >= 34 & ExactAgeinMonths < 35)
md34[['SampleID']] <- rownames(md34)
md34 <- plyr::ddply(md34, 'SubjectMaskId', head, n=1)
rownames(md34) <- md34[['SampleID']]

md35 <- subset(md, ExactAgeinMonths >= 35 & ExactAgeinMonths < 36)
md35[['SampleID']] <- rownames(md35)
md35 <- plyr::ddply(md35, 'SubjectMaskId', head, n=1)
rownames(md35) <- md35[['SampleID']]

md36 <- subset(md, ExactAgeinMonths >= 36 & ExactAgeinMonths < 37)
md36[['SampleID']] <- rownames(md36)
md36 <- plyr::ddply(md36, 'SubjectMaskId', head, n=1)
rownames(md36) <- md36[['SampleID']]

md37 <- subset(md, ExactAgeinMonths >= 37 & ExactAgeinMonths < 38)
md37[['SampleID']] <- rownames(md37)
md37 <- plyr::ddply(md37, 'SubjectMaskId', head, n=1)
rownames(md37) <- md37[['SampleID']]

md38 <- subset(md, ExactAgeinMonths >= 38 & ExactAgeinMonths < 39)
md38[['SampleID']] <- rownames(md38)
md38 <- plyr::ddply(md38, 'SubjectMaskId', head, n=1)
rownames(md38) <- md38[['SampleID']]

md39 <- subset(md, ExactAgeinMonths >= 39 & ExactAgeinMonths < 40)
md39[['SampleID']] <- rownames(md39)
md39 <- plyr::ddply(md39, 'SubjectMaskId', head, n=1)
rownames(md39) <- md39[['SampleID']]

md40 <- subset(md, ExactAgeinMonths >= 40 & ExactAgeinMonths < 41)
md40[['SampleID']] <- rownames(md40)
md40 <- plyr::ddply(md40, 'SubjectMaskId', head, n=1)
rownames(md40) <- md40[['SampleID']]

md41 <- subset(md, ExactAgeinMonths >= 41 & ExactAgeinMonths < 42)
md41[['SampleID']] <- rownames(md41)
md41 <- plyr::ddply(md41, 'SubjectMaskId', head, n=1)
rownames(md41) <- md41[['SampleID']]

md42 <- subset(md, ExactAgeinMonths >= 42 & ExactAgeinMonths < 43)
md42[['SampleID']] <- rownames(md42)
md42 <- plyr::ddply(md42, 'SubjectMaskId', head, n=1)
rownames(md42) <- md42[['SampleID']]

md43 <- subset(md, ExactAgeinMonths >= 43 & ExactAgeinMonths < 44)
md43[['SampleID']] <- rownames(md43)
md43 <- plyr::ddply(md43, 'SubjectMaskId', head, n=1)
rownames(md43) <- md43[['SampleID']]

md44 <- subset(md, ExactAgeinMonths >= 44 & ExactAgeinMonths < 45)
md44[['SampleID']] <- rownames(md44)
md44 <- plyr::ddply(md44, 'SubjectMaskId', head, n=1)
rownames(md44) <- md44[['SampleID']]

md45 <- subset(md, ExactAgeinMonths >= 45 & ExactAgeinMonths < 46)
md45[['SampleID']] <- rownames(md45)
md45 <- plyr::ddply(md45, 'SubjectMaskId', head, n=1)
rownames(md45) <- md45[['SampleID']]

md46 <- subset(md, ExactAgeinMonths >= 46 & ExactAgeinMonths < 47)
md46[['SampleID']] <- rownames(md46)
md46 <- plyr::ddply(md46, 'SubjectMaskId', head, n=1)
rownames(md46) <- md46[['SampleID']]

md47 <- subset(md, ExactAgeinMonths >= 47 & ExactAgeinMonths < 48)
md47[['SampleID']] <- rownames(md47)
md47 <- plyr::ddply(md47, 'SubjectMaskId', head, n=1)
rownames(md47) <- md47[['SampleID']]

md48 <- subset(md, ExactAgeinMonths >= 48 & ExactAgeinMonths < 49)
md48[['SampleID']] <- rownames(md48)
md48 <- plyr::ddply(md48, 'SubjectMaskId', head, n=1)
rownames(md48) <- md48[['SampleID']]


#Bind Comparisons
md2to3 <- rbind(md2, md3) 
md3to4 <- rbind(md3, md4) 
md4to5 <- rbind(md4, md5) 
md5to6 <- rbind(md5, md6)
md6to7 <- rbind(md6, md7)
md7to8 <- rbind(md7, md8)
md8to9 <- rbind(md8, md9)
md9to10 <- rbind(md9, md10)
md10to11 <- rbind(md10, md11)
md11to12 <- rbind(md11, md12)
md12to13 <- rbind(md12, md13)
md13to14 <- rbind(md13, md14)
md14to15 <- rbind(md14, md15)
md15to16 <- rbind(md15, md16)
md16to17 <- rbind(md16, md17)
md17to18 <- rbind(md17, md18)
md18to19 <- rbind(md18, md19)
md19to20 <- rbind(md19, md20)
md20to21 <- rbind(md20, md21)
md21to22 <- rbind(md21, md22)
md22to23 <- rbind(md22, md23)
md23to24 <- rbind(md23, md24)
md24to25 <- rbind(md24, md25)
md25to26 <- rbind(md25, md26)
md26to27 <- rbind(md26, md27)
md27to28 <- rbind(md27, md28)
md28to29 <- rbind(md28, md29)
md29to30 <- rbind(md29, md30)
md30to31 <- rbind(md30, md31)
md31to32 <- rbind(md31, md32)
md32to33 <- rbind(md32, md33)
md33to34 <- rbind(md33, md34)
md34to35 <- rbind(md34, md35)
md35to36 <- rbind(md35, md36)
md36to37 <- rbind(md36, md37)
md37to38 <- rbind(md37, md38)
md38to39 <- rbind(md38, md39)
md39to40 <- rbind(md39, md40)
md40to41 <- rbind(md40, md41)
md41to42 <- rbind(md41, md42)
md42to43 <- rbind(md42, md43)
md43to44 <- rbind(md43, md44)
md44to45 <- rbind(md44, md45)
md45to46 <- rbind(md45, md46)
md46to47 <- rbind(md46, md47)
md47to48 <- rbind(md47, md48)

###RUN COMPARISONS###
# Set metadata to intended comparison
mf <- md2to3
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"2") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md2to3.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md3to4
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"3") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md3to4.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md4to5
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"4") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md4to5.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md5to6
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"5") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md5to6.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md6to7
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"6") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md6to7.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md7to8
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"7") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md7to8.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md8to9
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"8") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md8to9.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md9to10
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"9") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md9to10.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md10to11
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"10") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md10to11.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md11to12
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"11") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md11to12.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md12to13
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"12") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md12to13.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md13to14
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"13") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md13to14.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md14to15
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"14") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md14to15.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md15to16
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"15") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md15to16.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md16to17
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"16") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md16to17.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md17to18
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"17") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md17to18.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md18to19
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"18") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md18to19.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md19to20
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"19") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md19to20.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md20to21
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"20") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md20to21.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md21to22
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"21") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md21to22.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md22to23
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"22") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md22to23.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md23to24
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"23") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md23to24.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md24to25
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"24") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md24to25.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md25to26
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"25") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md25to26.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md26to27
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"26") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md26to27.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md27to28
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"27") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md27to28.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md28to29
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"28") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md28to29.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md29to30
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"29") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md29to30.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md30to31
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"30") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md30to31.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md31to32
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"31") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md31to32.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md32to33
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"32") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md32to33.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md33to34
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"33") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md33to34.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md34to35
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"34") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md34to35.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md35to36
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"35") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md35to36.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md36to37
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"36") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md36to37.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md37to38
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"37") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md37to38.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md38to39
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"38") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md38to39.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md39to40
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"39") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md39to40.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md40to41
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"40") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md40to41.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md41to42
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"41") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md41to42.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md42to43
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"42") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md42to43.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md43to44
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"43") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md43to44.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md44to45
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"44") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md44to45.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md45to46
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"45") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md45to46.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md46to47
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"46") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md46to47.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md47to48
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"47") 
write.table(dffinal, file = "/Users/CMMR/Desktop/FungiBCdissimilarities/md47to48.txt", sep = "\t", row.names = TRUE, col.names = NA)

#cat all files together into one big file = all.txt
#Get rid of quotes by: cat all.txt | sed -e 's/"//g' > allnoquotes.txt 
#Get rid of titles by: cat allnoquotes.txt | grep -v "Distance" > allready.txt
#Manually add titles back to top line

#Calculating the mean/median/standard deviation of categories so can draw line
BC <- read.table(file="/Users/CMMR/Desktop/FungiBCdissimilarities/allready.txt", header = TRUE, sep = "\t")  
BC$NewCol <- paste(BC$Category, BC$Time, sep='') 
FungiWithinBetweenMeansMediansSDSE <- ddply(BC, .(NewCol), summarize, "Mean"= mean(Distance), "Median" = median(Distance), "SD" = sd(Distance), "SE" = sd(Distance)/sqrt(length(Distance)))
write.table(FungiWithinBetweenMeansMediansSDSE, "/Users/CMMR/Desktop/FungiBCdissimilarities/Stats.txt")
#Manually add back the Time and Category


############## BACTERIA ##########################################################################################
biom.FilePath <- "/Users/CMMR/Desktop/16S_OTU_Table_r3000.biom"

# Read in the OTU table
phylo <- import_biom( biom.FilePath, 
                      parallel = TRUE,
                      parseFunction = function(x){
                        x <- gsub("^[\ _]+", "", x)
                        x <- gsub("[\ ]+$",  "", x)
                        names(x) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:length(x)]
                        return(x)
                      })
# Rarefy (already done)
#phylo3000 <- rarefy_even_depth(phylo, 3000, replace=FALSE)

# Calculate the distance matrix and save
DistanceMatrixBray<- distance(phylo, method="bray", binary=FALSE)
write.table(as.matrix(DistanceMatrixBray), file = "/Users/CMMR/Desktop/16SDistanceMatrixBray.txt", sep = "\t", row.names = TRUE, col.names = NA)
#Get rid of quotes by: cat 16SDistanceMatrixBray.txt | sed -e 's/"//g' > 16SDistanceMatrixBrayNoQuotes.txt 

# File names to read from / write to
MappingFile <- "/Users/CMMR/Desktop/MP139_SUPERMETADATA.txt"
GroupingCol <- "SubjectMaskId"
DistMatrix  <- "/Users/CMMR/Desktop/16SDistanceMatrixBrayNoQuotes.txt" 

# Read in the mapping file and distance matrix
md <- read.table(MappingFile, header=TRUE, row.names=1, sep="\t")
dmfull <- as.matrix(read.table(DistMatrix, header=TRUE, row.names=1, as.is=TRUE))

# To subset by age
md2 <- subset(md, ExactAgeinMonths >= 2 & ExactAgeinMonths < 3)
md2[['SampleID']] <- rownames(md2)
md2 <- plyr::ddply(md2, 'SubjectMaskId', head, n=1)
rownames(md2) <- md2[['SampleID']]

md3 <- subset(md, ExactAgeinMonths >= 3 & ExactAgeinMonths < 4)
md3[['SampleID']] <- rownames(md3)
md3 <- plyr::ddply(md3, 'SubjectMaskId', head, n=1)
rownames(md3) <- md3[['SampleID']]

md4 <- subset(md, ExactAgeinMonths >= 4 & ExactAgeinMonths < 5)
md4[['SampleID']] <- rownames(md4)
md4 <- plyr::ddply(md4, 'SubjectMaskId', head, n=1)
rownames(md4) <- md4[['SampleID']]

md5 <- subset(md, ExactAgeinMonths >= 5 & ExactAgeinMonths < 6)
md5[['SampleID']] <- rownames(md5)
md5 <- plyr::ddply(md5, 'SubjectMaskId', head, n=1)
rownames(md5) <- md5[['SampleID']]

md6 <- subset(md, ExactAgeinMonths >= 6 & ExactAgeinMonths < 7)
md6[['SampleID']] <- rownames(md6)
md6 <- plyr::ddply(md6, 'SubjectMaskId', head, n=1)
rownames(md6) <- md6[['SampleID']]

md7 <- subset(md, ExactAgeinMonths >= 7 & ExactAgeinMonths < 8)
md7[['SampleID']] <- rownames(md7)
md7 <- plyr::ddply(md7, 'SubjectMaskId', head, n=1)
rownames(md7) <- md7[['SampleID']]

md8 <- subset(md, ExactAgeinMonths >= 8 & ExactAgeinMonths < 9)
md8[['SampleID']] <- rownames(md8)
md8 <- plyr::ddply(md8, 'SubjectMaskId', head, n=1)
rownames(md8) <- md8[['SampleID']]

md9 <- subset(md, ExactAgeinMonths >= 9 & ExactAgeinMonths < 10)
md9[['SampleID']] <- rownames(md9)
md9 <- plyr::ddply(md9, 'SubjectMaskId', head, n=1)
rownames(md9) <- md9[['SampleID']]

md10 <- subset(md, ExactAgeinMonths >= 10 & ExactAgeinMonths < 11)
md10[['SampleID']] <- rownames(md10)
md10 <- plyr::ddply(md10, 'SubjectMaskId', head, n=1)
rownames(md10) <- md10[['SampleID']]

md11 <- subset(md, ExactAgeinMonths >= 11 & ExactAgeinMonths < 12)
md11[['SampleID']] <- rownames(md11)
md11 <- plyr::ddply(md11, 'SubjectMaskId', head, n=1)
rownames(md11) <- md11[['SampleID']]

md12 <- subset(md, ExactAgeinMonths >= 12 & ExactAgeinMonths < 13)
md12[['SampleID']] <- rownames(md12)
md12 <- plyr::ddply(md12, 'SubjectMaskId', head, n=1)
rownames(md12) <- md12[['SampleID']]

md13 <- subset(md, ExactAgeinMonths >= 13 & ExactAgeinMonths < 14)
md13[['SampleID']] <- rownames(md13)
md13 <- plyr::ddply(md13, 'SubjectMaskId', head, n=1)
rownames(md13) <- md13[['SampleID']]

md14 <- subset(md, ExactAgeinMonths >= 14 & ExactAgeinMonths < 15)
md14[['SampleID']] <- rownames(md14)
md14 <- plyr::ddply(md14, 'SubjectMaskId', head, n=1)
rownames(md14) <- md14[['SampleID']]

md15 <- subset(md, ExactAgeinMonths >= 15 & ExactAgeinMonths < 16)
md15[['SampleID']] <- rownames(md15)
md15 <- plyr::ddply(md15, 'SubjectMaskId', head, n=1)
rownames(md15) <- md15[['SampleID']]

md16 <- subset(md, ExactAgeinMonths >= 16 & ExactAgeinMonths < 17)
md16[['SampleID']] <- rownames(md16)
md16 <- plyr::ddply(md16, 'SubjectMaskId', head, n=1)
rownames(md16) <- md16[['SampleID']]

md17 <- subset(md, ExactAgeinMonths >= 17 & ExactAgeinMonths < 18)
md17[['SampleID']] <- rownames(md17)
md17 <- plyr::ddply(md17, 'SubjectMaskId', head, n=1)
rownames(md17) <- md17[['SampleID']]

md18 <- subset(md, ExactAgeinMonths >= 18 & ExactAgeinMonths < 19)
md18[['SampleID']] <- rownames(md18)
md18 <- plyr::ddply(md18, 'SubjectMaskId', head, n=1)
rownames(md18) <- md18[['SampleID']]

md19 <- subset(md, ExactAgeinMonths >= 19 & ExactAgeinMonths < 20)
md19[['SampleID']] <- rownames(md19)
md19 <- plyr::ddply(md19, 'SubjectMaskId', head, n=1)
rownames(md19) <- md19[['SampleID']]

md20 <- subset(md, ExactAgeinMonths >= 20 & ExactAgeinMonths < 21)
md20[['SampleID']] <- rownames(md20)
md20 <- plyr::ddply(md20, 'SubjectMaskId', head, n=1)
rownames(md20) <- md20[['SampleID']]

md21 <- subset(md, ExactAgeinMonths >= 21 & ExactAgeinMonths < 22)
md21[['SampleID']] <- rownames(md21)
md21 <- plyr::ddply(md21, 'SubjectMaskId', head, n=1)
rownames(md21) <- md21[['SampleID']]

md22 <- subset(md, ExactAgeinMonths >= 22 & ExactAgeinMonths < 23)
md22[['SampleID']] <- rownames(md22)
md22 <- plyr::ddply(md22, 'SubjectMaskId', head, n=1)
rownames(md22) <- md22[['SampleID']]

md23 <- subset(md, ExactAgeinMonths >= 23 & ExactAgeinMonths < 24)
md23[['SampleID']] <- rownames(md23)
md23 <- plyr::ddply(md23, 'SubjectMaskId', head, n=1)
rownames(md23) <- md23[['SampleID']]

md24 <- subset(md, ExactAgeinMonths >= 24 & ExactAgeinMonths < 25)
md24[['SampleID']] <- rownames(md24)
md24 <- plyr::ddply(md24, 'SubjectMaskId', head, n=1)
rownames(md24) <- md24[['SampleID']]

md25 <- subset(md, ExactAgeinMonths >= 25 & ExactAgeinMonths < 26)
md25[['SampleID']] <- rownames(md25)
md25 <- plyr::ddply(md25, 'SubjectMaskId', head, n=1)
rownames(md25) <- md25[['SampleID']]

md26 <- subset(md, ExactAgeinMonths >= 26 & ExactAgeinMonths < 27)
md26[['SampleID']] <- rownames(md26)
md26 <- plyr::ddply(md26, 'SubjectMaskId', head, n=1)
rownames(md26) <- md26[['SampleID']]

md27 <- subset(md, ExactAgeinMonths >= 27 & ExactAgeinMonths < 28)
md27[['SampleID']] <- rownames(md27)
md27 <- plyr::ddply(md27, 'SubjectMaskId', head, n=1)
rownames(md27) <- md27[['SampleID']]

md28 <- subset(md, ExactAgeinMonths >= 28 & ExactAgeinMonths < 29)
md28[['SampleID']] <- rownames(md28)
md28 <- plyr::ddply(md28, 'SubjectMaskId', head, n=1)
rownames(md28) <- md28[['SampleID']]

md29 <- subset(md, ExactAgeinMonths >= 29 & ExactAgeinMonths < 30)
md29[['SampleID']] <- rownames(md29)
md29 <- plyr::ddply(md29, 'SubjectMaskId', head, n=1)
rownames(md29) <- md29[['SampleID']]

md30 <- subset(md, ExactAgeinMonths >= 30 & ExactAgeinMonths < 31)
md30[['SampleID']] <- rownames(md30)
md30 <- plyr::ddply(md30, 'SubjectMaskId', head, n=1)
rownames(md30) <- md30[['SampleID']]

md31 <- subset(md, ExactAgeinMonths >= 31 & ExactAgeinMonths < 32)
md31[['SampleID']] <- rownames(md31)
md31 <- plyr::ddply(md31, 'SubjectMaskId', head, n=1)
rownames(md31) <- md31[['SampleID']]

md32 <- subset(md, ExactAgeinMonths >= 32 & ExactAgeinMonths < 33)
md32[['SampleID']] <- rownames(md32)
md32 <- plyr::ddply(md32, 'SubjectMaskId', head, n=1)
rownames(md32) <- md32[['SampleID']]

md33 <- subset(md, ExactAgeinMonths >= 33 & ExactAgeinMonths < 34)
md33[['SampleID']] <- rownames(md33)
md33 <- plyr::ddply(md33, 'SubjectMaskId', head, n=1)
rownames(md33) <- md33[['SampleID']]

md34 <- subset(md, ExactAgeinMonths >= 34 & ExactAgeinMonths < 35)
md34[['SampleID']] <- rownames(md34)
md34 <- plyr::ddply(md34, 'SubjectMaskId', head, n=1)
rownames(md34) <- md34[['SampleID']]

md35 <- subset(md, ExactAgeinMonths >= 35 & ExactAgeinMonths < 36)
md35[['SampleID']] <- rownames(md35)
md35 <- plyr::ddply(md35, 'SubjectMaskId', head, n=1)
rownames(md35) <- md35[['SampleID']]

md36 <- subset(md, ExactAgeinMonths >= 36 & ExactAgeinMonths < 37)
md36[['SampleID']] <- rownames(md36)
md36 <- plyr::ddply(md36, 'SubjectMaskId', head, n=1)
rownames(md36) <- md36[['SampleID']]

md37 <- subset(md, ExactAgeinMonths >= 37 & ExactAgeinMonths < 38)
md37[['SampleID']] <- rownames(md37)
md37 <- plyr::ddply(md37, 'SubjectMaskId', head, n=1)
rownames(md37) <- md37[['SampleID']]

md38 <- subset(md, ExactAgeinMonths >= 38 & ExactAgeinMonths < 39)
md38[['SampleID']] <- rownames(md38)
md38 <- plyr::ddply(md38, 'SubjectMaskId', head, n=1)
rownames(md38) <- md38[['SampleID']]

md39 <- subset(md, ExactAgeinMonths >= 39 & ExactAgeinMonths < 40)
md39[['SampleID']] <- rownames(md39)
md39 <- plyr::ddply(md39, 'SubjectMaskId', head, n=1)
rownames(md39) <- md39[['SampleID']]

md40 <- subset(md, ExactAgeinMonths >= 40 & ExactAgeinMonths < 41)
md40[['SampleID']] <- rownames(md40)
md40 <- plyr::ddply(md40, 'SubjectMaskId', head, n=1)
rownames(md40) <- md40[['SampleID']]

md41 <- subset(md, ExactAgeinMonths >= 41 & ExactAgeinMonths < 42)
md41[['SampleID']] <- rownames(md41)
md41 <- plyr::ddply(md41, 'SubjectMaskId', head, n=1)
rownames(md41) <- md41[['SampleID']]

md42 <- subset(md, ExactAgeinMonths >= 42 & ExactAgeinMonths < 43)
md42[['SampleID']] <- rownames(md42)
md42 <- plyr::ddply(md42, 'SubjectMaskId', head, n=1)
rownames(md42) <- md42[['SampleID']]

md43 <- subset(md, ExactAgeinMonths >= 43 & ExactAgeinMonths < 44)
md43[['SampleID']] <- rownames(md43)
md43 <- plyr::ddply(md43, 'SubjectMaskId', head, n=1)
rownames(md43) <- md43[['SampleID']]

md44 <- subset(md, ExactAgeinMonths >= 44 & ExactAgeinMonths < 45)
md44[['SampleID']] <- rownames(md44)
md44 <- plyr::ddply(md44, 'SubjectMaskId', head, n=1)
rownames(md44) <- md44[['SampleID']]

md45 <- subset(md, ExactAgeinMonths >= 45 & ExactAgeinMonths < 46)
md45[['SampleID']] <- rownames(md45)
md45 <- plyr::ddply(md45, 'SubjectMaskId', head, n=1)
rownames(md45) <- md45[['SampleID']]

md46 <- subset(md, ExactAgeinMonths >= 46 & ExactAgeinMonths < 47)
md46[['SampleID']] <- rownames(md46)
md46 <- plyr::ddply(md46, 'SubjectMaskId', head, n=1)
rownames(md46) <- md46[['SampleID']]

md47 <- subset(md, ExactAgeinMonths >= 47 & ExactAgeinMonths < 48)
md47[['SampleID']] <- rownames(md47)
md47 <- plyr::ddply(md47, 'SubjectMaskId', head, n=1)
rownames(md47) <- md47[['SampleID']]

md48 <- subset(md, ExactAgeinMonths >= 48 & ExactAgeinMonths < 49)
md48[['SampleID']] <- rownames(md48)
md48 <- plyr::ddply(md48, 'SubjectMaskId', head, n=1)
rownames(md48) <- md48[['SampleID']]


#Bind Comparisons
md2to3 <- rbind(md2, md3)
md3to4 <- rbind(md3, md4) 
md4to5 <- rbind(md4, md5) 
md5to6 <- rbind(md5, md6)
md6to7 <- rbind(md6, md7)
md7to8 <- rbind(md7, md8)
md8to9 <- rbind(md8, md9)
md9to10 <- rbind(md9, md10)
md10to11 <- rbind(md10, md11)
md11to12 <- rbind(md11, md12)
md12to13 <- rbind(md12, md13)
md13to14 <- rbind(md13, md14)
md14to15 <- rbind(md14, md15)
md15to16 <- rbind(md15, md16)
md16to17 <- rbind(md16, md17)
md17to18 <- rbind(md17, md18)
md18to19 <- rbind(md18, md19)
md19to20 <- rbind(md19, md20)
md20to21 <- rbind(md20, md21)
md21to22 <- rbind(md21, md22)
md22to23 <- rbind(md22, md23)
md23to24 <- rbind(md23, md24)
md24to25 <- rbind(md24, md25)
md25to26 <- rbind(md25, md26)
md26to27 <- rbind(md26, md27)
md27to28 <- rbind(md27, md28)
md28to29 <- rbind(md28, md29)
md29to30 <- rbind(md29, md30)
md30to31 <- rbind(md30, md31)
md31to32 <- rbind(md31, md32)
md32to33 <- rbind(md32, md33)
md33to34 <- rbind(md33, md34)
md34to35 <- rbind(md34, md35)
md35to36 <- rbind(md35, md36)
md36to37 <- rbind(md36, md37)
md37to38 <- rbind(md37, md38)
md38to39 <- rbind(md38, md39)
md39to40 <- rbind(md39, md40)
md40to41 <- rbind(md40, md41)
md41to42 <- rbind(md41, md42)
md42to43 <- rbind(md42, md43)
md43to44 <- rbind(md43, md44)
md44to45 <- rbind(md44, md45)
md45to46 <- rbind(md45, md46)
md46to47 <- rbind(md46, md47)
md47to48 <- rbind(md47, md48)

###RUN COMPARISONS###
# Set metadata to intended comparison
mf <- md2to3
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"2") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md2to3.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md3to4
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"3") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md3to4.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md4to5
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"4") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md4to5.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md5to6
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"5") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md5to6.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md6to7
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"6") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md6to7.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md7to8
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"7") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md7to8.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md8to9
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"8") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md8to9.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md9to10
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"9") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md9to10.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md10to11
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"10") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md10to11.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md11to12
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"11") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md11to12.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md12to13
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"12") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md12to13.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md13to14
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"13") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md13to14.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md14to15
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"14") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md14to15.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md15to16
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"15") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md15to16.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md16to17
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"16") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md16to17.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md17to18
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"17") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md17to18.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md18to19
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"18") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md18to19.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md19to20
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"19") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md19to20.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md20to21
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"20") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md20to21.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md21to22
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"21") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md21to22.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md22to23
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"22") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md22to23.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md23to24
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"23") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md23to24.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md24to25
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"24") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md24to25.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md25to26
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"25") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md25to26.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md26to27
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"26") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md26to27.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md27to28
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"27") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md27to28.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md28to29
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"28") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md28to29.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md29to30
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"29") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md29to30.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md30to31
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"30") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md30to31.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md31to32
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"31") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md31to32.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md32to33
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"32") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md32to33.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md33to34
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"33") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md33to34.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md34to35
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"34") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md34to35.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md35to36
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"35") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md35to36.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md36to37
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"36") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md36to37.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md37to38
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"37") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md37to38.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md38to39
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"38") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md38to39.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md39to40
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"39") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md39to40.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md40to41
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"40") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md40to41.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md41to42
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"41") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md41to42.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md42to43
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"42") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md42to43.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md43to44
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"43") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md43to44.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md44to45
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"44") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md44to45.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md45to46
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"45") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md45to46.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md46to47
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"46") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md46to47.txt", sep = "\t", row.names = TRUE, col.names = NA)

# Set metadata to intended comparison
mf <- md47to48
# Drop any samples not found in the distance matrix
dm <- dmfull
mf <- mf[rownames(mf) %in% rownames(dm),,drop=FALSE]
dm <- dm[rownames(mf), rownames(mf)]
# Convert from wide to long format, lookup sample groups, and categorize pairs
df          <- melt(dm, varnames=c('Sample1', 'Sample2'), value.name='Distance', as.is=TRUE)
df          <- df[which(df$Sample1 > df$Sample2),]
df$Group1   <- mf[as.character(df$Sample1), GroupingCol]
df$Group2   <- mf[as.character(df$Sample2), GroupingCol]
df$Category <- ifelse(df$Group1 == df$Group2, 'Within', 'Between')
# Label values with time point before saving
dffinal <- data.frame(df[,1:6],"47") 
write.table(dffinal, file = "/Users/CMMR/Desktop/16SBCdissimilarities/md47to48.txt", sep = "\t", row.names = TRUE, col.names = NA)

####################################################################################################################################################
#cat all files together into one big file = all.txt
#Get rid of quotes by: cat all.txt | sed -e 's/"//g' > allnoquotes.txt 
#Get rid of titles by: cat allnoquotes.txt | grep -v "Distance" > allready.txt
#Manually add titles back to top line

#Calculating the mean/median/standard deviation of categories so can draw line
BC <- read.table(file="/Users/CMMR/Desktop/16SBCdissimilarities/allready.txt", header = TRUE, sep = "\t")  
BC$NewCol <- paste(BC$Category, BC$Time, sep='') 
SixteenSWithinBetweenMeansMediansSDSE <- ddply(BC, .(NewCol), summarize, "Mean"= mean(Distance), "Median" = median(Distance), "SD" = sd(Distance), "SE" = sd(Distance)/sqrt(length(Distance)))
write.table(SixteenSWithinBetweenMeansMediansSDSE, "/Users/CMMR/Desktop/16SBCdissimilarities/Stats.txt")
#Manually add back the Time and Category

#Plotting medians & means over time
Stats <- read.table(file="/Users/CMMR/Desktop/16SBCdissimilarities/StatsReady.txt", header = TRUE, sep = "\t")
ggplot(Stats, aes(x=Time, y=Median, color=Category))  +
  theme(axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
  labs(x = "Months", y="Median of Bray-Curtis Dissimilarity") + 
  geom_point(size=1.5) + geom_smooth(method="loess", size=2, se=FALSE) +
  scale_color_manual(values=c("tan1", "steelblue")) 
ggplot(Stats, aes(x=Time, y=Mean, color=Category))  +
  theme(axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
  labs(x = "Months", y="Mean of Bray-Curtis Dissimilarity") + 
  geom_point(size=1.5) + geom_smooth(method="loess", span=0.35, size=2, se=FALSE) + geom_crossbar(aes(ymin=lower, ymax=upper), width=0.2) +
  scale_color_manual(values=c("tan1", "steelblue")) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48))
####################################################################################################################################################
#Plotting 16S & ITS2 Dissimilarity together
  BothStats <- read.table(file="/Users/CMMR/Desktop/BothBetaDivStatsReady.txt", header = TRUE, sep = "\t")
  ggplot(BothStats)  +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line    = element_line(color='black'),
          plot.margin=margin(l=20,b=15,t=15,r=20,unit = "pt")) + 
    labs(x = "Age (Months)") + 
    geom_smooth(aes(x=Time, y=Mean, color=Category), method="loess", size=2, se=FALSE) +
    geom_errorbar(aes(x=Time, ymin=lower, ymax=upper,color=Category), width=0.2) +
    scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48)) +
    coord_cartesian(xlim = c(0, 48), ylim = c(0.3, 1.0), expand = FALSE) +
    geom_point(aes(x=Time, y=Mean, shape = BetweenWithin, color=Category), size = 5) +
    scale_shape_manual(values = c(19, 1))+
    scale_color_manual(values = c("red","red","blue","blue")) +
    guides(col = guide_legend(reverse = TRUE))

#Export image as 2000 x 1000  
  
#legend.position="bottom",
#, y="Mean of Bray-Curtis Dissimilarity"
# axis.title.y = element_text(size = 30),
# axis.title.x = element_text(size = 30),
#legend.text = element_text(size=31),
#legend.key.size = unit(2, "cm"),
  
####################################################################################################################################################
#Plotting 16S & ITS2 Similarity together
  BothStats <- read.table(file="/Users/CMMR/Desktop/BothBetaDivStatsReady.txt", header = TRUE, sep = "\t")
  ggplot(BothStats)  +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line    = element_line(color='black'),
          plot.margin=margin(l=20,b=15,t=15,r=20,unit = "pt")) + 
    labs(x = "Age (Months)") + 
    geom_smooth(aes(x=Time, y=OneMinusMean, color=Category), method="loess", size=2, se=FALSE) +
    geom_errorbar(aes(x=Time, ymin=OneMinuslower, ymax=OneMinusupper,color=Category), width=0.2) +
    scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48)) +
    coord_cartesian(xlim = c(0, 48), ylim = c(0.0, 0.7), expand = FALSE) +
    geom_point(aes(x=Time, y=OneMinusMean, shape = BetweenWithin, color=Category), size = 5) +
    scale_shape_manual(values = c(19, 1))+
    scale_color_manual(values = c("red","red","blue","blue")) +
    guides(col = guide_legend(reverse = TRUE))
  
  #Export image as 2000 x 1000  
  