library("plyr")
library(phyloseq)
library(beepr)

biom.FilePath <- "/Users/tommyauchtung/Desktop/OTU_Table.Merged.99.justFungi.NoBleed.r3K.biom"
#OTU_Table.Merged.99.justFungi.NoBleed.r3K.biom generated in steps 1-5 of GeneratingFungiSpeciesTable.sh 
meta.FilePath <- "/Users/tommyauchtung/Desktop/HMPTEDDYmetadataUpdated.txt"

# Read in the OTU table
#----------------------
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
phylo2 <- phyloseq(otu_table(phylo), tax_table(phylo), sample_data(MappingFile))
ds <- transform_sample_counts(phylo2, function(x) x / sum(x))
long <- psmelt(tax_glom(ds, "Species"))
beep("wilhelm")

# Calculating p-values
# ----------------------
st <- ddply(long, 'Species', function (x) {
  kt <- kruskal.test(Abundance ~ Study, x)
  return(data.frame(pVal=kt$p.value))
})

# Exporting p-values
# ----------------------
pvals1  <- aggregate(pVal~Species, st, sum)
pvals2 <- pvals1[order(pvals1$pVal),]
write.table(pvals2, "/Users/tommyauchtung/Desktop/TEDDYHMPFungalSpeciesDifferencesPvalues.txt", sep="\t")

# Calculating and Exporting q-values
# ----------------------
st$qVal <- signif(p.adjust(st$pVal, "fdr"), 3)
qvals1  <- aggregate(qVal~Species, st, sum)
qvals2 <- qvals1[order(qvals1$qVal),]
write.table(qvals2, "/Users/tommyauchtung/Desktop/TEDDYHMPFungalSpeciesDifferences_FDR_Qvalues.txt", sep="\t")